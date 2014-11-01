#include <cstdio>
#include <cassert>
#include <string>
#include <sstream>
#include <vector>
#include <chrono>
#include <cstring>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <boost/program_options.hpp>

#include "ProgramOptions.hpp"
#include "File.hpp"
#include "FastIdx.hpp"

namespace po = boost::program_options;
using namespace std;

#define UsageFunction( function )                                                                     \
	do {                                                                                               \
			function( "Usage: ", argv[0], " [options] -* input-file\n", requiredOptions, baseOption );   \
	} while ( false )

void IndexCreation( FILE * inFile, FILE * outFile ) {
		auto proteinsName = vector< string >{};
		auto proteinsSeq  = vector< string >{};
		auto ProcessSequence = [ & ]( string && name, string && seq ) {
				proteinsName.push_back( move( name ) );
				proteinsSeq.push_back( move( seq ) );
		};

		printf( "Indexing input file...\n" );

		enum class ReadingState { Start, NewLineFromSeq, Header, Sequence };
		ReadingState readingState = ReadingState::Start;
		size_t lineNumber = 0, nbSeqProcessed = 0;
		string text, name;
		auto startTimer = chrono::high_resolution_clock::now();
		for ( char c; (c = fgetc_unlocked( inFile )) != EOF; ) {
				if ( readingState == ReadingState::Start ) {
						if ( c != '>' ) {
								fprintf( stderr, "Invalid '%c' character at start of file, was expecting '>'\n", c );
						} else {
								readingState = ReadingState::Header;
						}
				} else if ( c == '>' && readingState == ReadingState::NewLineFromSeq ) {
						readingState = ReadingState::Header;

						ProcessSequence( move( name ), move( text ) );
						++nbSeqProcessed;
						text = string{};
				} else if ( c == '\n' ) {
						++lineNumber;
						switch ( readingState ) {
							case ReadingState::NewLineFromSeq: {
									fprintf( stderr, "      warning: suspicious new line at line number %zu; continuing\n", lineNumber + 1 );
								} break;
	
							case ReadingState::Header: {
									readingState = ReadingState::Sequence;
									name = move( text );
									text = string{};
								} break;
	
							case ReadingState::Sequence: {
									readingState = ReadingState::NewLineFromSeq;
								} break;
	
							default: {
								} assert( !"Should never get here!" );
						}
				} else {
						if ( readingState == ReadingState::NewLineFromSeq ) {
								readingState = ReadingState::Sequence;
						}
						text += c;
				}
		}
		if ( readingState != ReadingState::Header && !text.empty() ) {
				ProcessSequence( move( name ), move( text ) );
				++nbSeqProcessed;
		}
		auto finishTimer = chrono::high_resolution_clock::now();
		auto elapsed1 = finishTimer - startTimer;
		printf( "   ...indexed %zu sequence%s from input file in %ld seconds.\n"
		      , nbSeqProcessed, nbSeqProcessed > 1 ? "s" : ""
		      , chrono::duration_cast< chrono::seconds >( elapsed1 ).count()
		      );

		try {
				printf( "Writting protein index structure...\n" );
				startTimer = chrono::high_resolution_clock::now();

				MemFastIdx idx( proteinsName, proteinsSeq );
				idx.Write( outFile );

				finishTimer = chrono::high_resolution_clock::now();
				auto elapsed2 = finishTimer - startTimer;
				printf( "   ...written in %ld seconds.\n"
				      , chrono::duration_cast< chrono::seconds >( elapsed2 ).count()
				      );
		} catch( std::exception & e ) {
				fprintf( stderr, "%s\n", e.what() );
				exit( 1 );
		}
}

void IndexSizePrinting( MMappedFastIdx const & idx, FILE * outFile ) {
		fprintf( outFile, "Number of proteins in index: %zu\n", idx.Size() );
}

void IndexPrinting( MMappedFastIdx const & idx, FILE * outFile ) {
		auto size      = idx.GetIndicesSize();
		auto indices   = idx.GetIndicesData();
		auto names     = idx.GetNamesData();
		auto sequences = idx.GetSequencesData();

		auto nodeNumber = size_t{ 0 };
		for ( size_t i = 0; i < size; i += 2 ) {
				fprintf( outFile, "(%05zu) %06zX: %06X | %06X\n", nodeNumber++, i, indices[ i ], indices[ i+1 ] );
		}

		fprintf( outFile, "~~~~~~\n" );

		nodeNumber = 0;
		for ( size_t i = 0; i < size; i += 2 ) {
				fprintf( outFile, "(%05zu) %06X: %s\n", nodeNumber++, indices[ i ], &names[ indices[ i ] ] );
		}

		fprintf( outFile, "~~~~~~\n" );

		nodeNumber = 0;
		for ( size_t i = 1; i < size; i += 2 ) {
				fprintf( outFile, "(%05zu) %06X: %s\n", nodeNumber++, indices[ i ], &sequences[ indices[ i ] ] );
		}
}

int main( int argc, char * argv[] ) try {
		po::options_description baseOption( "Options", 999, 999 );
		baseOption.add_options()
			( "help,h"  , "Print this help message and exit" )
			( "output,o", po::value< string >(), "Output file (default: -c -> input-file+\".fastIdx\"\n"
			                                     "                    else -> stdout)" )
			;

		po::options_description requiredOptions( " where -* is one of", 999, 999 );
		requiredOptions.add_options()
			( "create,c", "create the protein index from input fasta file" )
			( "size,s"  , "print the number of proteins in the index" )
			( "print,p" , "print the index in human 'interpretable' format" )
			;

		po::options_description hiddenOptions( "Required options", 999, 999 );
		hiddenOptions.add_options()
			( "input-file", po::value< string >()->required(), "Input .fasta file" )
			;
		po::positional_options_description posOptions;
		posOptions.add( "input-file", 1 );

		po::options_description cmdLineOptions( "Command line options", 999, 999 );
		cmdLineOptions.add( baseOption ).add( requiredOptions ).add( hiddenOptions );

		po::variables_map vm;
		try {
				po::store( po::command_line_parser( argc, argv ).options( cmdLineOptions )
				                                                .positional( posOptions )
				                                                .run(), vm );
				po::notify( vm );
		} catch ( po::error & e ) {
				cerr << "Invalid command line: " << e.what() << endl;
				UsageFunction( UsageError );
		}
		if ( vm.count( "help" ) ) {
				UsageFunction( Usage );
		}
		if ( !DetectExclusiveOptions( vm, "create", "size", "print" ) ) {
				cerr << "Mutually exclusive options error" << endl;
				UsageFunction( UsageError );
		}

		if ( vm.count( "create" ) ) {
				auto inFile = OpenFile( vm["input-file"].as< string >().c_str(), "rb" );
				string outFileName = vm.count( "output" ) ? vm["output"].as< string >()
				                                          : vm["input-file"].as< string >()+".fastIdx";
				auto outFile = OpenFile( outFileName.c_str(), "wb" );
				IndexCreation( inFile, outFile );
		} else {
				MMappedFastIdx idx( vm["input-file"].as< string >().c_str() );
				auto outFile = vm.count( "output" ) ? OpenFile( vm["output"].as< string >().c_str(), "wb" )
				                                    : AutoFILE( stdout, []( auto ) {   return 0;   } );
				if ( vm.count( "size" ) ) {
						IndexSizePrinting( idx, outFile );
				} else { //if ( vm.count( "print" ) ) {
						IndexPrinting( idx, outFile );
				}
		}

		return 0;
} catch ( std::exception & e ) {
		fprintf( stderr, "Unhandled exception reached main: %s, application will now exit\n", e.what() );
		return 271828;
} catch ( ... ) {
		fprintf( stderr, "Unhandled object reached main, application will now exit\n" );
		return 314159;
}
