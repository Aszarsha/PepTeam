#include <cstdio>
#include <cassert>
#include <string>
#include <sstream>
#include <vector>
#include <chrono>
#include <cstring>
#include <cstdint>
#include <stdexcept>

#include "FastIdx.hpp"

using namespace std;

void UsageError( char * argv[] ) {
		fprintf( stderr
		       , "Usage: %s -* input-file\n"
		         "   where * is one of:\n"
		         "     c -> create the protein index from input fasta file\n"
		         "     s -> print the number of proteins in the index\n"
		         "     p -> print the index in human 'interpretable' format\n"
		       , argv[ 0 ]
		       );
		exit( 1 );
}

FILE * OpenInputFile( char const * filename ) {
		FILE * inputFile = fopen( filename, "rb" );
		if ( !inputFile ) {
				fprintf( stderr, "Unable to open input file \"%s\"\n", filename );
				exit( 1 );
		}
		return inputFile;
}

void IndexCreation( char * argv[] ) {
		char const * filename = argv[ 2 ];
		FILE * inputFile = OpenInputFile( filename );

		auto proteinsName = vector< string >{};
		auto proteinsSeq  = vector< string >{};
		auto ProcessSequence = [ & ]( string && name, string && seq ) {
				proteinsName.push_back( move( name ) );
				proteinsSeq.push_back( move( seq ) );
		};

		printf( "Indexing \"%s\"...\n", filename );

		enum class ReadingState { Start, NewLineFromSeq, Header, Sequence };
		ReadingState readingState = ReadingState::Start;
		size_t lineNumber = 0, nbSeqProcessed = 0;
		string text, name;
		auto startTimer = chrono::high_resolution_clock::now();
		for ( char c; (c = fgetc_unlocked( inputFile )) != EOF; ) {
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
		      , nbSeqProcessed, nbSeqProcessed > 1 ? "s" : "", chrono::duration_cast< chrono::seconds >( elapsed1 ).count()
		      );
		fclose( inputFile );

		try {
				printf( "Writting protein index structure...\n" );
				startTimer = chrono::high_resolution_clock::now();

				stringstream outputProtIdxFilenameStream;
				outputProtIdxFilenameStream << filename << ".fastIdx";
				FILE * outputProtIdxFile = fopen( outputProtIdxFilenameStream.str().c_str(), "wb" );
				if ( !outputProtIdxFile ) {
						fprintf( stderr, "Unable to open output file \"%s\"\n", outputProtIdxFilenameStream.str().c_str() );
						exit( 1 );
				}
				WriteLinearizedProtIdx( outputProtIdxFile, LinearizeProtIdx( proteinsName, proteinsSeq ) );
				fclose( outputProtIdxFile );

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

void IndexSizePrinting( char * argv[] ) {
		auto inputFile = OpenInputFile( argv[ 2 ] );
		printf( "Number of proteins in index: %zu\n", ReadProteinIndexSize( inputFile ) );
		fclose( inputFile );
}

void IndexPrinting( char * argv[] ) {
		auto inputFile = OpenInputFile( argv[ 2 ] );
		auto idx = ReadProteinIndex( inputFile );
		fclose( inputFile );

		auto const & indices   = GetProtIdxIndices  ( idx );
		auto const & names     = GetProtIdxNames    ( idx );
		auto const & sequences = GetProtIdxSequences( idx );

		auto nodeNumber = size_t{ 0 };
		for ( size_t i = 0, e = indices.size(); i < e; i += 2 ) {
				printf( "(%05zu) %06zX: %06X | %06X\n", nodeNumber++, i, indices[ i ], indices[ i+1 ] );
		}

		printf( "~~~~~~\n" );

		nodeNumber = 0;
		for ( size_t i = 0, e = indices.size(); i < e; i += 2 ) {
				printf( "(%05zu) %06X: %s\n", nodeNumber++, indices[ i ], &names[ indices[ i ] ] );
		}

		printf( "~~~~~~\n" );

		nodeNumber = 0;
		for ( size_t i = 1, e = indices.size(); i < e; i += 2 ) {
				printf( "(%05zu) %06X: %s\n", nodeNumber++, indices[ i ], &sequences[ indices[ i ] ] );
		}
}

int main( int argc, char * argv[] ) {
		if ( argc != 3 || argv[ 1 ][ 0 ] != '-' || strlen( argv[ 1 ] ) != 2 ) {
				UsageError( argv );
		}

		switch ( argv[ 1 ][ 1 ] ) {
			case 'c': {   IndexCreation    ( argv );   } break;
			case 's': {   IndexSizePrinting( argv );   } break;
			case 'p': {   IndexPrinting    ( argv );   } break;
			default: UsageError( argv );
		}

		return 0;
}
