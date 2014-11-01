#include <cassert>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cstring>
#include <chrono>
#include <limits>
#include <cstdint>

#include "ProgramOptions.hpp"
#include "File.hpp"
#include "Fasta.hpp"
#include "PepTree.hpp"
#include "FastIdx.hpp"

using namespace std;

#define UsageFunction( function )                                                                          \
	do {                                                                                                    \
			function( "Usage: ", argv[0], " [options] -c fastIdx-file fragments-size\n"                       \
			                "  -c [ --create ]       create the pepTree file from the input fastIdx file\n"   \
			        , "or   : ", argv[0], " [options] -* pepTree-file\n"                                      \
			        , requiredOptions, baseOption                                                             \
			        );                                                                                        \
	} while ( false )

void PepTreeCreation( MMappedFastIdx const & idx, size_t fragSize, FILE * outFile ) {
		printf( "Creating PepTree of depth %zu...\n", fragSize );
		auto startTimer = chrono::high_resolution_clock::now();
		auto indices   = idx.GetIndicesData();
		auto names     = idx.GetNamesData();
		auto sequences = idx.GetSequencesData();
		auto seqSize   = idx.GetSequencesSize();

		Trie trie( fragSize );

		auto proteinIndex  = size_t{ 0 };
		auto fragmentStart = size_t{ 0 }, sequenceStart = size_t{ 0 };
		for ( size_t i = 0; i < seqSize; ++i ) {
				while ( !Fasta::IsValidAA( sequences[ i ] ) ) {   // require a loop to detect multiple successive invalid chars
						if ( sequences[ i ] == '\0' ) {
								++proteinIndex;
								sequenceStart = i+1;
						} else {
								auto startStr = sequences + (i < fragSize ? 0 : i-fragSize+1);
								fprintf( stderr
								       , "      warning: ignoring fragments "
								         "in %s (centered at position %zu) "
								         "in \"%s\" [%zu]: fragments contain "
								       , string( startStr, startStr + 2*fragSize-1 ).c_str()
								       , i - sequenceStart, names + indices[ proteinIndex*2 ], proteinIndex
								       );
								switch ( sequences[ i ] ) {
									case 'U': {   fprintf( stderr, "selenocysteine (U)\n" );                       } break;
									case '-': {   fprintf( stderr, "gap of indeterminate length (-)\n" );          } break;
									default:  {   fprintf( stderr, "invalid character (%c)\n", sequences[ i ] );   } break;
								}
						}

						i = fragmentStart = i+1;
						if ( i >= seqSize ) {   // last character was invalid
								goto End_Outer_Loop;
						}
				}
				if ( i - fragmentStart == fragSize-1 ) {
						auto leafIndex = trie.GetLeafCreatePath( sequences + fragmentStart );
						trie.GetLeaf( leafIndex )->positions[ proteinIndex ].push_back( fragmentStart - sequenceStart );

						++fragmentStart;
				}
		}
	End_Outer_Loop:
		auto finishTimer = chrono::high_resolution_clock::now();
		auto elapsed1 = finishTimer - startTimer;
		printf( "   ...PepTree trie created in %ld seconds.\n"
		      , chrono::duration_cast< chrono::seconds >( elapsed1 ).count()
		      );

		try {
				auto nbNodes = trie.NumLeaves() + trie.NumNodes();
				printf( "Linearizing tree structure (%zu node%s: %zu internal%s, %zu lea%s)...\n"
				      , nbNodes, nbNodes > 1 ? "s" : ""
				      , trie.NumNodes(), trie.NumNodes() > 1 ? "s" : ""
				      , trie.NumLeaves(), trie.NumLeaves() > 1 ? "ves" : "f"
				      );
				startTimer = chrono::high_resolution_clock::now();

				auto tree = trie.LinearizeTree();

				finishTimer = chrono::high_resolution_clock::now();
				auto elapsed2 = finishTimer - startTimer;
				printf( "   ...linearized in %ld seconds .\n"
				      , chrono::duration_cast< chrono::seconds >( elapsed2 ).count()
				      );

				tree.Write( outFile );
		} catch( std::exception & e ) {
				fprintf( stderr, "%s\n", e.what() );
				exit( 1 );
		}
}

int main( int argc, char * argv[] ) try {
		po::options_description baseOption( "Options", 999, 999 );
		baseOption.add_options()
			( "help,h"  , "Print this help message and exit" )
			( "output,o", po::value< string >(), "Output file (default: -c -> input-file+\".pepTree.\"+fragments-size\n"
			                                     "                    else -> stdout)" )
			;

		po::options_description requiredOptions( " where -* is one of", 999, 999 );
		requiredOptions.add_options()
		   ( "depth,d"  , "print the tree depth" )
		   ( "vector,v" , "print the tree vector in human 'interpretable' format" )
		   ( "nodes,n"  , "print the tree nodes in human 'interpretable' format" )
		   ( "leaves,l" , "print the tree leaves in human 'interpretable' format" )
		   ( "offsets,p", "print the tree leaf offsets in human 'interpretable' format" )
			;
		po::options_description hiddenRequiredOptions( " -c", 999, 999 );
		hiddenRequiredOptions.add_options()
			( "create,c", "create the pepTree file from the input fastIdx file" )
			;

		po::options_description hiddenOptions( "Required options", 999, 999 );
		hiddenOptions.add_options()
			( "input-file"    , po::value< string >()->required(), "Input .fasta file" )
			( "fragments-size", po::value< size_t >(), "Peptide size in number of AA" )
			;
		po::positional_options_description posOptions;
		posOptions.add( "input-file"    , 1 );
		posOptions.add( "fragments-size", 1 );

		po::options_description cmdLineOptions( "Command line options", 999, 999 );
		cmdLineOptions.add( baseOption ).add( requiredOptions ).add( hiddenRequiredOptions ).add( hiddenOptions );

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
		if ( !DetectExclusiveOptions( vm, "create", "depth", "vector", "nodes", "leaves", "offsets" ) ) {
				cerr << "Mutually exclusive options error" << endl;
				UsageFunction( UsageError );
		}

		if ( vm.count( "create" ) ) {
				if ( !vm.count( "fragments-size" ) ) {
						cerr << "Missing fragments-size" << endl;
						UsageFunction( UsageError );
				}
				MMappedFastIdx idx( vm["input-file"].as< string >().c_str() );
				size_t fragSize = vm["fragments-size"].as< size_t >();
				string outFileName = vm.count( "output" ) ? vm["output"].as< string >()
				                                          : vm["input-file"].as< string >()+".pepTree."+to_string(fragSize);
				auto outFile = OpenFile( outFileName.c_str(), "wb" );
				PepTreeCreation( idx, fragSize , outFile );
		} else {
				MMappedPepTree tree( vm["input-file"].as< string >().c_str() );
				auto outFile = vm.count( "output" ) ? OpenFile( vm["output"].as< string >().c_str(), "wb" )
				                                    : AutoFILE( stdout, []( auto ) {   return 0;   } );
				
				if ( vm.count( "depth" ) ) {
   					fprintf( outFile, "Depth: %u\n", tree.Depth() );
				} else if ( vm.count( "vector" ) ) {
						tree.WriteReadableTree( outFile );
				} else if ( vm.count( "nodes" ) ) {
						tree.WriteReadableNodes( outFile );
				} else if ( vm.count( "leaves" ) ) {
						tree.WriteReadableLeaves( outFile );
				} else { //if ( vm.count( "offsets" ) {
						tree.WriteReadableLeafPos( outFile );
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
