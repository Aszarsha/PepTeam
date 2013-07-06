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

#include "Matrices.hpp"
#include "Fasta.hpp"
#include "PepTree.hpp"
#include "FastIdx.hpp"

using namespace std;

void UsageError( char * argv[] ) {
		fprintf( stderr
		       , "Usage: %s -c input-FastIdx-file fragments-size\n"
		         "          create the PepTree file from the input FastIdx file\n"
		         "   or: %s -%% pepTree-file\n"
		         "   where %% is one of:\n"
		         "     d -> print the tree depth\n"
		         "     v -> print the tree vector in human 'interpretable' format\n"
		         "     n -> print the tree nodes in human 'interpretable' format\n"
		         "     l -> print the tree leaves in human 'interpretable' format\n"
		         "     p -> print the tree leaf positions in human 'interpretable' format\n"
		       , argv[ 0 ], argv[ 0 ]
		       );
		exit( 1 );
}

void PepTreeCreation( char * argv[] ) {
		auto fragSize = static_cast< size_t >( atoi( argv[ 3 ] ) );
		MMappedFastIdx idx( argv[2] );

		printf( "Creating PepTree of depth %zu from \"%s\"...\n", fragSize, argv[ 2 ] );
		auto indices   = idx.GetIndices();
		auto names     = idx.GetNames();
		auto sequences = idx.GetSequences();
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
						trie.GetLeaf( leafIndex )->positions[ proteinIndex ].push_back( i - sequenceStart );

						++fragmentStart;
				}
		}
	End_Outer_Loop:

		try {
				printf( "Linearizing tree structure...\n" );
				auto startTimer = chrono::high_resolution_clock::now();

				auto tree = trie.LinearizeTree();

				auto finishTimer = chrono::high_resolution_clock::now();
				auto elapsed2 = finishTimer - startTimer;
				auto nbNodes = trie.NumLeaves() + trie.NumNodes();
				printf( "   ...linearized in %ld seconds (%zu node%s: %zu internal%s, %zu lea%s).\n"
				      , chrono::duration_cast< chrono::seconds >( elapsed2 ).count()
				      , nbNodes, nbNodes > 1 ? "s" : ""
				      , trie.NumNodes(), trie.NumNodes() > 1 ? "s" : ""
				      , trie.NumLeaves(), trie.NumLeaves() > 1 ? "ves" : "f"
				      );

				ostringstream outputPepTreeFilenameStream;
				outputPepTreeFilenameStream << argv[ 2 ] << ".pepTree." << fragSize;
				auto outputPepTreeFile = fopen( outputPepTreeFilenameStream.str().c_str(), "wb" );
				if ( !outputPepTreeFile ) {
						fprintf( stderr, "Unable to open output file \"%s\"\n", outputPepTreeFilenameStream.str().c_str() );
						exit( 1 );
				}
				WriteLinearizedTree( outputPepTreeFile, GetTreeData( tree ) );
				fclose( outputPepTreeFile );
		} catch( std::exception & e ) {
				fprintf( stderr, "%s\n", e.what() );
				exit( 1 );
		}
}

int main( int argc, char * argv[] ) {
		if ( argc < 3 || argc > 4 || argv[ 1 ][ 0 ] != '-' || strlen( argv[ 1 ] ) != 2 ) {
				UsageError( argv );
		}

		if ( argv[ 1 ][ 1 ] == 'c' ) {
				if ( argc != 4 ) {
						UsageError( argv );
				}
				PepTreeCreation( argv );
		} else {
				if ( argc != 3 ) {
						UsageError( argv );
				}
				auto treeFile = fopen( argv[ 2 ], "rb" );
				if ( !treeFile ) {
						fprintf( stderr, "Unable to open PepTree file \"%s\"\n", argv[ 2 ] );
						exit( 1 );
				}
				auto tree = ReadLinearizedTree( treeFile );

				switch ( argv[ 1 ][ 1 ] ) {
					case 'd': {   fprintf( stdout, "Depth: %d\n", GetTreeDepth( tree ) );                       } break;
					case 'v': {   WriteReadableLinearizedTree   ( stdout, GetTreeData( tree ) );                } break;
					case 'n': {   WriteReadableLinearizedNodes  ( stdout, GetTreeData( tree ) );                } break;
					case 'l': {   WriteReadableLinearizedLeaves ( stdout, GetTreeData( tree ) );                } break;
					case 'p': {   WriteReadableLinearizedLeafPos( stdout, GetTreeData( tree ) );                } break;
					default: UsageError( argv );
				}
				fclose( treeFile );
		}

		return 0;
}
