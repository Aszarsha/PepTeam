#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <vector>
#include <set>
#include <sstream>
#include <utility>
#include <algorithm>
#include <cstring>
#include <chrono>
#include <climits>
#include <boost/range/algorithm/for_each.hpp>

#include "ProgramOptions.hpp"
#include "File.hpp"
#include "Matrices.hpp"
#include "Fasta.hpp"
#include "PepTree.hpp"

using namespace std;
using boost::range::for_each;

#define UsageFunction( function )                                                                   \
	do {                                                                                             \
			function( "Usage: ", argv[0], " pepTree-query-file pepTree-subject-file cutoff-homology"   \
			        , baseOption                                                                       \
			        );                                                                                 \
	} while ( false )

static size_t fragSize;
static double cutoffHomology;
static int homologyMatrix[24][24];
static int maxHomology = INT_MIN;
static int minHomology = INT_MAX;

namespace {

	inline void InitHomology() {
			memcpy( homologyMatrix, Matrix::Pam30, sizeof(homologyMatrix) );
			maxHomology = *max_element( &homologyMatrix[0][0] + 0, &homologyMatrix[23][23] + 1 );
			minHomology = *min_element( &homologyMatrix[0][0] + 0, &homologyMatrix[23][23] + 1 );
			printf( "MaxHomology: %d, MinHomology: %d\n", maxHomology, minHomology );
	}

	typedef pair< int, int > SimilarityScore;
	inline int GetScoreNum( SimilarityScore const & s ) {   return get< 0 >( s );   }
	inline int GetScoreDen( SimilarityScore const & s ) {   return get< 1 >( s );   }

	inline double WordsSimilarityFunction( char const * q, char const * s ) {
			int subjectCost = 0, queryCost = 0, homologyCost = 0;
			while ( *q != '\0' ) {
					Fasta::AAIndex qIdx = Fasta::Char2Index( *q );
					Fasta::AAIndex sIdx = Fasta::Char2Index( *s );

					queryCost    += homologyMatrix[qIdx][qIdx];
					subjectCost  += homologyMatrix[sIdx][sIdx];
					homologyCost += homologyMatrix[qIdx][sIdx];

					++q;
					++s;
			}
			return 2.0*homologyCost/static_cast< double >(queryCost + subjectCost);
	}

	inline SimilarityScore SimilarityFunction( char qChar, char sChar, SimilarityScore const & s ) {
			int n = 2*homologyMatrix[Fasta::Char2Index( qChar )][Fasta::Char2Index( sChar )];
			int d = homologyMatrix[Fasta::Char2Index( qChar )][Fasta::Char2Index( qChar )]
			      + homologyMatrix[Fasta::Char2Index( sChar )][Fasta::Char2Index( sChar )];
			return { GetScoreNum( s ) + n, GetScoreDen( s ) + d };
	}

	bool Refuse( SimilarityScore const & s, size_t depth ) {
			double k = 2.0*((fragSize-depth)*(double)maxHomology);
			return (GetScoreNum( s ) + k) / (GetScoreDen( s ) + k) < cutoffHomology;   // early refuse
	}

	bool Accept( SimilarityScore const & s, size_t depth ) {
			double k = 2.0*((fragSize-depth)*(double)maxHomology);
			double l = 2.0*((fragSize-depth)*(double)minHomology);
			return (GetScoreNum( s ) + l) / (GetScoreDen( s ) + k) >= cutoffHomology;   // early accept
	}

	size_t nbStringSimilarity = 0;
#if defined( PROFILE_PERF )
	static vector< tuple< size_t, size_t, size_t > > refuseStats;
	static vector< tuple< size_t, size_t, size_t > > acceptStats;
#endif

	inline size_t GetRangeNumLeaves( uint32_t start, uint32_t stop ) {
			return (stop - start) / LeavesLinkSize( fragSize );
	}

	template< typename F >
	void ResolveMapping( FILE * file
	                   , MMappedPepTree const & query  , uint32_t queryStartIndex  , uint32_t queryStopIndex
	                   , MMappedPepTree const & subject, uint32_t subjectStartIndex, uint32_t subjectStopIndex
	                   , F && scoreFunc
	                   ) {
			for ( uint32_t qIdx = queryStartIndex; qIdx != queryStopIndex; ++qIdx ) {
					for ( uint32_t sIdx = subjectStartIndex; sIdx != subjectStopIndex; ++sIdx ) {
							query.ForLeaf( qIdx, [=]( char const * qStr, uint32_t ) {
									subject.ForLeaf( sIdx, [=]( char const * sStr, uint32_t ) {
											fprintf( file, "%u %u %g\n", qIdx, sIdx, scoreFunc( qStr, sStr ) );
									});
							});
							++nbStringSimilarity;
					}
			}
	}

	void MapTrees( FILE * file
	             , MMappedPepTree const & query  , uint32_t queryIndex
	             , MMappedPepTree const & subject, uint32_t subjectIndex
	             , SimilarityScore curScore, size_t depth
	             ) {
			query.ForNodeChildren( queryIndex
			               , [&,subjectIndex,depth,curScore]( size_t
			                                                , char queryChar, uint32_t queryChildIndex
			                                                , uint32_t queryStartLeaf, uint32_t queryStopLeaf
			                                                ) {
					subject.ForNodeChildren( subjectIndex
					               , [&,depth,curScore]( size_t
					                                   , char subjectChar, uint32_t subjectChildIndex
					                                   , uint32_t subjectStartLeaf, uint32_t subjectStopLeaf
					                                   ) {
							auto newScore = SimilarityFunction( queryChar, subjectChar, curScore );
							if ( Refuse( newScore, depth ) ) {
#if defined( PROFILE_PERF )
									get< 0 >( refuseStats[depth-1] ) += 1;
									get< 1 >( refuseStats[depth-1] ) += GetRangeNumLeaves( queryStartLeaf  , queryStopLeaf   );
									get< 2 >( refuseStats[depth-1] ) += GetRangeNumLeaves( subjectStartLeaf, subjectStopLeaf );
#endif
							} else if ( Accept( newScore, depth ) ) {
									if ( depth == fragSize ) {
											double scoreVal = GetScoreNum( newScore ) / (double)GetScoreDen( newScore );
											ResolveMapping( file
											              , query  , queryStartLeaf  , queryStopLeaf
											              , subject, subjectStartLeaf, subjectStopLeaf
											              , [scoreVal]( char const *, char const * ) {   return scoreVal;   }
											              );
									} else {
											ResolveMapping( file
											              , query  , queryStartLeaf  , queryStopLeaf
											              , subject, subjectStartLeaf, subjectStopLeaf
											              , &WordsSimilarityFunction
											              );
									}
#if defined( PROFILE_PERF )
									get< 0 >( acceptStats[depth-1] ) += 1;
									get< 1 >( acceptStats[depth-1] ) += GetRangeNumLeaves( queryStartLeaf  , queryStopLeaf   );
									get< 2 >( acceptStats[depth-1] ) += GetRangeNumLeaves( subjectStartLeaf, subjectStopLeaf );
#endif
							} else if ( depth < fragSize ) {
									MapTrees( file
									        , query, queryChildIndex, subject, subjectChildIndex
									        , newScore, depth + 1
									        );
							}
					});
			});
	}

}

void MapTrees( FILE * file, MMappedPepTree const & query, MMappedPepTree const & subject ) {
		fragSize = query.Depth();
		size_t d = subject.Depth();
		if ( d != fragSize ) {
				ostringstream s;
				s << "Unable to map query over subject, different fragments sizes (" << d << " vs. " << fragSize << ')';
				throw std::runtime_error{ s.str() };
		}

#if defined( PROFILE_PERF )
		refuseStats.resize( fragSize );
		acceptStats.resize( fragSize );
#endif

		MapTrees( file, query, 0, subject, 0, { 0.0, 0.0 }, 1 );
}

#if defined( PROFILE_PERF )
void PrintExecutionStats( FILE * file ) {
		fprintf( file, "Refuses:\n" );
		size_t height = 0;
		for_each( refuseStats, [&]( tuple< size_t, size_t, size_t > const & v ) {
				size_t size = get< 0 >( v );
				size_t totalQueryCut   = get< 1 >( v );
				size_t totalSubjectCut = get< 2 >( v );
				fprintf( file, "%zu: %10zu | %10zu / %10zu\n"
				       , height, get< 0 >( v ), totalQueryCut, totalSubjectCut
				       );
				++height;
		});
		fprintf( file, "Accepts:\n" );
		height = 0;
		for_each( acceptStats, [&]( tuple< size_t, size_t, size_t > const & v ) {
				size_t size = get< 0 >( v );
				size_t totalQueryCut   = get< 1 >( v );
				size_t totalSubjectCut = get< 2 >( v );
				fprintf( file, "%zu: %10zu | %10zu / %10zu\n"
				       , height, get< 0 >( v ), totalQueryCut, totalSubjectCut
				       );
				++height;
		});
}
#endif

int main( int argc, char * argv[] ) try {
		po::options_description baseOption( "Options", 999, 999 );
		baseOption.add_options()
			( "help,h"  , "Print this help message and exit" )
			( "output,o", po::value< string >(), "Output file (default: input-query-file+\".mapping.\"+cutoff-homology)" )
			;

		po::options_description hiddenOptions( "Required options", 999, 999 );
		hiddenOptions.add_options()
			( "query-file"  , po::value< string >()->required(), "Input .pepTree query file" )
			( "subject-file", po::value< string >()->required(), "Input .pepTree subject file" )
			( "cutoff-homology", po::value< double >()->required(), "Target minimal homology threshold" )
			;
		po::positional_options_description posOptions;
		posOptions.add( "query-file"  , 1 );
		posOptions.add( "subject-file", 2 );
		posOptions.add( "cutoff-homology", 3 );

		po::options_description cmdLineOptions( "Command line options", 999, 999 );
		cmdLineOptions.add( baseOption ).add( hiddenOptions );

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

		cutoffHomology = vm["cutoff-homology"].as< double >();
		printf( "Similarity threshold: %f\n", cutoffHomology );


		string outFileName = vm.count( "output" ) ? vm["output"].as< string >()
		                                          : vm["query-file"].as< string >() + ".mapping."
		                                               + to_string( (int)cutoffHomology ) + "_"
		                                               + to_string( (int)floor( 100*(cutoffHomology - (int)cutoffHomology) ) );
		auto outFile = OpenFile( outFileName.c_str(), "wb" );

		MMappedPepTree query( vm["query-file"].as< string >().c_str() );
		MMappedPepTree subject( vm["subject-file"].as< string >().c_str() );

		InitHomology();

		try {
				printf( "Intersecting peptides and proteins fragments sets...\n" );
				auto startTimer = chrono::high_resolution_clock::now();

				MapTrees( outFile, query, subject );

				auto finishTimer = chrono::high_resolution_clock::now();
				auto elapsed2 = finishTimer - startTimer;
				printf( "   ...%zu mappings found in %ld seconds.\n"
				      , nbStringSimilarity, chrono::duration_cast< chrono::seconds >( elapsed2 ).count()
				      );
#if defined( PROFILE_PERF )
				PrintExecutionStats( stdout );
#endif
		} catch( std::exception & e ) {
				fprintf( stderr, "%s\n", e.what() );
				return 1;
		}
		return 0;
} catch ( std::exception & e ) {
		fprintf( stderr, "Unhandled exception reached main: %s, application will now exit\n", e.what() );
		return 271828;
} catch ( ... ) {
		fprintf( stderr, "Unhandled object reached main, application will now exit\n" );
		return 314159;
}
