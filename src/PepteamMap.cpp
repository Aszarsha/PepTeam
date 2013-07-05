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

#include "Matrices.hpp"
#include "Fasta.hpp"
#include "PepTree.hpp"

using namespace std;
using boost::range::for_each;

typedef string   fragment_t;
typedef uint32_t protein_t;
typedef string   peptide_t;

typedef map< peptide_t, vector< pair< size_t, double > > > pep2posAndScore_map;
typedef map< fragment_t, pep2posAndScore_map > frag2pep_map;
typedef map< protein_t, frag2pep_map > prot2pep_map;

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

	struct MappingResolvingFunctor {
			static size_t nbStringSimilarityResolved;
			static size_t nbMappingsResolved;

			double score;
			string const & query;
			string const & subject;
			uint32_t curProt;
			prot2pep_map & m;

			MappingResolvingFunctor( string const & q, string const & s, double sco, prot2pep_map & mm )
				: score( sco )
				, query( q )
				, subject( s )
				, m( mm ) {
						++nbStringSimilarityResolved;
			}

			void ListSize( uint16_t ) {   }

			void AddHeader( uint32_t protNumber, uint16_t ) {
					curProt = protNumber;
			}
			void StopHeader() {   }

			void AddPos( uint16_t p ) {
#if !defined( PROFILE_PERF )
					m[curProt][subject][query].push_back( make_pair( p, score ) );
#endif
					++nbMappingsResolved;
			}
			void StopPos() {   }
	};
	size_t MappingResolvingFunctor::nbMappingsResolved;
	size_t MappingResolvingFunctor::nbStringSimilarityResolved;

#if defined( PROFILE_PERF )
	static vector< tuple< size_t, size_t, size_t > > refuseStats;
	static vector< tuple< size_t, size_t, size_t > > acceptStats;
#endif

	inline size_t GetRangeNumLeaves( uint32_t start, uint32_t stop ) {
			return (stop - start) / LeavesLinkSize( fragSize );
	}

	void ResolveMapping( prot2pep_map & results
	                   , LinearizedTreeData const & query  , uint32_t queryStartIndex  , uint32_t queryStopIndex
	                   , LinearizedTreeData const & subject, uint32_t subjectStartIndex, uint32_t subjectStopIndex
	                   , double score
	                   ) {
			ForLeaf( query, queryStartIndex, [&,subjectStartIndex,score]( char const * qStr, uint32_t qOffset ) {
					ForLeaf( subject, subjectStartIndex, [&,score]( char const * sStr, uint32_t sOffset ) {
							string qString( qStr ), sString( sStr );
							ForLeafPos( subject, sOffset, MappingResolvingFunctor{ qString, sString, score, results } );
					});
			});
	}

	void MapTrees( prot2pep_map & results
	             , LinearizedTreeData const & query  , uint32_t queryIndex
	             , LinearizedTreeData const & subject, uint32_t subjectIndex
	             , SimilarityScore curScore, size_t depth
	             ) {
			ForNodeChildren( query, queryIndex
			               , [&,subjectIndex,depth,curScore]( size_t
			                                                , char queryChar, uint32_t queryChildIndex
			                                                , uint32_t queryStartLeaf, uint32_t queryStopLeaf
			                                                ) {
					ForNodeChildren( subject, subjectIndex
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
											double scoreVal = GetScoreNum( newScore ) / (double)GetScoreDen( newScore );
											ResolveMapping( results
											              , query  , queryStartLeaf  , queryStopLeaf
											              , subject, subjectStartLeaf, subjectStopLeaf
											              , scoreVal
											              );
#if defined( PROFILE_PERF )
											get< 0 >( acceptStats[depth-1] ) += 1;
											get< 1 >( acceptStats[depth-1] ) += GetRangeNumLeaves( queryStartLeaf  , queryStopLeaf   );
											get< 2 >( acceptStats[depth-1] ) += GetRangeNumLeaves( subjectStartLeaf, subjectStopLeaf );
#endif
							} else if ( depth < fragSize ) {
									MapTrees( results
									        , query, queryChildIndex, subject, subjectChildIndex
									        , newScore, depth + 1
									        );
							}
					});
			});
	}

}

prot2pep_map MapTrees( LinearizedTreeData const & query, LinearizedTreeData const & subject ) {
		prot2pep_map results;

		fragSize = GetTreeDepth( query );
		size_t d = GetTreeDepth( subject );
		if ( d != fragSize ) {
				ostringstream s;
				s << "Unable to map query over subject, different fragments sizes (" << d << " vs. " << fragSize << ')';
				throw std::runtime_error{ s.str() };
		}

#if defined( PROFILE_PERF )
		refuseStats.resize( fragSize );
		acceptStats.resize( fragSize );
#endif

		MapTrees( results, query, 0, subject, 0, { 0.0, 0.0 }, 1 );

		return results;
}

void PrintProteins2Peptides( FILE * file, prot2pep_map const & mappings ) {
		for_each( mappings, [&]( prot2pep_map::value_type const & prot ) {
				fprintf( file, "%u\t", prot.first );
				for_each( prot.second, [&]( frag2pep_map::value_type const & frag ) {
						for_each( frag.second, [&]( pep2posAndScore_map::value_type const & pep ) {
								for_each( pep.second, [&]( pair< size_t, double > const & pos ) {
										fprintf( file, "%s/%s[%zu](%g)\t"
										       , pep.first.c_str(), frag.first.c_str(), pos.first, pos.second
										       );
								});
						});
				});
				fprintf( file, "\n" );
		});
}

void UsageError( char * argv[] ) {
		fprintf( stderr, "Usage: %s pepTree-query-file pepTree-subject-file cutoff-homology\n", argv[0] );
		exit( 1 );
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

int main( int argc, char * argv[] ) {
		if ( argc != 4 ) {
				UsageError( argv );
		}

		cutoffHomology = atof( argv[3] );

		printf( "Similarity threshold: %f\n", cutoffHomology );

		FILE * subjectTreeFile = fopen( argv[2], "rb" );
		if ( !subjectTreeFile ) {
				fprintf( stderr, "Unable to open subject pepTree file \"%s\"\n", argv[2] );
				UsageError( argv );   // exit here to avoid creation of the other files if input file is invalid
		}
		string queryTreeFilename( argv[1] );
		FILE * queryTreeFile = fopen( argv[1], "rb" );
		if ( !queryTreeFile ) {
				fprintf( stderr, "Unable to open query pepTree file \"%s\"\n", argv[1] );
				UsageError( argv );   // exit here to avoid creation of the other files if input file is invalid
		}
		ostringstream prot2pepStream;
		prot2pepStream << queryTreeFilename << ".mapping.prot2pep." << (int)(100*cutoffHomology);
		FILE * prot2pepFile = fopen( prot2pepStream.str().c_str(), "w" );
		if ( !prot2pepFile ) {
				fprintf( stderr, "Unable to open output file \"%s\"\n", prot2pepStream.str().c_str() );
				UsageError( argv );   // exit here to avoid creation of the other files if input file is invalid
		}

		auto query   = ReadLinearizedTree( queryTreeFile );
		fclose( queryTreeFile );
		auto subject = ReadLinearizedTree( subjectTreeFile );
		fclose( subjectTreeFile );

		InitHomology();

		try {
				printf( "Intersecting peptides and proteins fragments sets...\n" );
				auto startTimer = chrono::high_resolution_clock::now();

				prot2pep_map mappings = MapTrees( GetTreeData( query ), GetTreeData( subject ) );

				auto finishTimer = chrono::high_resolution_clock::now();
				auto elapsed2 = finishTimer - startTimer;
				printf( "   ...%zu mappings found in %zu proteins in %ld seconds.\n"
				      , MappingResolvingFunctor::nbMappingsResolved, mappings.size()
				      , chrono::duration_cast< chrono::seconds >( elapsed2 ).count()
				      );

				cerr << "Printing proteins to peptides..." << endl;
				startTimer = chrono::high_resolution_clock::now();

				PrintProteins2Peptides( prot2pepFile, mappings );

				finishTimer = chrono::high_resolution_clock::now();
				auto elapsed4 = finishTimer - startTimer;
				cerr << "   ...done in " << chrono::duration_cast< chrono::seconds >( elapsed4 ).count() << " seconds." << endl;

				cerr << "Total execution time: "
				     << chrono::duration_cast< chrono::seconds >( elapsed2 + elapsed4 ).count()
				     << " seconds." << endl;
#if defined( PROFILE_PERF )
				PrintExecutionStats( stdout );
#endif
		} catch( std::exception & e ) {
				fprintf( stderr, "%s\n", e.what() );
				return 1;
		}
		return 0;
}
