#include <cstdio>
#include <cstdint>
#include <cstring>
#include <boost/range/algorithm/for_each.hpp>

#include "FastIdx.hpp"
#include "PepTree.hpp"

using namespace std;
using boost::range::for_each;

typedef map< uint32_t, vector< unsigned int > > protProfile_map;
protProfile_map protProfiles;

struct ProtFunctor {
	public:
		ProtFunctor( MMappedFastIdx const & idx_, size_t pepSz )
			: idx( idx_ )
			, pepSize( pepSz ) {
		}

	public:
		void ListSize( uint16_t ) {   }

		void AddHeader( uint32_t protNumber, uint16_t ) {
				auto curProt = protProfiles.find( protNumber );
				if ( curProt == protProfiles.end() ) {
						size_t seqLength = strlen( idx.GetSequence( protNumber ) );
						curProt = protProfiles.insert( make_pair( protNumber
						                                        , vector< unsigned int >( seqLength )
						                                        ) ).first;
				}
				curVec = &curProt->second;
		}
		void StopHeader() {   }

		void AddPos( uint16_t p ) {
				for ( size_t i = p, e = p + pepSize; i < e; ++i ) {
						++(*curVec)[i];
				}
		}
		void StopPos() {   }

	private:
		MMappedFastIdx const & idx;
		size_t pepSize;

		vector< unsigned int > * curVec;
};

int main( int argc, char * argv[] ) {
		if ( argc != 6 ) {
				printf( "Usage: %s mapping-file query-fastIdx-file query-pepTree-file subject-fastIdx-file subject-pepTree-file\n"
				      , argv[0]
				      );
				return 1;
		}

		MMappedFastIdx queryFastIdx( argv[2] );
		MMappedPepTree queryPepTree( argv[3] );
		MMappedFastIdx subjectFastIdx( argv[4] );
		MMappedPepTree subjectPepTree( argv[5] );

		size_t szQuery, szSubject;
		queryPepTree  .ForLeaf( 0, [&]( char const * s, uint32_t ) {   szQuery   = strlen( s );   } );
		subjectPepTree.ForLeaf( 0, [&]( char const * s, uint32_t ) {   szSubject = strlen( s );   } );
		if ( szQuery != szSubject ) {
				printf( "Invalid pepTree files, not same words' size (query: %zu and subject: %zu)\n", szQuery, szSubject );
				return 1;
		}
		printf( "Words' size: %zu\n", szQuery );

		FILE * mappingFile = fopen( argv[1], "rb" );
		while( !feof( mappingFile ) ) {
				size_t queryIndex, subjectIndex;
				double score;
				fscanf( mappingFile, "%zu %zu %lg\n", &queryIndex, &subjectIndex, &score );

				subjectPepTree.ForLeaf( subjectIndex, [&]( char const *, uint32_t offset ) {
						subjectPepTree.ForLeafPos( offset, ProtFunctor( subjectFastIdx, szQuery ) );
				});
		}
		fclose( mappingFile );

		ostringstream ostr;
		ostr << argv[1] << ".profiles";
		FILE * outputFile = fopen( ostr.str().c_str(), "wb" );
		for_each( protProfiles, [=]( protProfile_map::value_type const & p ) {
				fprintf( outputFile, "%s\t", subjectFastIdx.GetName( p.first ) );
				for_each( p.second, [=]( unsigned int v ) {
						fprintf( outputFile, "%u ", v );
				});
				fprintf( outputFile, "\n" );
		});
		fclose( outputFile );

		return 0;
}
