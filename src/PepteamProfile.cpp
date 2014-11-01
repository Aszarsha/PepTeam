#include <cstdint>
#include <cstring>
#include <iostream>
#include <boost/range/algorithm/for_each.hpp>

#include "ProgramOptions.hpp"
#include "File.hpp"
#include "FastIdx.hpp"
#include "PepTree.hpp"

using namespace std;
using boost::range::for_each;

#define UsageFunction( function )                                            \
	do {                                                                      \
			function( "Usage: ", argv[0]                                        \
			        , " mapping-file query-fastIdx-file query-pepTree-file"     \
			          " subject-fastIdx-file subject-pepTree-file\n"            \
			        , baseOption                                                \
			        );                                                          \
	} while ( false )

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

int main( int argc, char * argv[] ) try {
		po::options_description baseOption( "Options", 999, 999 );
		baseOption.add_options()
			( "help,h"  , "Print this help message and exit" )
			( "output,o", po::value< string >(), "Output file (default: mapping-file+\".profiles\")" )
			;

		po::options_description hiddenOptions( "Required options", 999, 999 );
		hiddenOptions.add_options()
			( "mapping-file"   , po::value< string >()->required(), "Input mapping file" )
			( "query-fastIdx"  , po::value< string >()->required(), "Query repertoire .fastIdx input file" )
			( "query-pepTree"  , po::value< string >()->required(), "Query repertoire .pepTree input file" )
			( "subject-fastIdx", po::value< string >()->required(), "Subject proteome .fastIdx input file" )
			( "subject-pepTree", po::value< string >()->required(), "Subject proteome .pepTree input file" )
			;
		po::positional_options_description posOptions;
		posOptions.add( "mapping-file"   , 1 );
		posOptions.add( "query-fastIdx"  , 1 );
		posOptions.add( "query-pepTree"  , 1 );
		posOptions.add( "subject-fastIdx", 1 );
		posOptions.add( "subject-pepTree", 1 );

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

		MMappedFastIdx queryFastIdx( vm["query-fastIdx"].as< string >().c_str() );
		MMappedPepTree queryPepTree( vm["query-pepTree"].as< string >().c_str() );
		MMappedFastIdx subjectFastIdx( vm["subject-fastIdx"].as< string >().c_str() );
		MMappedPepTree subjectPepTree( vm["subject-pepTree"].as< string >().c_str() );

		size_t szQuery, szSubject;
		queryPepTree  .ForLeaf( 0, [&]( char const * s, uint32_t ) {   szQuery   = strlen( s );   } );
		subjectPepTree.ForLeaf( 0, [&]( char const * s, uint32_t ) {   szSubject = strlen( s );   } );
		if ( szQuery != szSubject ) {
				printf( "Invalid pepTree files, not same words' size (query: %zu and subject: %zu)\n", szQuery, szSubject );
				return 1;
		}
		printf( "Words' size: %zu\n", szQuery );

		auto mappingFile = OpenFile( vm["mapping-file"].as< string >().c_str(), "rb" );
		while( !feof( mappingFile ) ) {
				size_t queryIndex, subjectIndex;
				double score;
				fscanf( mappingFile, "%zu %zu %lg\n", &queryIndex, &subjectIndex, &score );

				subjectPepTree.ForLeaf( subjectIndex, [&]( char const *, uint32_t offset ) {
						subjectPepTree.ForLeafPos( offset, ProtFunctor( subjectFastIdx, szQuery ) );
				});
		}
		mappingFile.close();

		string outFileName = vm.count( "output" ) ? vm["output"].as< string >()
		                                          : vm["input-file"].as< string >() + ".profiles";
		auto outFile = OpenFile( outFileName.c_str(), "wb" );
		for_each( protProfiles, [=,&outFile]( protProfile_map::value_type const & p ) {
				fprintf( outFile, "%s\t", subjectFastIdx.GetName( p.first ) );
				for_each( p.second, [=,&outFile]( unsigned int v ) {
						fprintf( outFile, "%u ", v );
				});
				fprintf( outFile, "\n" );
		});
		outFile.close();

		return 0;
} catch ( std::exception & e ) {
		fprintf( stderr, "Unhandled exception reached main: %s, application will now exit\n", e.what() );
		return 271828;
} catch ( ... ) {
		fprintf( stderr, "Unhandled object reached main, application will now exit\n" );
		return 314159;
}
