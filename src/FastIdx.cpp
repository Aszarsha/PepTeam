#include <cstdio>
#include <string>
#include <sstream>
#include <vector>
#include <cstring>
#include <cstdint>
#include <limits>
#include <stdexcept>

#include "FastIdx.hpp"

using namespace std;

LinearizedProtIdx LinearizeProtIdx( vector< string > const & proteinsName
                                  , vector< string > const & proteinsSeq
                                  ) {
		LinearizedProtIdx idx;
		auto & indices   = GetProtIdxIndices( idx );
		auto & names     = GetProtIdxNames( idx );
		auto & sequences = GetProtIdxSequences( idx );

		auto ConcatInVector = []( vector< char > & vec, string const & seq ) -> size_t {
				auto curIndex = vec.size();
				vec.insert( end( vec ), begin( seq ), end( seq ) );
				vec.push_back( '\0' );
				return curIndex;
		};

		for ( size_t i = 0, e = proteinsName.size(); i != e; ++i ) {
				auto snPos  = ConcatInVector( names, proteinsName[i] );
				if ( snPos > numeric_limits<uint32_t>::max() ) {
						throw std::runtime_error{ "Protein name index overflow, abording" };
				}
				indices.push_back( static_cast< uint32_t >( snPos ) );
				auto seqPos = ConcatInVector( sequences, proteinsSeq[i] );
				if ( seqPos > numeric_limits<uint32_t>::max() ) {
						throw std::runtime_error{ "Protein sequence index overflow, abording" };
				}
				indices.push_back( static_cast< uint32_t >(  seqPos ) );
		}

		return idx;
}

static auto const indicesOffset = static_cast< uint32_t >( 2 * sizeof( uint32_t ) );

void WriteLinearizedProtIdx( FILE * file, LinearizedProtIdx const & idx ) {
		auto const & indices   = GetProtIdxIndices( idx );
		auto const & names     = GetProtIdxNames( idx );
		auto const & sequences = GetProtIdxSequences( idx );

		auto snOffset  = static_cast< uint32_t >( indicesOffset + indices.size() * sizeof( indices.front() ) );
		auto seqOffset = static_cast< uint32_t >( snOffset      + names.size()   * sizeof( names.front() ) );

		uint32_t arr[] = { snOffset, seqOffset };
		fwrite( arr, sizeof( arr[0] ), sizeof( arr ) / sizeof( arr[0] ), file );
		fwrite( indices.data(), sizeof( indices.front() ), indices.size(), file );
		fwrite( names.data(), sizeof( names.front() ), names.size(), file );
		fwrite( sequences.data(), sizeof( sequences.front() ), sequences.size(), file );
}

size_t ReadProteinIndexSize( FILE * file ) {
		uint32_t arr[2];
		fread( arr, sizeof( arr[0] ), sizeof( arr ) / sizeof( arr[0] ), file );

		return (arr[0] - indicesOffset)/(2 * sizeof( uint32_t ));
}

LinearizedProtIdx ReadProteinIndex( FILE * file ) {
		fseek( file, 0, SEEK_END );
		auto fileSize = ftell( file );
		fseek( file, 0, SEEK_SET );

		uint32_t arr[2];
		fread( arr, sizeof( arr[0] ), sizeof( arr ) / sizeof( arr[0] ), file );

		auto indices = vector< uint32_t >( (arr[0] - indicesOffset)/sizeof( uint32_t ) );
		fread( indices.data(), sizeof( indices.front() ), indices.size(), file );

		auto names = vector< char >( arr[1] - arr[0] );
		fread( names.data(), sizeof( names.front() ), names.size(), file );

		auto sequences = vector< char >( fileSize - arr[1] );
		fread( sequences.data(), sizeof( sequences.front() ), sequences.size(), file );

		return LinearizedProtIdx{ indices, names, sequences };
}
