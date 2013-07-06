#include <cstdio>
#include <string>
#include <sstream>
#include <vector>
#include <cstring>
#include <cstdint>
#include <limits>
#include <stdexcept>

#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>

#include "FastIdx.hpp"

using namespace std;

MemFastIdx::MemFastIdx( vector< string > const & proteinsName
                      , vector< string > const & proteinsSeq
                      ) {
		auto & indices   = GetIndices();
		auto & names     = GetNames();
		auto & sequences = GetSequences();

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
}

static auto const indicesOffset = static_cast< uint32_t >( 2 * sizeof( uint32_t ) );

void MemFastIdx::Write( FILE * file ) const {
		auto const & indices   = GetIndices();
		auto const & names     = GetNames();
		auto const & sequences = GetSequences();

		auto snOffset  = static_cast< uint32_t >( indicesOffset + indices.size() * sizeof( indices.front() ) );
		auto seqOffset = static_cast< uint32_t >( snOffset      + names.size()   * sizeof( names.front() ) );

		uint32_t arr[] = { snOffset, seqOffset };
		fwrite( arr, sizeof( arr[0] ), sizeof( arr ) / sizeof( arr[0] ), file );
		fwrite( indices.data(), sizeof( indices.front() ), indices.size(), file );
		fwrite( names.data(), sizeof( names.front() ), names.size(), file );
		fwrite( sequences.data(), sizeof( sequences.front() ), sequences.size(), file );
}

MMappedFastIdx::MMappedFastIdx( const char * filename )
	: fd( open( filename, O_RDONLY ) ) {
		if ( !fd ) {
				throw std::runtime_error{ "Unable to open input FastIdx file\n" };
		}
		struct stat fStat;
		fstat( fd, &fStat );
		fileSize = static_cast< uint32_t >( fStat.st_size );
		ptr = static_cast< char const * >( mmap( nullptr, fileSize, PROT_READ, MAP_SHARED, fd, 0 ) );
}

MMappedFastIdx::~MMappedFastIdx() {
		close( fd );
}


size_t MMappedFastIdx::Size() const {
		return GetIndicesSize()/2;
}

size_t MMappedFastIdx::GetIndicesSize() const {
		return (reinterpret_cast< uint32_t const * >( ptr )[0] - indicesOffset)/sizeof( uint32_t );
}

size_t MMappedFastIdx::GetNamesSize() const {
		uint32_t const * tmp = reinterpret_cast< uint32_t const * >( ptr );
		return tmp[1] - tmp[0];
}

size_t MMappedFastIdx::GetSequencesSize() const {
		return fileSize - reinterpret_cast< uint32_t const * >( ptr )[1];
}

uint32_t const * MMappedFastIdx::GetIndices() const {
		return reinterpret_cast< uint32_t const * >( ptr + indicesOffset );
}

char const * MMappedFastIdx::GetNames() const {
		return ptr + reinterpret_cast< uint32_t const * >( ptr )[0];
}

char const * MMappedFastIdx::GetSequences() const {
		return ptr + reinterpret_cast< uint32_t const * >( ptr )[1];
}
