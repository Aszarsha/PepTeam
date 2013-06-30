#ifndef FASTIDX_HPP
#define FASTIDX_HPP

#include <cstdint>
#include <string>
#include <vector>
#include <tuple>

typedef std::tuple< std::vector< uint32_t >
                  , std::vector< char >
                  , std::vector< char >
                  > LinearizedProtIdx;

inline std::vector< uint32_t > & GetProtIdxIndices( LinearizedProtIdx & idx ) {
		return std::get< 0 >( idx );
}
inline std::vector< uint32_t > const & GetProtIdxIndices( LinearizedProtIdx const & idx ) {
		return std::get< 0 >( idx );
}

inline std::vector< char > & GetProtIdxNames( LinearizedProtIdx & idx ) {
		return std::get< 1 >( idx );
}
inline std::vector< char > const & GetProtIdxNames( LinearizedProtIdx const & idx ) {
		return std::get< 1 >( idx );
}

inline std::vector< char > & GetProtIdxSequences( LinearizedProtIdx & idx ) {
		return std::get< 2 >( idx );
}
inline std::vector< char > const & GetProtIdxSequences( LinearizedProtIdx const & idx ) {
		return std::get< 2 >( idx );
}

LinearizedProtIdx LinearizeProtIdx( std::vector< std::string > const & proteinsName
                                  , std::vector< std::string > const & proteinsSeq
                                  );

void WriteLinearizedProtIdx( FILE * file, LinearizedProtIdx const & idx );

size_t ReadProteinIndexSize( FILE * file );
LinearizedProtIdx ReadProteinIndex( FILE * file );

#endif
