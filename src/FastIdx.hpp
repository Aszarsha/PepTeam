#ifndef FASTIDX_HPP
#define FASTIDX_HPP

#include <cstdint>
#include <string>
#include <vector>
#include <tuple>

class MemFastIdx {
	public:
		MemFastIdx( std::vector< std::string > const & proteinsName
                , std::vector< std::string > const & proteinsSeq
                );

	public:
		void Write( FILE * file ) const;

	private:
		typedef std::tuple< std::vector< uint32_t >
		                  , std::vector< char >
		                  , std::vector< char >
		                  > LinearizedProtIdx;

	private:
		LinearizedProtIdx idx;

	private:
		std::vector< uint32_t >       & GetIndices()       {   return std::get< 0 >( idx );   }
		std::vector< uint32_t > const & GetIndices() const {   return std::get< 0 >( idx );   }

		std::vector< char >       & GetNames()       {   return std::get< 1 >( idx );   }
		std::vector< char > const & GetNames() const {   return std::get< 1 >( idx );   }

		std::vector< char >       & GetSequences()       {   return std::get< 2 >( idx );   }
		std::vector< char > const & GetSequences() const {   return std::get< 2 >( idx );   }
};

class MMappedFastIdx {
	public:
		explicit MMappedFastIdx( const char * filename );

		~MMappedFastIdx();

	public:
		size_t           Size() const;

		uint32_t const * GetIndicesData() const;
		size_t           GetIndicesSize() const;

		char const     * GetNamesData() const;
		size_t           GetNamesSize() const;

		char const     * GetSequencesData() const;
		size_t           GetSequencesSize() const;

		char const * GetName    ( size_t index ) const {   return GetNamesData    () + GetIndicesData()[ 2*index   ];   }
		char const * GetSequence( size_t index ) const {   return GetSequencesData() + GetIndicesData()[ 2*index+1 ];   }

	private:
		int fd;
		uint32_t fileSize;
		char const * ptr;
};

#endif
