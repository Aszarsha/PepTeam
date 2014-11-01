#ifndef FILE_HPP
#define FILE_HPP

#include <cstdio>
#include <cstdlib>
#include <memory>

struct AutoFILE {
	public:
		AutoFILE( FILE * file ): ptr( file, fclose ) {   }
		AutoFILE( FILE * file, int (*deleter)( FILE * ) ): ptr( file, deleter ) {   }

		void close() {   ptr.reset();   }

		operator FILE *() const {   return ptr.get();   }

	private:
		using unique_FILE_ptr = std::unique_ptr< FILE, int (*)( FILE * ) >;

	private:
		unique_FILE_ptr ptr;
};

auto OpenFile( char const * fileName, char const * mode ) {
		FILE * inputFile = fopen( fileName, mode );
		if ( !inputFile ) {
				fprintf( stderr, "Unable to open input file \"%s\"\n", fileName );
				exit( 1 );
		}
		return AutoFILE( inputFile );
}

#endif
