#ifndef PROGRAMOPTIONS_HPP
#define PROGRAMOPTIONS_HPP

#include <cstdio>
#include <ostream>
#include <iostream>
#include <boost/program_options.hpp>

namespace po = boost::program_options;

namespace impl {
	[[ noreturn ]] void PrintUsageAndExit( std::ostream & ostr, int exitCode ) {
			exit( exitCode );
	}

	template< typename PrintableType, typename... PrintableTypes >
	[[ noreturn ]] void PrintUsageAndExit( std::ostream & ostr, int exitCode
	                                     , PrintableType && printable, PrintableTypes &&... printables
	                      ) {
			ostr << printable;
			impl::PrintUsageAndExit( ostr, exitCode, std::forward< PrintableTypes >( printables )... );
	}

	bool DetectExclusiveOptions( po::variables_map const & vm, size_t optionsCount ) {
			return optionsCount == 1;
	}

	template< typename... ConstChar >
	bool DetectExclusiveOptions( po::variables_map const & vm, size_t optionsCount
	                           , char const * opt, ConstChar &&... options
	                           ) {
			return impl::DetectExclusiveOptions( vm, optionsCount + (vm.count( opt ) && !vm[opt].defaulted())
			                                   , std::forward< ConstChar >( options )...
			                                   );
	}
}

template< typename... Options >
bool DetectExclusiveOptions( po::variables_map const & vm, Options &&... options ) {
		return impl::DetectExclusiveOptions( vm, 0, std::forward< Options >( options )... );
}

template< typename... Types >
[[ noreturn ]] void Usage( Types &&... printables ) {
		impl::PrintUsageAndExit( std::cout, 0, std::forward< Types >( printables )... );
}

template< typename... Types >
[[ noreturn ]] void UsageError( Types &&... printables ) {
		impl::PrintUsageAndExit( std::cerr, 1, std::forward< Types >( printables )... );
}
#endif
