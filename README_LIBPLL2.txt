

Before installing libpll and pll-modules, make sure you have the following dependencies/tools installed:
apt-get install flex bison autotools-dev autoconf libtool

Clone and build pll-modules and libpll-2:
./build_libpll2.sh

Then add the libraries and includes to -DCMAKE_LIBRARY_PATH and -DCMAKE_INCLUDE_PATH when calling cmake
