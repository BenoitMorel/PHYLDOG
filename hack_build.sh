mkdir -p build
cd build
cmake .. -DCMAKE_LIBRARY_PATH="/home/morelbt/github/PHYLDOG/dependencies/lib/" -DCMAKE_INCLUDE_PATH="/home/morelbt/github/PHYLDOG/dependencies/include" -DBOOST_LIBRARYDIR="/hits/sw/shared/apps/Boost/1.58.0-goolf-1.7.20/lib" -DBOOST_ROOT=/hits/sw/shared/apps/Boost/1.58.0-goolf-1.7.20/
sed -i '/depends/c\ echo "benoit hack: skipping dependencies update"' src/CMakeFiles/phyldog.dir/build.make
sed -i 's#/hits/sw/shared/apps/GCCcore/4.9.3/bin/c++#skin /hits/sw/shared/apps/GCCcore/4.9.3/bin/c++#g' src/CMakeFiles/phyldog.dir/build.make
sed -i 's#/hits/sw/shared/apps/GCCcore/4.9.3/bin/c++#skin /hits/sw/shared/apps/GCCcore/4.9.3/bin/c++#g' src/CMakeFiles/phyldog.dir/link.txt

make
