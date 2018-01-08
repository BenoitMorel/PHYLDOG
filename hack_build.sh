mkdir -p build
cd build
cmake .. -DCMAKE_LIBRARY_PATH="/home/morelbt/github/PHYLDOG/dependencies/lib/" -DCMAKE_INCLUDE_PATH="/home/morelbt/github/PHYLDOG/dependencies/include" -DBOOST_LIBRARYDIR="/hits/sw/shared/apps/Boost/1.58.0-goolf-1.7.20/lib" -DBOOST_ROOT=/hits/sw/shared/apps/Boost/1.58.0-goolf-1.7.20/ -DCMAKE_CXX_STANDARD=11

filesToReplace=( "src/CMakeFiles/benchreconciliation.dir/build.make" "src/CMakeFiles/benchreconciliation.dir/link.txt" "src/CMakeFiles/phyldog.dir/build.make" "src/CMakeFiles/phyldog.dir/link.txt")

#gccPrefix="skin "
#gccSuffix=" -g -DSCOREP"
#gccPrefix=""
gccSuffix=" -g "

for fileToReplace in "${filesToReplace[@]}"
do
  sed -i "s#/hits/sw/shared/apps/GCCcore/5.4.0/bin/c++#$gccPrefix /hits/sw/shared/apps/GCCcore/5.4.0/bin/c++ $gccSuffix #g" $fileToReplace
done

sed -i '/depends/c\ echo "benoit hack: skipping dependencies update"' src/CMakeFiles/phyldog.dir/build.make
#sed -i '/depends/c\ echo "benoit hack: skipping dependencies update"' src/CMakeFiles/phyldog_light.dir/build.make
sed -i '/depends/c\ echo "benoit hack: skipping dependencies update"' src/CMakeFiles/benchreconciliation.dir/build.make


make
