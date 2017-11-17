
repo=libpll2_repo
dep=libpll2_dependencies
libs=${dep}/libs
include=${dep}/include
pllinclude=$include/pll
pllmodulesinclude=$include/pllmodules

#### CLONE #####
rm -rf $repo
mkdir $repo
cd $repo
git clone --recursive git@github.com:ddarriba/pll-modules.git
cd ..
# small hack before we fix this in pllmodules
# otherwise NNI moves rollback crash in a wrong assert..
# todobenoit.
sed -i 's/assert(UNIMPLEMENTED);//g' libpll2_repo/pll-modules/src/tree/pll_tree.c



#### BUILD ####
cd $repo/pll-modules/libs/libpll/
./autogen.sh
./configure
make

cd ../../
./install-with-libpll.sh libs/libpll/

cd ../../
rm -rf $dep
mkdir $dep
mkdir $libs
mkdir $include
mkdir $pllinclude

cp $repo/pll-modules/libs/libpll/src/.libs/*.so*  $libs
cp $repo/pll-modules/src/*/.libs/*so* $libs

cp $repo/pll-modules/libs/libpll/src/*.h $pllinclude
cp $repo/pll-modules/src/*/*.h $pllinclude

echo "Libraries were copied to $libs"
echo "Headers were copied to $pllinclude"
