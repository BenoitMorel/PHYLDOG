

phyldogpath=$('pwd')/..
cd /home/morelbt/github/raxml-ng/libs/pll-modules/
make
cp src/*/.libs/*.so* ${phyldogpath}/dependencies/lib/
cp src/*/*.h $(phyldogpath)/dependencies/include/pllmodules/

