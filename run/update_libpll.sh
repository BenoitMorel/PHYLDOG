
phyldogpath=$('pwd')/..
cd /home/morelbt/github/raxml-ng/libs/pll-modules/libs/libpll
make
cp src/.libs/*.so* ${phyldogpath}/dependencies/lib/
cp src/pll.h $(phyldogpath)/dependencies/include/pllmodules


