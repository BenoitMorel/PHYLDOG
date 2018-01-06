import reccommon as rc


iterations = "1"
for i in range(0, 10):
  rc.bench_gene("HBG01100" + str(i), iterations)

for i in range(10, 40):
  rc.bench_gene("HBG0110" + str(i), iterations)


