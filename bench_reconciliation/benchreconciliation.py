import subprocess






def bench_gene(gene):
  executable = "/home/morelbt/github/PHYLDOG/build/bin/benchreconciliation"
  speciesTree = "/home/morelbt/github/PHYLDOG/benoitdata/DataExample/speciesTree.newick"
  genesTree = "/home/morelbt/github/PHYLDOG/benoitdata/DataExample/RaxmlTrees/" + gene + ".raxml.bestTree"
  link = "/home/morelbt/github/PHYLDOG/benoitdata/DataExample/LinkFiles/" + gene + ".link"
  print(executable + " " + speciesTree + " " + genesTree + " " + link)
  try:
    subprocess.check_call([executable, speciesTree, genesTree, link])
  except:
    pass

bench_gene("HBG011000")

