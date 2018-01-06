import subprocess






def bench_gene(gene, iterations):
  executable = "/home/morelbt/github/PHYLDOG/build/bin/benchreconciliation"
  speciesTree = "/home/morelbt/github/PHYLDOG/benoitdata/DataExample/speciesTree.newick"
  genesTree = "/home/morelbt/github/PHYLDOG/benoitdata/DataExample/RaxmlTrees/" + gene + ".raxml.bestTree"
  link = "/home/morelbt/github/PHYLDOG/benoitdata/DataExample/LinkFiles/" + gene + ".link"
  useFastRec = "1"


  #print(executable + " " + speciesTree + " " + genesTree + " " + link + " " + useFastRec + " " + iterations)
  try:
    subprocess.check_call([executable, speciesTree, genesTree, link, useFastRec, iterations])
  except:
    pass
  
  

