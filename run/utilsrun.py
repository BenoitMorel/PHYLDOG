import os
import sys
import shutil
import shlex



def resizeFile(inputFile, newNumberOfLines):
  with  open(inputFile) as readFile:
    lines = readFile.readlines()
  with open(inputFile,'w') as w:
    w.writelines([item for item in lines[:newNumberOfLines]])

def replaceLinesStartingWith(inputFile, start, newLine):
  with open(inputFile) as readFile:
    lines = readFile.readlines()
  with open(inputFile,'w') as w:
    for line in lines:
      if line.startswith(start):
        w.write(newLine + "\n")
      else:
        w.write(line)
  
def fixListGenes(inputFile, newOptionsDir):
  with  open(inputFile) as readFile:
    lines = readFile.readlines()
  with open(inputFile,'w') as w:
    for line in lines:
      f = line.split('/')[-1]
      w.write(os.path.join(newOptionsDir,f)) 


def preparePhydlogFiles(prefix, suffix, dataset, seed, speciesNumber, genesNumber, method, useBestTrees, geneTreeSuffix):
  """ 
  Create a new directory, extract and modify an original dataset to this directory, and return the directory name
  prefix and suffix are appended at the begining and the end of the new directory name
  """
  # get paths
  path = os.path.dirname(os.path.realpath(__file__))
  phyldogDir = os.path.dirname(path)
  outputDir = prefix + "_" + dataset + "_" + seed + "_" + speciesNumber + "_" + genesNumber + "_" + method + suffix
  outputDir = os.path.join(path, outputDir)
  originDataDir = os.path.join(phyldogDir, "benoitdata", dataset)
  originOptionsDir = os.path.join(originDataDir, "OptionFiles")
  newDataDir = os.path.join(outputDir, dataset)
  newOptionsDir = os.path.join(newDataDir, "OptionFiles")
  newGeneralOptionsFile = os.path.join(newOptionsDir, "GeneralOptions.txt")
  geneTreesDir = os.path.join(originDataDir, "RaxmlTrees")
  resultsDir = os.path.join(outputDir, "results")
  # create directories
  shutil.rmtree(outputDir, True)
  os.mkdir(outputDir)
  os.mkdir(newDataDir)
  os.mkdir(resultsDir)
  # build the options files
  shutil.copytree(originOptionsDir, newOptionsDir)
  resizeFile(os.path.join(newOptionsDir, "listSpecies.txt"), int(speciesNumber))
  resizeFile(os.path.join(newOptionsDir, "listGenes.txt"), int(genesNumber))
  fixListGenes(os.path.join(newOptionsDir, "listGenes.txt"), newOptionsDir)
  replaceLinesStartingWith(newGeneralOptionsFile, "RESULT=", "RESULT=" + outputDir + "/") 
  replaceLinesStartingWith(newGeneralOptionsFile, "OPT=", "OPT=" + newOptionsDir + "/") 
  with open(newGeneralOptionsFile, "a") as w:
    w.write("rearrangement.gene.tree=spr\n")
    w.write("reset.gene.trees=no\n")
    w.write("likelihood.evaluator=" + method + "\n")
    w.write("seed=" + seed + "\n")
  # build per gene files
  if (useBestTrees):
    optionsFilenames = (opt for opt in os.listdir(newOptionsDir) if opt.endswith(".opt"))
    for opt in optionsFilenames:
      geneName = os.path.splitext(opt)[0]
      print(geneName)
      optFile = os.path.join(newOptionsDir, opt)
      geneTreeFile = os.path.join(geneTreesDir, geneName + geneTreeSuffix)
      replaceLinesStartingWith(optFile, "init.gene.tree=bionj", "init.gene.tree=user")
      replaceLinesStartingWith(optFile, "RESULT=", "RESULT=" + resultsDir + "/")
      replaceLinesStartingWith(optFile, "gene.tree.file=", "gene.tree.file=" + geneTreeFile)
  return outputDir
