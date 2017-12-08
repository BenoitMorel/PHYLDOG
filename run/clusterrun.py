import sys
import os
import shutil
import subprocess
import shlex
import utilsrun
import utilscluster

# parse command line
if (len(sys.argv) < 8):
  print("Invalid syntax")
prefix=sys.argv[1]
dataset=sys.argv[2]
seed=sys.argv[3]
speciesNumber=sys.argv[4]
genesNumber=sys.argv[5]
method=sys.argv[6]
cores=sys.argv[7]
useBestTrees = True
geneTreeSuffix = ".raxml.startTree"

suffix = "_" + cores

#prepare files
outputDir = utilsrun.preparePhydlogFiles(prefix, suffix, dataset, seed, speciesNumber, genesNumber, method, useBestTrees, geneTreeSuffix)
newGeneralOptionsFile = os.path.join(outputDir, dataset, "OptionFiles", "GeneralOptions.txt")
logFile = os.path.join(outputDir, "logs.txt")
executable = "../../build/bin/phyldog"
submitName = "phyldog.sh"
submit = utilscluster.createPhyldogSubmit(outputDir, submitName, executable, newGeneralOptionsFile, cores)
utilscluster.launchSubmit(submit)




# gogogo
#os.chdir(outputDir)
#runCommand = "mpirun -np 4 " + executable + " param=" + newGeneralOptionsFile   
#printCommand = './show_summary.sh ' + outputDir
#print("\n" + runCommand + "\n")
#with open(logFile,"wb") as f:
#  subprocess.call(shlex.split(runCommand), stdout=f, stderr=f)   
#os.chdir(path)
#os.system(printCommand)


