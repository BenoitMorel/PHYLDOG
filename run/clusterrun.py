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
startingTrees = sys.argv[7]
cores=sys.argv[8]

suffix = "_" + cores

#prepare files
outputDir = utilsrun.preparePhydlogFiles(prefix, suffix, dataset, seed, speciesNumber, genesNumber, method, startingTrees)
newGeneralOptionsFile = os.path.join(outputDir, dataset, "OptionFiles", "GeneralOptions.txt")
logFile = os.path.join(outputDir, "logs.txt")
executable = "../../build/bin/phyldog"
submitName = "phyldog.sh"
submit = utilscluster.createPhyldogSubmit(outputDir, submitName, executable, newGeneralOptionsFile, cores)
utilscluster.launchSubmit(submit)





