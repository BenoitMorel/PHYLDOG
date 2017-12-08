import sys
import os
import subprocess
import shlex

def createSubmit(fullPath, cores):
  nodes = str((int(cores) - 1)//16 + 1)
  with open(fullPath, "w") as w:
    w.write("#!/bin/bash\n")
    w.write("#SBATCH -o ng_%j.out\n")
    w.write("#SBATCH -N " + nodes + "\n")
    w.write("#SBATCH -n " + cores + "\n")
    w.write("#SBATCH -B 2:8:1\n")
    w.write("#SBATCH --threads-per-core=1\n")
    w.write("#SBATCH --cpus-per-task=1\n")
    w.write("#SBATCH -t 24:00:00\n")
    w.write("\n")
    w.write("export SCOREP_PROFILING_MAX_CALLPATH_DEPTH=100\n")
    w.write("export SCOREP_PROFILING_ENABLE_CORE_FILES=1\n")
    w.write("\n")

def createPhyldogSubmit(directory, name, executable, generalOptionsFile, cores):
  fullPath = os.path.join(directory, name)
  logFile = os.path.join(directory, "logs.txt")
  createSubmit(fullPath, cores)
  with open(fullPath, "a") as w:
    w.write("mpirun -np " + cores + " " + executable)
    w.write(" param=" + generalOptionsFile)
    w.write(" &> " + logFile)
  with open(fullPath) as r:
    lines = r.readlines()
    for line in lines:
      print(line)
  return fullPath

def launchSubmit(fullPath):
 command = "sbatch " + os.path.basename(fullPath)
 os.chdir(os.path.dirname(fullPath))
 subprocess.call(shlex.split(command))

