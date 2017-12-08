
if [ "$#" -ne 5 ]; then
  echo "Illegal number of parameters"
  exit
fi

prefix=$1
seed=$2
species=$3
genes=$4
method=$5

outputdir=${prefix}_${seed}_${species}_${genes}_${method}
rm -rf $outputdir
mkdir $outputdir
cd $outputdir
fullpath=$('pwd')
optionspath="../../ExampleData/OptionFiles/"
fullpathdata=${fullpath}/ExampleData
mkdir $fullpathdata
cp ${optionspath}/GeneralOptions.txt $fullpathdata
head -n $genes  ${optionspath}/listGenes.txt > $fullpathdata/listGenes.txt
head -n $species ${optionspath}/listSpecies.txt > $fullpathdata/listSpecies.txt
sed -i "s#RESULT=.*#RESULT=${fullpath}/#g" ${fullpathdata}/GeneralOptions.txt
sed -i "s#OPT=.*/#OPT=${fullpathdata}/#g" ${fullpathdata}/GeneralOptions.txt
echo "rearrangement.gene.tree=spr" >> ${fullpathdata}/GeneralOptions.txt
echo "reset.gene.trees=no" >> ${fullpathdata}/GeneralOptions.txt
echo "likelihood.evaluator=${method}" >> ${fullpathdata}/GeneralOptions.txt
echo "seed=${seed}" >> ${fullpathdata}/GeneralOptions.txt
export SCOREP_PROFILING_MAX_CALLPATH_DEPTH=50



mpirun -np 4 ../../build/bin/phyldog param=${fullpathdata}/GeneralOptions.txt &> ${fullpath}/logs.txt
cd ..
./show_summary.sh $outputdir



