
if [ "$#" -ne 4 ]; then
  echo "Illegal number of parameters"
  exit
fi

prefix=$1
species=$2
genes=$3
hackmode=$4

outputdir=${prefix}_${species}_${genes}_${hackmode}
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


sed -i "s#RESULT=/home/morelbt/github/PHYLDOG/ExampleData/ResultFiles/#RESULT=${fullpath}/#g" ${fullpathdata}/GeneralOptions.txt
sed -i "s#OPT=/home/morelbt/github/PHYLDOG/ExampleData/OptionFiles/#OPT=${fullpathdata}/#g" ${fullpathdata}/GeneralOptions.txt

echo "rearrangement.gene.tree=nni" >> ${fullpathdata}/GeneralOptions.txt
echo "reset.gene.trees=no" >> ${fullpathdata}/GeneralOptions.txt

mpirun -np 4 ../../build/bin/phyldog param=${fullpathdata}/GeneralOptions.txt &> ${fullpath}/logs.txt




