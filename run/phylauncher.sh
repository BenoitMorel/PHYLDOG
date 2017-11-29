#!/bin/bash


submit_file=
output_dir=test_launcher


write()
{
  echo $1 >> $submit_file
}

generate_submit()
{
  seed=$1
  prefix=$2
  species=$3
  genes=$4
  method=$5
  threads=$6
  dataset=$7
  nodes=$(echo "($threads - 1)/16 + 1" | bc)
  submit_file="s_${prefix}_${dataset}_${seed}_${species}_${genes}_${method}_${threads}"
  rm -f  $submit_file
  write  "#!/bin/bash"
  write "#SBATCH -o ng_$submit_file_%j.out"
  write "#SBATCH -N $nodes"
  write "#SBATCH -n $threads"
  write "#SBATCH -B 2:8:1"
  write "#SBATCH --threads-per-core=1"
  write "#SBATCH --cpus-per-task=1"
  write "#SBATCH -t 24:00:00"

  write "seed=$seed"
  write "prefix=$prefix"
  write "species=$species"
  write "genes=$genes"
  write "method=$method"
  write "threads=$threads"
  write "dataset=$dataset"

  echo $threads
  
  cat submitprefix >> $submit_file

  sbatch $submit_file
  rm $submit_file
}


libpll_only()
{

  prefix=$1
  species=$2
  genes=$3
  threads=$4
  dataset=$5
  generate_submit 10 $prefix $species $genes LIBPLL2 $threads $dataset
  #generate_submit 20 $prefix $species $genes LIBPLL2 $threads $dataset
  #generate_submit 30 $prefix $species $genes LIBPLL2 $threads $dataset
}

pll_only()
{

  prefix=$1
  species=$2
  genes=$3
  threads=$4
  dataset=$5
  generate_submit 10 $prefix $species $genes PLL $threads $dataset
  #generate_submit 20 $prefix $species $genes PLL $threads $dataset
  #generate_submit 30 $prefix $species $genes PLL $threads $dataset
}

compare_methods()
{

  prefix=$1
  species=$2
  genes=$3
  threads=$4
  dataset=$5
  generate_submit 10 $prefix $species $genes LIBPLL2 $threads $dataset
  #generate_submit 20 $prefix $species $genes LIBPLL2 $threads $dataset
  #generate_submit 30 $prefix $species $genes LIBPLL2 $threads $dataset
  generate_submit 10 $prefix $species $genes PLL $threads $dataset
  #generate_submit 20 $prefix $species $genes PLL $threads $dataset
  #generate_submit 30 $prefix $species $genes PLL $threads $dataset
}

#templates
#compare_methods smallClusterFullOptNoReset 15 50 4 DataCarine
#compare_methods mediumClusterFullOptNoReset 15 500 128 DataCarine
#compare_methods fullClusterFullOptNoReset 15 8880 512 DataCarine
#pll_only smallClusterNoResetBastien 15 50 4 DataCarine
#pll_only mediumClusterNoResetBastien 15 500 128 DataCarine
#pll_only fullClusterNoResetBastien 15 8880 512 DataCarine

# for tests
#libpll_only smallClusterFullOptNoResetScalasca 15 150 16 DataCarine

#to analyse
#libpll_only fullClusterNoResetScalasca 15 8880 128 DataCarine
#libpll_only fullClusterNoResetScalasca 15 8880 256 DataCarine
#libpll_only fullClusterNoResetScalasca 15 8880 512 DataCarine
#libpll_only fullClusterNoResetScalasca 15 8880 1024 DataCarine
#running
#libpll_only fullClusterNoResetScalasca 15 8880 64 DataCarine
#libpll_only fullClusterNoResetScalasca 15 8880 32 DataCarine
#libpll_only fullClusterNoResetScalasca 15 8880 16 DataCarine
libpll_only fullClusterNoResetScalascaLikelihoodIt 15 8880 256 DataCarine
libpll_only fullClusterNoResetScalascaLikelihoodIt 15 8880 1024 DataCarine1024

