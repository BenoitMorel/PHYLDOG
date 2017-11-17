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
  nodes=$(echo "($threads - 1)/16 + 1" | bc)
  submit_file="physub_${prefix}_${species}_${genes}_${hackmode}_${threads}"
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


  
  cat submitprefix >> $submit_file

  sbatch $submit_file
  rm $submit_file
}



generate_small()
{
  prefix=smallCluster
  generate_submit 42 $prefix 10 3 PLL 4
  generate_submit 42 $prefix 10 3 LIBPLL2 4
  generate_submit 20 $prefix 10 3 PLL 4
  generate_submit 20 $prefix 10 3 LIBPLL2 4
  generate_submit 30 $prefix 10 3 PLL 4
  generate_submit 30 $prefix 10 3 LIBPLL2 4
}

generate_medium()
{
  prefix=mediumCluster
  generate_submit 42 $prefix 20 10 PLL 4
  generate_submit 42 $prefix 20 10 LIBPLL2 4
  generate_submit 20 $prefix 20 10 PLL 4
  generate_submit 20 $prefix 20 10 LIBPLL2 4
  generate_submit 30 $prefix 20 10 PLL 4
  generate_submit 30 $prefix 20 10 LIBPLL2 4
}

generate_big()
{
  prefix=bigCluster
  generate_submit 42 $prefix 30 20 PLL 10
  generate_submit 42 $prefix 30 20 LIBPLL2 10
  generate_submit 20 $prefix 30 20 PLL 10
  generate_submit 20 $prefix 30 20 LIBPLL2 10
  generate_submit 30 $prefix 30 20 PLL 10
  generate_submit 30 $prefix 30 20 LIBPLL2 10
}

generate_verybig()
{
  prefix=veryBigCluster
  generate_submit 42 $prefix 50 40 PLL 16
  generate_submit 42 $prefix 50 40 LIBPLL2 16
  generate_submit 20 $prefix 50 40 PLL 16
  generate_submit 20 $prefix 50 40 LIBPLL2 16
  generate_submit 30 $prefix 50 40 PLL 16
  generate_submit 30 $prefix 50 40 LIBPLL2 16
}

generate_full()
{
  prefix=fullCluster
  #generate_submit 42 $prefix 55 40 PLL 16
  generate_submit 42 $prefix 55 40 LIBPLL2 16
  #generate_submit 20 $prefix 55 40 PLL 16
  generate_submit 20 $prefix 55 40 LIBPLL2 16
  #generate_submit 30 $prefix 55 40 PLL 16
  generate_submit 30 $prefix 55 40 LIBPLL2 16
}

#generate_small
#generate_medium
#generate_big
#generate_verybig
generate_full

#
