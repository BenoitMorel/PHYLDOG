#!/bin/bash


submit_file=
output_dir=test_launcher


write()
{
  echo $1 >> $submit_file
}

generate_submit()
{
  prefix=$1
  species=$2
  genes=$3
  method=$4
  threads=$5
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

  write "prefix=$prefix"
  write "species=$species"
  write "genes=$genes"
  write "method=$method"
  write "threads=$threads"


  
  cat submitprefix >> $submit_file

  sbatch $submit_file
  rm $submit_file
}




generate_submit small 10 3 BPP 4
generate_submit small 10 3 PLL 4
generate_submit small 10 3 LIBPLL2 4

#generate_submit medium 20 10 PLL 4
#generate_submit medium 20 10 LIBPLL2 4

#generate_submit all 55 40 PLL 16
#generate_submit all 55 40 LIBPLL2 16

#generate_submit all 55 40 PLL 40
#generate_submit all 55 40 LIBPLL2 40



