#!/bin/bash

default_dir="../CAMI_experiments"

usage="
USAGE:\n
                $0 [path/to/CAMI_experiments/=$default_dir] [-t num_threads=16] [-o overwrite='n']
\n"

if [ $1 == "-h" ]; then
	echo -e $usage
	exit
fi

# Set arguments
if [ ! -z $1 ]; then
  analysis_dir=$1
else
  analysis_dir=$default_dir
fi

while getopts ":t:o:" opt; do
  case $opt in
  t) n_cpu="$OPTARG"
    ;;
  o) overwrite="$OPTARG"
    ;;
  \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done
if [ -z $n_cpu ]; then
  n_cpu=16
fi
if [ -z $overwrite ]; then
  overwrite="n"
fi

analysis_refpkgs=${analysis_dir}/refpkgs/
file_array=("/mnt/sdb/Hallam_projects/Hallam_Databases/formatted/CAMI_I_HIGH/gold_standard_high_single.fasta.gz"
"${analysis_dir}/RH_S001_merged.fq"
"${analysis_dir}/RH_S001_forward.fq")
outputs_array=()

if [ ! -d $analysis_refpkgs ]; then
  echo "ERROR: Unable to find directory with reference packages '${analysis_refpkgs}'"
  exit
fi

for f in "${file_array[@]}"
do
  if [ ! -f $f ]; then
    echo "ERROR: Unable to find input file '${f}'"
    exit
  fi

  assign_output=$analysis_dir/$( basename $f | sed 's/.fq.*//g' | sed 's/.fasta.*//g' )
  if [ -d ${assign_output} ] && [ $overwrite == 'n' ]; then
    echo "Skipping classification of '${f}'"
  else
    treesapp assign \
    -i $f \
    -o ${assign_output} \
    --refpkg_dir ${analysis_refpkgs} \
    --num_procs $n_cpu \
    --trim_align \
    --overwrite --delete
  fi
  outputs_array+=(${assign_output})
done

for f in ${analysis_refpkgs}/*pkl
do
	rp=$( basename $f | sed 's/_build.pkl//g' )
	for d in "${outputs_array[@]}"
	do
	  treesapp phylotu \
	  --refpkg_path $f \
	  --assign_output $d \
	  -o ${d}/phylotu_out_${rp}/
	done
done
