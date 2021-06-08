#!/bin/bash

file_array=("/mnt/sdb/Hallam_projects/Hallam_Databases/formatted/CAMI_I_HIGH/gold_standard_high_single.fasta.gz"
"../CAMI_experiments/RH_S001_merged.fq"
"../CAMI_experiments/RH_S001_forward.fq")

analysis_dir="../CAMI_experiments/"
analysis_refpkgs=${analysis_dir}/refpkgs/
n_cpu=16
outputs_array=()

for f in "${file_array[@]}"
do
  assign_output=$analysis_dir/$( basename $f | sed 's/.fq.*//g' | sed 's/.fasta.*//g' )
  treesapp assign \
  -i $f \
  -o ${assign_output} \
  --refpkg_dir ${analysis_refpkgs} \
  --num_procs $n_cpu \
  --trim_align \
  --overwrite --delete
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
