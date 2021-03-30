#!/bin/bash

file_array=("/mnt/sdb/Hallam_projects/Hallam_Databases/formatted/CAMI_I_HIGH/gold_standard_high_single.fasta.gz"
"CAMI_experiments/RH_S001_merged.fq"
"CAMI_experiments/RH_S001_forward.fq")

outputs_array=()

for f in "${file_array[@]}"
do
  assign_output=$( basename $f | sed 's/.fq.*//g' | sed 's/.fasta.*//g' )
  treesapp assign \
-i $f \
-o CAMI_experiments/${assign_output} \
--refpkg_dir CAMI_experiments/refpkgs/ \
--hmm_coverage 10 \
--num_procs 16 \
--trim_align --overwrite
  outputs_array+=("CAMI_experiments/${assign_output}")
done

for f in CAMI_experiments/refpkgs/*pkl
do
	rp=$( basename $f | sed 's/_build.pkl//g' )
	for d in "${outputs_array[@]}"
	do
	  treesapp phylotu --refpkg_path $f --assign_output $d -o ${d}/phylotu_out_${rp}/
	done
done
