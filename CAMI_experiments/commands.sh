#!/bin/bash

treesapp assign \
-i /mnt/sdb/Hallam_projects/Hallam_Databases/formatted/CAMI_I_HIGH/gold_standard_high_single.fasta.gz \
-o CAMI_experiments/gold_standard_high_single/ \
--refpkg_dir CAMI_experiments/refpkgs/ \
--hmm_coverage 10 \
--num_procs 16 \
--trim_align --overwrite

treesapp assign \
-i CAMI_experiments/RH_S001_merged.fq \
-o CAMI_experiments/RH_S001_merged/ \
--refpkg_dir CAMI_experiments/refpkgs/ \
--hmm_coverage 10 \
--num_procs 16 \
--trim_align --overwrite

treesapp assign \
-i CAMI_experiments/RH_S001.1.fq \
-o CAMI_experiments/RH_S001_forward/ \
--refpkg_dir CAMI_experiments/refpkgs/ \
--hmm_coverage 10 \
--num_procs 16 \
--trim_align --overwrite

for f in CAMI_experiments/refpkgs/*pkl
do
	rp=$( basename $f | sed 's/_build.pkl//g' )
	treesapp phylotu --refpkg_path $f --assign_output CAMI_experiments/gold_standard_high_single/ -o CAMI_experiments/gold_standard_high_single/phylotu_out_${rp}/
	treesapp phylotu --refpkg_path $f --assign_output CAMI_experiments/RH_S001_merged/ -o CAMI_experiments/RH_S001_merged/phylotu_out_${rp}/
	treesapp phylotu --refpkg_path $f --assign_output CAMI_experiments/RH_S001_forward/ -o CAMI_experiments/RH_S001_forward/phylotu_out_${rp}/
done

