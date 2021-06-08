#!/bin/bash

treesapp_install=/home/connor/Bioinformatics/Hallam_projects/TreeSAPP/treesapp/
eval_fasta=/mnt/sdb/Hallam_projects/Hallam_Databases/raw/Protein/EggNOG_v5.0/prokaryotic_e5.proteomes.faa
annotation_class_file=/mnt/sdb/Hallam_projects/Hallam_Databases/raw/Protein/EggNOG_v5.0/EggNOG_refpkg_OG_map.tsv
analysis_root=../tax_summary/
n_cpu=16

$treesapp_install/mcc_calculator.py \
-i $eval_fasta \
--tool graftm \
--output $analysis_root/MCC_graftm_0.13.0_prok \
--refpkg_dir $analysis_root/tax_summary_refpkgs/ \
--gpkg_dir $analysis_root/GraftM_gpkgs/ \
--annot_map $annotation_class_file \
--num_procs $n_cpu \
--overwrite

$treesapp_install/mcc_calculator.py \
-i $eval_fasta \
--placement_summary aelw \
--output $analysis_root/MCC_treesapp_0.11.1_aelw_prok/ \
--refpkg_dir $analysis_root/tax_summary_refpkgs/ \
--trim_align \
--annot_map $annotation_class_file \
--num_procs $n_cpu \
--overwrite

$treesapp_install/mcc_calculator.py \
-i $eval_fasta \
--placement_summary max_lwr \
--output $analysis_root/MCC_treesapp_0.11.1_aelw_prok/ \
--refpkg_dir $analysis_root/tax_summary_refpkgs/ \
--trim_align \
--annot_map $annotation_class_file \
--num_procs $n_cpu \
--overwrite
