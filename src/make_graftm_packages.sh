#!/bin/bash

refpkg_path=../tax_summary/tax_summary_refpkgs/
gpkg_dir=../tax_summary/GraftM_gpkgs/

if [ ! -d $refpkg_path ]; then
	echo "ERROR: Unable to find path to directory containing reference packages for taxonomic summary analysis."
	exit
fi

if [ ! -d $gpkg_dir ]; then
	mkdir $gpkg_dir
fi

treesapp info

for rp in $refpkg_path/*pkl
do
	prefix=$( basename $rp | sed 's/_build.pkl//g' )
	tmp_dir=$gpkg_dir/$prefix/
	treesapp package view msa -r $rp | seqkit seq -g >$tmp_dir/$prefix.fasta
	treesapp package view lineage_ids -r $rp | gawk -F"\t" '{ OFS="\t"; print $1,$3 }' | sed "s/\t/_$prefix\t/g" >$tmp_dir/lineage_ids.tsv
	graftM create \
		--sequences $tmp_dir/$prefix.fasta \
		--taxonomy $tmp_dir/lineage_ids.tsv \
		--output $gpkg_dir/$prefix.gpkg
done

