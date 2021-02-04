#!/bin/bash

usage="USAGE:\n
                $0 genomes.csv path/to/RefPkgs/ [threads=8]"

if [ $# -eq 0 ]; then
	echo "ERROR: No arguments provided."
	echo -e $usage
	exit
elif [ $1 == "-h" ]; then
	echo -e $usage
	exit
elif [ $# -lt 2 ]; then
	echo "ERROR: Insufficient number of arguments provided."
	echo $usage
fi

# Set arguments
genomes_list=$1
refpkgs_repo_dir=$2
if [ -z $3 ]; then
  n_threads=8
else
  n_threads=$3
fi

proteins_fa=proteome.fasta
refpkg_dir=clustering_refpkgs/
ts_assign_out=proteome_classified
phylotu_out=phylotu_outputs/

if [ -f $proteins_fa ]; then
	rm $proteins_fa
fi
if [ ! -d $refpkg_dir ]; then
	mkdir $refpkg_dir
fi
if [ ! -d $phylotu_out ]; then
	mkdir $phylotu_out
fi

# Download genomes
while read line
do
	echo "Downloading sequences for $( echo "$line" | gawk -F, '{ print $1 }')"
	taxid=$( echo "$line" | gawk -F, '{ print $2 }')
	faa_url=$( echo "$line" | gawk -F, '{ print $4 }')
	wget -O $taxid.faa.gz $faa_url 1>/dev/null 2>&1
	gunzip -c $taxid.faa.gz >>$proteins_fa
	rm $taxid.faa.gz
done<$genomes_list

echo "Copying reference packages"
cp \
	$refpkgs_repo_dir/Methanogenesis/McrA/seed_refpkg/final_outputs/McrA_build.pkl \
	$refpkgs_repo_dir/Nitrogen_metabolism/Fixation/NifH/seed_refpkg/final_outputs/NifH_build.pkl \
	$refpkgs_repo_dir/Sulfur_metabolism/SoxY/seed_refpkg/final_outputs/SoxY_build.pkl \
	$refpkgs_repo_dir/Translation/RpoB/seed_refpkg/final_outputs/RpoB_build.pkl \
	$refpkgs_repo_dir/Translation/PF01655/seed_refpkg/final_outputs/PF01655_build.pkl \
	$refpkg_dir

for r in 50 100 200 400 "full"
do
	prefix=length_$r
	if [ -d $prefix ]; then
	       rm -r $prefix
        fi
	mkdir $prefix
	mkdir $prefix/$phylotu_out

	if [ $r == "full" ]; then
	  queries=$proteins_fa
	else
	  queries=$prefix/$proteins_fa
	  overlap=$( echo $r | awk '{ print $1/2 }' )
	  echo "Creating subsequences of length $r with overlap of $overlap from $proteins_fa"
	  cat $proteins_fa | seqkit sliding --greedy --step $overlap --window $r | seqkit seq --min-len 30 >$queries
	fi

  # TODO: transit across different taxonomic ranks
	# Classify the sequences
	if [ ! -d $prefix/$ts_assign_out ]; then 
		treesapp assign \
		--fastx_input $queries \
		--output $prefix/$ts_assign_out \
		--refpkg_dir $refpkg_dir \
		-m prot --num_procs $n_threads --delete --overwrite
	fi

  # Cluster the classified query sequences
	for f in $refpkg_dir/*pkl
	do
	  # Reference-guided based on placement edges
		treesapp phylotu \
			--refpkg_path $f \
			--assign_output $prefix/$ts_assign_out \
			-o $prefix/$phylotu_out/$( basename $f | sed 's/_build.pkl//g' )\_phylotus_rg_s
		# De novo clusters by recreating the phylogeny
		treesapp phylotu \
			--refpkg_path $f \
			--assign_output $prefix/$ts_assign_out \
			-o $prefix/$phylotu_out/$( basename $f | sed 's/_build.pkl//g' )\_phylotus_dn_s \
			--mode de_novo
	done
done

