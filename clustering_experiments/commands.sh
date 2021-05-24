#!/bin/bash

usage="USAGE:\n
                $0 -g [e5.refpkg_proteome.faa] -r [path/to/RefPkgs/] [-t num_threads=8] [-o overwrite='y'|'n']"

if [ $# -eq 0 ]; then
	echo "ERROR: No arguments provided."
	echo -e $usage
	exit
elif [ $1 == "-h" ]; then
	echo -e $usage
	exit
elif [ $# -lt 4 ]; then
	echo "ERROR: Insufficient number of arguments provided."
	echo -e $usage
	exit
fi

# Set arguments
while getopts ":g:r:t:o:" opt; do
  case $opt in
  g) proteins_fa="$OPTARG"
    ;;
  r) refpkgs_repo_dir="$OPTARG"
    ;;
  t) n_threads="$OPTARG"
    ;;
  o) overwrite="$OPTARG"
    ;;
  \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done
if [ -z $n_threads ]; then
  n_threads=8
fi
if [ -z $overwrite ]; then
  overwrite="y"
fi

log_file=TreeSAPP_clustering_experiments_log.txt
refpkg_dir=clustering_refpkgs/
ts_assign_out=proteome_classified
phylotu_out=phylotu_outputs/
tmp_input=$( basename $proteins_fa | sed 's/.gz$//g' ).tmp

printf "Running with the following configuration:
        Query sequence FASTA: $proteins_fa
        RefPkgs path: $refpkgs_repo_dir
        Overwrite: $overwrite
        Threads: $n_threads
        Log file: $log_file\n"

if [ ! -d $refpkg_dir ]; then
	mkdir $refpkg_dir
fi
if [ ! -d $phylotu_out ]; then
	mkdir $phylotu_out
fi
if [ $(file --mime-type $proteins_fa | grep -c gzip ) -eq 1 ]; then
  gunzip -c $proteins_fa >$tmp_input
else
  cp $proteins_fa $tmp_input
fi

echo "Copying reference packages"
cp \
	$refpkgs_repo_dir/Methanogenesis/McrA/seed_refpkg/final_outputs/McrA_build.pkl \
	$refpkgs_repo_dir/Nitrogen_metabolism/Fixation/NifH/seed_refpkg/final_outputs/NifH_build.pkl \
	$refpkgs_repo_dir/Sulfur_metabolism/SoxY/seed_refpkg/final_outputs/SoxY_build.pkl \
	$refpkgs_repo_dir/Translation/RpoB/seed_refpkg/final_outputs/RpoB_build.pkl \
	$refpkgs_repo_dir/Translation/RecA/seed_refpkg/final_outputs/RecA_build.pkl \
	$refpkgs_repo_dir/Translation/PF01655/seed_refpkg/final_outputs/PF01655_build.pkl \
	$refpkg_dir

for r in $(seq 100 50 600)
do
	prefix=length_$r
	if [ -d $prefix ] && [ $overwrite == 'y' ]; then
	  rm -r $prefix
    mkdir $prefix
    mkdir $prefix/$phylotu_out
	elif [ ! -d $prefix ]; then
    mkdir $prefix
    mkdir $prefix/$phylotu_out
  fi

  queries=$( echo $prefix/$tmp_input | sed 's/.tmp$//g' )
  overlap=$( echo $r | awk '{ print $1/2 }' )
  echo "Creating subsequences of length $r with overlap of $overlap from $proteins_fa"
  cat $tmp_input | \
  seqkit sliding --greedy --step $overlap --window $r | \
  seqkit seq  --remove-gaps --min-len $overlap >$queries

	if [ ! -d $prefix/$ts_assign_out ]; then
    echo "Classifying the sequences"
		treesapp assign \
		--fastx_input $queries \
		--output $prefix/$ts_assign_out \
		--refpkg_dir $refpkg_dir \
		-m prot --num_procs $n_threads --delete --overwrite
	else
	  echo "Using treesapp assign classifications in $prefix/$ts_assign_out"
	fi

  echo "Clustering the classified query sequences"
	for f in $refpkg_dir/*pkl
	do
	  for rank in "species" "family" "class"
	  do
	    # Output directories
	    potu_rg_out=$prefix/$phylotu_out/$( basename $f | sed 's/_build.pkl//g' )\_phylotus_rg_$rank
	    potu_dn_psc_out=$prefix/$phylotu_out/$( basename $f | sed 's/_build.pkl//g' )\_phylotus_dn_psc_$rank
	    potu_dn_aln_out=$prefix/$phylotu_out/$( basename $f | sed 's/_build.pkl//g' )\_phylotus_dn_aln_$rank

	    # Reference-guided based on placement edges
	    if [ -d $potu_rg_out ] && [ $overwrite == 'n' ]; then
	      echo "Skipping $potu_rg_out"
	    else
        treesapp phylotu \
          --refpkg_path $f \
          --assign_output $prefix/$ts_assign_out \
          --tax_rank $rank \
          -o $potu_rg_out \
          --mode ref_guided
      fi

      # De novo clusters by recreating the phylogeny
      if [ -d $potu_dn_psc_out ] && [ $overwrite == 'n' ]; then
	      echo "Skipping $potu_dn_psc_out"
	    else
        treesapp phylotu \
          --refpkg_path $f \
          --assign_output $prefix/$ts_assign_out \
          --tax_rank $rank \
          -o $potu_dn_psc_out \
          --mode de_novo --pre_cluster psc
      fi

      # De novo clusters by recreating the phylogeny
      if [ -d $potu_dn_aln_out ] && [ $overwrite == 'n' ]; then
	      echo "Skipping $potu_dn_aln_out"
	    else
        treesapp phylotu \
          --refpkg_path $f \
          --assign_output $prefix/$ts_assign_out \
          --tax_rank $rank \
          -o $potu_dn_aln_out \
          --mode de_novo --pre_cluster align
      fi
    done
  done
done

rm $tmp_input
