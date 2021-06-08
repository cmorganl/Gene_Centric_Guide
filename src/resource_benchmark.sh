#!/bin/bash

project_path=$( dirname $PWD )

usage="USAGE:\n
                $0 -g [genomes.csv] [-p project_path=$project_path] [-x tmp_dir=$project_path/tmp] [-r diamond_db=/mnt/nfs/sharknado/LimsData/Hallam_Databases/formatted/diamond/refseq.dmnd]"

if [ $# -eq 0 ]; then
	echo "ERROR: No arguments provided."
	echo -e $usage
	exit
elif [ $1 == "-h" ]; then
	echo -e $usage
	exit
elif [ $# -lt 1 ]; then
	echo "ERROR: Insufficient number of arguments provided."
	echo -e $usage
	exit
fi

# Set arguments
while getopts ":g:p:r:x:" opt; do
  case $opt in
  g) genomes_list="$OPTARG"
    ;;
  p) project_path="$OPTARG"
    ;;
  x) tmp_dir="$OPTARG"
    ;;
  r) refseq_db="$OPTARG"
    ;;
  \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done

if [ -z $tmp_dir ]; then
  tmp_dir="${project_path}/tmp/"
fi
if [ -z $refseq_db ]; then
  refseq_db="/mnt/nfs/sharknado/LimsData/Hallam_Databases/formatted/diamond/refseq.dmnd"
fi

gpkg_root=${project_path}/tax_summary/GraftM_gpkgs/
refpkg_root=${project_path}/tax_summary/tax_summary_refpkgs/
time_log=${project_path}/runtime_log.csv
thread_scale=(4 8 16)

if [ ! -d $gpkg_root ] | [ ! -d $refpkg_root ]; then
  echo "ERROR: unable to locate reference package directories $refpkg_root and $gpkg_root"
  exit
fi

printf "Running with the following configuration:
        Project path: $project_path
        Genomes list: $genomes_list
        TreeSAPP RefPkgs path: $refpkg_root
        GraftM RefPkgs path: $gpkg_root
        DIAMOND database path: $refseq_db
        Threads: $thread_scale
        Log file: $time_log\n"
exit

if [ ! -d $tmp_dir ]; then
	mkdir $tmp_dir
fi

if [ -f $time_log ]; then
	rm $time_log
fi
touch $time_log

echo "Sample.Name,Fasta.Length (bp),Software,Molecule,Threads,Time (mm:ss),Memory.Max (kbytes)" >>$time_log

# Download genomes
while read line
do
  taxid=$( echo "$line" | gawk -F, '{ print $2 }')
  fa_url=$( echo "$line" | gawk -F, '{ print $4 }')
  genome_file=${tmp_dir}/${taxid}.fasta.gz
  if [ -f $genome_file ]; then
    continue
  fi
  if [ $fa_url == "PrefixURL" ]; then
    continue
  fi
  echo "Downloading sequences for $( echo "$line" | gawk -F, '{ print $1 }')"
  wget -O $genome_file $fa_url 1>/dev/null 2>&1
#  gunzip -c $taxid.faa.gz | seqkit replace -p "^" -r "${taxid}." >>$proteins_fa
#  rm $taxid.faa.gz
done<$genomes_list

for fa_nuc in ${tmp_dir}/*.fasta.gz
do
  for n_procs in "${thread_scale[@]}"
  do
    fa_nuc_len=$(seqkit stats $fa_nuc | tail -n 1 | gawk '{ print $5 }')
    fa_aa=$tmp_dir/$sid/final_outputs/${sid}_ORFs.faa
	  sid=$( basename $fa_nuc | sed 's/.fasta.gz//g' )

# Run DIAMOND blastx on nucleotide sequences
    /usr/bin/time -v --output=tmp.txt bash -c "diamond blastx \
--db $refseq_db \
--query $fa_nuc \
--out $tmp_dir/diamond-refseq.tsv --outfmt 6 \
--index-chunks 1 --block-size 10 \
--taxonmap /mnt/nfs/sharknado/LimsData/Hallam_Databases/raw/Taxonomy/prot.accession2taxid \
--taxonnodes /usr/local/share/centrifuge_DBs/nodes.dmp --taxonnames /usr/local/share/centrifuge_DBs/names.dmp \
--threads $n_procs --max-target-seqs 10 --strand both"
    mmem=$(grep -w 'Maximum resident set' tmp.txt | gawk '{ print $NF }')
    rtime=$(grep -w 'Elapsed' tmp.txt | gawk '{ print $NF }')
    echo "$sid,$fa_nuc_len,DIAMOND,nuc,$n_procs,$rtime,$mmem" >>$time_log

# Run GraftM on nucleotide sequences
    for gpkg in $gpkg_root/*gpkg
    do
		  /usr/bin/time -v --output=tmp.txt bash -c "graftM graft \
--search_method hmmsearch --assignment_method pplacer \
--output_directory $tmp_dir/graftm \
--forward $fa_nuc --graftm_package $gpkg \
--verbosity 2 --threads $n_procs \
--input_sequence_type nucleotide --force"
      mmem=$(grep -w 'Maximum resident set' tmp.txt | gawk '{ print $NF }')
      rtime=$(grep -w 'Elapsed' tmp.txt | gawk '{ print $NF }')
      echo "$sid,$fa_nuc_len,GraftM,nuc,$n_procs,$rtime,$mmem" >>$time_log
	  done

# Run TreeSAPP on nucleotide sequences
	  /usr/bin/time -v --output=tmp.txt bash -c "treesapp assign \
-i $fa_nuc \
-o $tmp_dir/$sid \
--refpkg_dir $refpkg_root \
-n $n_procs \
--overwrite --trim_align"
    mmem=$(grep -w 'Maximum resident set' tmp.txt | gawk '{ print $NF }')
    rtime=$(grep -w 'Elapsed' tmp.txt | gawk '{ print $NF }')
    echo "$sid,$fa_nuc_len,TreeSAPP,nuc,$n_procs,$rtime,$mmem" >>$time_log

# Find the total number of characters in the amino acid sequences
if [ ! -f $fa_aa ]; then
  echo "ERROR: $fa_aa doesn't exist."
  exit
else
	fa_aa_len=$(seqkit stats $fa_aa | tail -n 1 | gawk '{ print $5 }')
fi

# Run TreeSAPP on amino acid sequences
	/usr/bin/time -v --output=tmp.txt bash -c "treesapp assign \
-i $fa_aa \
-o $tmp_dir/$sid_prot_treesapp \
--refpkg_dir $refpkg_root \
-n $n_procs \
--overwrite --trim_align"
	mmem=$(grep -w 'Maximum resident set' tmp.txt | gawk '{ print $NF }')
	rtime=$(grep -w 'Elapsed' tmp.txt | gawk '{ print $NF }')
	echo "$sid,$fa_aa_len,TreeSAPP,aa,$n_procs,$rtime,$mmem" >>$time_log

# Run DIAMOND blastp on amino acid sequences
		/usr/bin/time -v --output=tmp.txt bash -c "diamond blastp \
--db $refseq_db \
--query $fa_aa \
--out $tmp_dir/diamond-refseq.tsv --outfmt 6 \
--index-chunks 1 --block-size 10 \
--taxonmap /mnt/nfs/sharknado/LimsData/Hallam_Databases/raw/Taxonomy/prot.accession2taxid \
--taxonnodes /usr/local/share/centrifuge_DBs/nodes.dmp --taxonnames /usr/local/share/centrifuge_DBs/names.dmp \
--threads $n_procs --max-target-seqs 10 --strand both"
		mmem=$(grep -w 'Maximum resident set' tmp.txt | gawk '{ print $NF }')
		rtime=$(grep -w 'Elapsed' tmp.txt | gawk '{ print $NF }')
		echo "$sid,$fa_aa_len,DIAMOND,aa,$n_procs,$rtime,$mmem" >>$time_log

# Run GraftM on amino acid sequences
	for gpkg in $gpkg_root
    /usr/bin/time -v --output=tmp.txt bash -c "graftM graft \
--search_method hmmsearch --assignment_method pplacer --input_sequence_type aminoacid \
--output_directory $tmp_dir/graftm \
--forward $fa_aa \
--graftm_package $gpkg \
--verbosity 2 --threads $n_procs --force"
    mmem=$(grep -w 'Maximum resident set' tmp.txt | gawk '{ print $NF }')
    rtime=$(grep -w 'Elapsed' tmp.txt | gawk '{ print $NF }')
    echo "$sid,$fa_aa_len,GraftM,aa,$n_procs,$rtime,$mmem" >>$time_log
  done
done

if [ -f $tmp_dir ]
	then rm -r $tmp_dir
fi

