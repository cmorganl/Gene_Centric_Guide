#!/bin/bash

project_path=$( dirname $PWD )
default_tmp_dir="${project_path}/tmp"
default_refseq_db="/mnt/sdb/Hallam_projects/Hallam_Databases/formatted/Protein/refseq.dmnd"
default_tax_root="/mnt/sdb/Hallam_projects/Hallam_Databases/raw/Taxonomy/"

usage="
USAGE:\n
                $0 -g [genomes.csv] \
                [-p project_path=$project_path] \
                [-x tmp_dir=$default_tmp_dir] \
                [-r diamond_db=$default_refseq_db] \
                [-t taxonomy_root=$default_tax_root]
\n\n
        'project_path' is the path to the Gene_Centric_Guide directory\n
        'tmp_dir' is the path to write all intermediate files and is removed at the end of the run\n
        'diamond_db' is the path to a formatted reference database for DIAMOND\n
        'taxonomy_root' is a directory containing the prot.accession2taxid, taxdump/nodes.dmp and taxdump/names.dmp
\n\n"

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
while getopts ":g:p:r:x:t:" opt; do
  case $opt in
  g) genomes_list="$OPTARG"
    ;;
  p) project_path="$OPTARG"
    ;;
  x) tmp_dir="$OPTARG"
    ;;
  r) refseq_db="$OPTARG"
    ;;
  t) tax_root="$OPTARG"
    ;;
  \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done

if [ -z $tmp_dir ]; then
  tmp_dir=$default_tmp_dir
fi
if [ -z $refseq_db ]; then
  refseq_db=$default_refseq_db
fi
if [ -z $tax_root ]; then
  tax_root=$default_tax_root
fi

gpkg_root=${project_path}/tax_summary/GraftM_gpkgs/
refpkg_root=${project_path}/tax_summary/tax_summary_refpkgs/
cat_genomes=${tmp_dir}/genomes.fasta
time_log=${project_path}/runtime_log.csv
thread_scale=(4 8 16)
read_len=8000

if [ ! -d $gpkg_root ] || [ ! -d $refpkg_root ]; then
  echo "ERROR: unable to locate reference package directories $refpkg_root and $gpkg_root"
  exit
fi

printf "Running with the following configuration:
        Project path: $project_path
        Genomes list: $genomes_list
        TreeSAPP RefPkgs path: $refpkg_root
        GraftM RefPkgs path: $gpkg_root
        DIAMOND database path: $refseq_db
        Taxonomy root directory: $tax_root
        Simulated read lengths: $read_len
        Threads: $thread_scale
        Log file: $time_log\n\n"

##
# Clean up old outputs if necessary
if [ ! -d $tmp_dir ]; then
	mkdir $tmp_dir
fi

if [ -f $cat_genomes ]; then
  rm $cat_genomes
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
  fa_url=$( echo "$line" | gawk -F, '{ print $4 }' | sed 's/$/_genomic.fna.gz/g')
  genome_file=${tmp_dir}/${taxid}.fasta.gz
  if [ -f $genome_file ]; then
    continue
  fi
  if [ $taxid == "TaxID" ]; then
    continue
  fi
  echo "Downloading sequences for $( echo "$line" | gawk -F, '{ print $1 }')"
  wget -O $genome_file $fa_url 1>/dev/null 2>&1
  gunzip -c $genome_file | seqkit replace -p "^" -r "${taxid}." >>$cat_genomes
  rm $genome_file
done<$genomes_list

printf "Simulating metagenomes... "
factor=10
while [ $factor -lt 10000 ]
do
  len=$(echo "$factor*100" | bc)
  if [ ! -f $tmp_dir/sim_${len}.1.fna ] || [ ! -f $tmp_dir/sim_${len}.1.fna ]; then
    wgsim -1 $read_len -2 $read_len -S 2021 -N $len $cat_genomes $tmp_dir/sim_${len}.1.fq $tmp_dir/sim_${len}.2.fq 1>/dev/null 2>&1
    # Convert the FASTQ files to FASTA
    for fq in $tmp_dir/sim_${len}.?.fq
    do
      fna=$( echo $fq | sed 's/.fq/.fna/g' )
      seqkit fq2fa $fq >$fna
      rm $fq
    done
  fi
  factor=$(echo "$factor*2" | bc)
done
printf "done.\n"

for fa_nuc in ${tmp_dir}/*.fna
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
--taxonmap $tax_root/prot.accession2taxid \
--taxonnodes $tax_root/taxdump/nodes.dmp --taxonnames $tax_root/taxdump/names.dmp \
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

#if [ -d $tmp_dir ]; then
#  rm $cat_genomes
#  rm -r $tmp_dir/graftm
#  rm -r $tmp_dir/$sid_prot_treesapp $tmp_dir/$sid
#  rm $tmp_dir/diamond-refseq.tsv
#  rm $tmp_dir/sim_*.?.fna
#fi

