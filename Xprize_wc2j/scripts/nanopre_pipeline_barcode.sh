#!/bin/bash
set -e
set -u
set -o
##########################################################################################################
##########################################################################################################
# script for nanopore metabarcoding: chopper to blastn pipeline
##########################################################################################################
###########################################################################################################
echo "starting analysis" 
#DAME="/usr/local/bin/DAMe/bin/"
HOMEFOLDER="/Users/caiwang/Documents/projects/MinKNOW_output/bigtree"
script="/Users/caiwang/Documents/projects/MinKNOW_output/scripts/"
cd ${HOMEFOLDER}
echo "Home folder is" ${HOMEFOLDER}
Q=9
minlength=150
OTUSIM=97
luluSIM=97.1
##########merge data##########
echo "merge raw data..."
#merge_seuqnceses
sample_info=samplelist.txt 
sample_names=($(cut -f 1 "${sample_info}" | uniq)) 

for sample in "${sample_names[@]}" 
do
gunzip fastq_pass/$sample/*.gz
cat fastq_pass/$sample/*.fastq > all_$sample.fastq
#rm fastq_pass/$sample/FAY*.fastq 
done
cat all_*.fastq > merge_all.fastq
rm all_*.fastq
#1.filer qulity >Q10, mini length is 100,  "-c, --contam <CONTAM>, Filter contaminants against a fasta"
echo "filering sequences with QC lower than" ${Q} 
echo "filering sequences shorter than" ${minlength} 
chopper -q ${Q} -l ${minlength} -t 10  -i merge_all.fastq > filtered_chopper_${Q}_${minlength}.fastq 

##########trim by primer##########
#16SMAM:CGGTTGGGGTGACCTCGGA GCTGTTATCCCTAGGGTAACT #97bp, here is about 120~
#12SV5:ACTGGGATTAGATACCCC  TAGAACAGGCTCCTCTAG -r 19:125 #80-110bp, here is about 105~
#12SU:GTGCCAGCNRCCGCGGTYANAC ATAGTRGGGTATCTAATCCYAGT -r 23:236 #207bp, here is about 211~
primer=12sv
mkdir -p ${primer}
echo "triming sequences by primer 12sv"
seqkit amplicon -F ACTGGGATTAGATACCCC -R TAGAACAGGCTCCTCTAG  -w 0 -j 8  -m 4  filtered_chopper_${Q}_${minlength}.fastq | seqkit seq -m 110 -M 141 > ${primer}/filtered_seqkit_${primer}.fastq
#seqkit subseq -r 19:-19| head ${primer}/filtered_seqkit_${primer}.fastq > ${primer}/filtered_seqkit_trim_${primer}.fastq
seqkit stat ${primer}/filtered_seqkit_${primer}.fastq
seqtk seq -a ${primer}/filtered_seqkit_${primer}.fastq > ${primer}/filtered_seqkit_${primer}.fasta
Rscript ${script}rename ${primer}/filtered_seqkit_${primer}.fasta ${primer}/filtered_seqkit_rename_${primer}.fasta

primer=12su 
mkdir -p ${primer}
echo "triming sequences by primer 12su"
seqkit amplicon -F GTGCCAGCNRCCGCGGTYANAC -R ATAGTRGGGTATCTAATCCYAGT  -r 23:236 -w 0 -j 8  -m 4  filtered_chopper_${Q}_${minlength}.fastq | seqkit seq -m 150 > $primer/filtered_seqkit_${primer}.fastq 
seqkit stat ${primer}/filtered_seqkit_${primer}.fastq
seqtk seq -a ${primer}/filtered_seqkit_${primer}.fastq > ${primer}/filtered_seqkit_${primer}.fasta
Rscript ${script}rename_single ${primer}/filtered_seqkit_${primer}.fasta ${primer}/filtered_seqkit_rename_${primer}.fasta

primer=16s 
mkdir -p ${primer}
echo "triming sequences by primer 16smam2"
seqkit amplicon -F CGGTTGGGGTGACCTCGGA -R GCTGTTATCCCTAGGGTAACT  -r 20:150 -w 0 -j 8  -m 4  filtered_chopper_${Q}_${minlength}.fastq | seqkit seq -m 80 > $primer/filtered_seqkit_${primer}.fastq 
seqkit stat ${primer}/filtered_seqkit_${primer}.fastq
seqtk seq -a ${primer}/filtered_seqkit_${primer}.fastq > ${primer}/filtered_seqkit_${primer}.fasta
Rscript ${script}rename_single ${primer}/filtered_seqkit_${primer}.fasta ${primer}/filtered_seqkit_rename_${primer}.fasta

##########vsearch pipeline##############
primer=12sv
input_vsearch=filtered_seqkit_rename_${primer}.fasta
cd ${primer}
echo "vsearch pipeline" 
# --minuniquesize 2
vsearch --derep_fulllength  ${input_vsearch}  --fasta_width 0 --sizeout --output uniques_min2.fasta -–threads 10
#vsearch --derep_fulllength  filtered_seqkit_rename.fasta  --fasta_width 0 --sizeout --output uniques_all.fasta 
vsearch --sortbysize uniques_min2.fasta --output FilteredReads.forvsearch_sorted.fasta -–threads 10
vsearch --uchime_denovo FilteredReads.forvsearch_sorted.fasta --nonchimeras FilteredReads.forvsearch_sorted_nochimeras.fasta -–threads 10
#rm  filtered_chopper_${Q}_${minlength}.fastq filtered_seqkit.fastq 
#rm uniques_min2.fasta FilteredReads.forvsearch_sorted.fasta filtered_seqkit.fasta
#######clusting###############
echo "clusting by sumaclust with" ${OTUSIM}
gsed 's/;size/ count/' FilteredReads.forvsearch_sorted_nochimeras.fasta > FilteredReads.forcluster_forsumaclust.fasta
sumaclust -t .${OTUSIM} -e FilteredReads.forcluster_forsumaclust.fasta > OTUs_${OTUSIM}_sumaclust.fna 
python ${script}tabulateSumaclust.py -i OTUs_${OTUSIM}_sumaclust.fna -o OTUs_${OTUSIM}_sumaclust_forblast_${primer}.txt -blast
mv OTUs_${OTUSIM}_sumaclust_forblast_${primer}.txt.blast.txt  OTUs_${OTUSIM}_sumaclust_forblast_${primer}.fasta
seqkit stat OTUs_${OTUSIM}_sumaclust_forblast_${primer}.fasta

########lulu pipeline###############
vsearch --usearch_global OTUs_${OTUSIM}_sumaclust_forblast.fasta --db OTUs_${OTUSIM}_sumaclust_forblast.fasta --strand plus --self --id .80 --iddef 1 --userout matchlist.txt --userfields query+target+id --maxaccepts 0 --query_cov .9 --maxhits 10
Rscript  ${script}lulu OTUs_${OTUSIM}_sumaclust_forblast_${primer}.txt ${luluSIM} ${primer}
rm lulu.log_* matchlist.txt
echo "the output fasta file is OTUs_sumaclust_lulu${luluSIM}_forblast_${primer}.fasta, and otutable is table_${OTUSIM}_lulu${luluSIM}_${primer}.txt"
seqkit stat OTUs_sumaclust_lulu${luluSIM}_forblast_${primer}.fasta
########local blast###############
echo "local blast against nt database"
blastn -query OTUs_sumaclust_lulu97.1_forblast_12sv.fasta -db /Users/caiwang/src/BLCA/db/nt -out out_min2_${OTUSIM}_sumaclust.xml -evalue 1e-5 -outfmt 5 -mt_mode 2 -num_threads 10

