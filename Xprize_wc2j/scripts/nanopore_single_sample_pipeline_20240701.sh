#!/bin/bash
set -e 
set -u 
set -o 
##########################################################################################################
##########################################################################################################
#script for nanopore metabarcoding: chopper to blastn pipeline
#this script works for only single sample there is another script for barcodes
##########################################################################################################
##########################################################################################################
echo "starting analysis" 
echo "working_dir is" $1
echo "output_dir is" $2
cd $2
mkdir -p 12su 
mkdir -p 12sv 
mkdir -p 16s 
script="/Users/caiwang/Documents/projects/MinKNOW_output/scripts/"
Q=9
minlength=150
OTUSIM=97
luluSIM=97.1
##############################1. merge data##############################
echo "merge raw data..."
#merge_seuqnceses
gunzip $1/fastq_pass/*.fastq.gz
cat $1/fastq_pass/*.fastq > merge_all.fastq
#rm fastq_pass/$sample/FAY*.fastq 
##########2.filter qulity >Q9, mini length is minlength####################
echo "filering sequences with QC lower than" ${Q} 
echo "filering sequences shorter than" ${minlength} 
chopper -q ${Q} -l ${minlength} -t 10  -i merge_all.fastq > filtered_chopper_${Q}_${minlength}.fastq 

##############################3.trim by primer##############################
#16SMAM:CGGTTGGGGTGACCTCGGA GCTGTTATCCCTAGGGTAACT #97bp, here is about 120~
#12SV5:ACTGGGATTAGATACCCC  TAGAACAGGCTCCTCTAG -r 19:125 #80-110bp, here is about 105~
#12SU:GTGCCAGCNRCCGCGGTYANAC ATAGTRGGGTATCTAATCCYAGT -r 23:236 #207bp, here is about 211~
echo "triming sequences by primer 12sv5"
seqkit amplicon -F ACTGGGATTAGATACCCC -R TAGAACAGGCTCCTCTAG  -r 19:125 -w 0 -j 8  -m 4  filtered_chopper_${Q}_${minlength}.fastq | seqkit seq -m 80 > 12sv/filtered_seqkit_12sv.fastq 
seqkit stat 12sv/filtered_seqkit_12sv.fastq
echo "triming sequences by primer 12su"
seqkit amplicon -F GTGCCAGCNRCCGCGGTYANAC -R ATAGTRGGGTATCTAATCCYAGT  -r 23:236 -w 0 -j 8  -m 4  filtered_chopper_${Q}_${minlength}.fastq | seqkit seq -m 150 > 12su/filtered_seqkit_12su.fastq 
seqkit stat 12su/filtered_seqkit_12su.fastq
echo "triming sequences by primer 16smam2"
seqkit amplicon -F CGGTTGGGGTGACCTCGGA -R GCTGTTATCCCTAGGGTAACT  -r 20:150 -w 0 -j 8  -m 4  filtered_chopper_${Q}_${minlength}.fastq | seqkit seq -m 80 > 16s/filtered_seqkit_16s.fastq 
seqkit stat 16s/filtered_seqkit_16s.fastq
seqtk seq -a 12su/filtered_seqkit_12su.fastq > 12su/filtered_seqkit_12su.fasta
seqtk seq -a 12sv/filtered_seqkit_12sv.fastq > 12sv/filtered_seqkit_12sv.fasta
seqtk seq -a 16s/filtered_seqkit_16s.fastq > 16s/filtered_seqkit_16s.fasta

Rscript ${script}rename_single 12sv/filtered_seqkit_12sv.fasta 12sv/filtered_seqkit_rename_12sv.fasta
Rscript ${script}rename_single 12su/filtered_seqkit_12sU.fasta 12su/filtered_seqkit_rename_12su.fasta
Rscript ${script}rename_single 16s/filtered_seqkit_16s.fasta 16s/filtered_seqkit_rename_16s.fasta

##############################4.vsearch pipeline###############################
echo "vsearch pipeline to derep and de novo"
vsearch --derep_fulllength  12sv/filtered_seqkit_rename_12sv.fasta  --minuniquesize 2 --fasta_width 0 --sizeout --output 12sv/uniques_min2_12sv.fasta 
vsearch --derep_fulllength  12su/filtered_seqkit_rename_12su.fasta  --minuniquesize 2 --fasta_width 0 --sizeout --output 12su/uniques_min2_12su.fasta 
vsearch --derep_fulllength  16s/filtered_seqkit_rename_16s.fasta  --minuniquesize 2 --fasta_width 0 --sizeout --output 16s/uniques_min2_16s.fasta 

vsearch --sortbysize 12sv/uniques_min2_12sv.fasta --output 12sv/FilteredReads.forvsearch_sorted_12sv.fasta 
vsearch --sortbysize 12su/uniques_min2_12su.fasta --output 12su/FilteredReads.forvsearch_sorted_12su.fasta
vsearch --sortbysize 16s/uniques_min2_16s.fasta --output 16s/FilteredReads.forvsearch_sorted_16s.fasta  

vsearch --uchime_denovo 12sv/FilteredReads.forvsearch_sorted_12sv.fasta  --nonchimeras 12sv/FilteredReads.forvsearch_sorted_nochimeras_12sv.fasta
vsearch --uchime_denovo 12su/FilteredReads.forvsearch_sorted_12su.fasta  --nonchimeras 12su/FilteredReads.forvsearch_sorted_nochimeras_12su.fasta
vsearch --uchime_denovo 16s/FilteredReads.forvsearch_sorted_16s.fasta  --nonchimeras 16s/FilteredReads.forvsearch_sorted_nochimeras_16s.fasta

#######clusting-12su###############
echo "clusting by sumaclust with" ${OTUSIM}
gsed 's/;size/ count/' 12su/FilteredReads.forvsearch_sorted_nochimeras_12su.fasta > 12su/FilteredReads.forcluster_forsumaclust_12su.fasta
sumaclust -t .${OTUSIM} -e 12su/FilteredReads.forcluster_forsumaclust_12su.fasta > 12su/OTUs_${OTUSIM}_sumaclust_12su.fna
python ${script}tabulateSumaclust.py -i 12su/OTUs_${OTUSIM}_sumaclust_12su.fna -o 12su/OTUs_${OTUSIM}_sumaclust_forblast_12su.txt -blast
mv 12su/OTUs_${OTUSIM}_sumaclust_forblast_12su.txt.blast.txt  12su/OTUs_${OTUSIM}_sumaclust_forblast_12su.fasta
seqkit stat 12su/OTUs_${OTUSIM}_sumaclust_forblast_12su.fasta
#######clusting-12sv###############
echo "clusting by sumaclust with" ${OTUSIM}
gsed 's/;size/ count/' 12sv/FilteredReads.forvsearch_sorted_nochimeras_12sv.fasta > 12sv/FilteredReads.forcluster_forsumaclust_12sv.fasta
sumaclust -t .${OTUSIM} -e 12sv/FilteredReads.forcluster_forsumaclust_12sv.fasta > 12sv/OTUs_${OTUSIM}_sumaclust_12sv.fna
python ${script}tabulateSumaclust.py -i 12sv/OTUs_${OTUSIM}_sumaclust_12sv.fna -o 12sv/OTUs_${OTUSIM}_sumaclust_forblast_12sv.txt -blast
mv 12sv/OTUs_${OTUSIM}_sumaclust_forblast_12sv.txt.blast.txt  12sv/OTUs_${OTUSIM}_sumaclust_forblast_12sv.fasta
seqkit stat 12sv/OTUs_${OTUSIM}_sumaclust_forblast_12sv.fasta
#######clusting-16s###############
gsed 's/;size/ count/' 16s/FilteredReads.forvsearch_sorted_nochimeras_16s.fasta > 16s/FilteredReads.forcluster_forsumaclust_16s.fasta
sumaclust -t .${OTUSIM} -e 16s/FilteredReads.forcluster_forsumaclust_16s.fasta > 16s/OTUs_${OTUSIM}_sumaclust_16s.fna
python ${script}tabulateSumaclust.py -i 16s/OTUs_${OTUSIM}_sumaclust_16s.fna -o 16s/OTUs_${OTUSIM}_sumaclust_forblast_16s.txt -blast
mv 16s/OTUs_${OTUSIM}_sumaclust_forblast_16s.txt.blast.txt  16s/OTUs_${OTUSIM}_sumaclust_forblast_16s.fasta
seqkit stat 16s/OTUs_${OTUSIM}_sumaclust_forblast_16s.fasta

########lulu pipeline###############
primer=12sv
vsearch --usearch_global ${primer}/OTUs_${OTUSIM}_sumaclust_forblast_12sv.fasta --db ${primer}/OTUs_${OTUSIM}_sumaclust_forblast_12sv.fasta --strand plus --self --id .80 --iddef 1 --userout matchlist.txt --userfields query+target+id --maxaccepts 0 --query_cov .9 --maxhits 10
Rscript  ${script}lulu ${primer}/OTUs_${OTUSIM}_sumaclust_forblast_12sv.txt ${luluSIM} ${primer}
rm lulu.log_* matchlist.txt
echo "the output fasta file is OTUs_sumaclust_lulu${luluSIM}_forblast_${primer}.fasta, and otutable is table_${OTUSIM}_lulu${luluSIM}_${primer}.txt"
seqkit stat OTUs_sumaclust_lulu${luluSIM}_forblast_${primer}.fasta

primer=12su
vsearch --usearch_global ${primer}/OTUs_${OTUSIM}_sumaclust_forblast_${primer}.fasta --db ${primer}/OTUs_${OTUSIM}_sumaclust_forblast_${primer}.fasta --strand plus --self --id .80 --iddef 1 --userout matchlist.txt --userfields query+target+id --maxaccepts 0 --query_cov .9 --maxhits 10
Rscript  ${script}lulu ${primer}/OTUs_${OTUSIM}_sumaclust_forblast_${primer}.txt ${luluSIM} ${primer}
rm lulu.log_* matchlist.txt
echo "the output fasta file is OTUs_sumaclust_lulu${luluSIM}_forblast_${primer}.fasta, and otutable is table_${OTUSIM}_lulu${luluSIM}_${primer}.txt"
seqkit stat OTUs_sumaclust_lulu${luluSIM}_forblast_${primer}.fasta

primer=16s
vsearch --usearch_global ${primer}/OTUs_${OTUSIM}_sumaclust_forblast_${primer}.fasta --db ${primer}/OTUs_${OTUSIM}_sumaclust_forblast_${primer}.fasta --strand plus --self --id .80 --iddef 1 --userout matchlist.txt --userfields query+target+id --maxaccepts 0 --query_cov .9 --maxhits 10
Rscript  ${script}lulu ${primer}/OTUs_${OTUSIM}_sumaclust_forblast_${primer}.txt ${luluSIM} ${primer}
rm lulu.log_* matchlist.txt
echo "the output fasta file is OTUs_sumaclust_lulu${luluSIM}_forblast_${primer}.fasta, and otutable is table_${OTUSIM}_lulu${luluSIM}_${primer}.txt"
seqkit stat OTUs_sumaclust_lulu${luluSIM}_forblast_${primer}.fasta
########local blast###############
#echo "local blast against nt database"
#blastn -query OTUs_sumaclust_lulu${luluSIM}forblast.fasta -db /Users/caiwang/Documents/projects/MARTi/blastdb/nt -out out_min2_${OTUSIM}_sumaclust.xml -evalue 1e-5 -outfmt 5 -mt_mode 2 -num_threads 8
#subsample if need
#seqkit sample -s 100 -n 500 OTUs_${OTUSIM}_sumaclust_forblast_16s.fasta >OTUs_${OTUSIM}_sumaclust_forblast_16s_500.fasta
#seqkit stat otus.vsearch97_sampling_5000.fasta
