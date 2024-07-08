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

script="/Users/caiwang/Documents/projects/MinKNOW_output/scripts/"
Q=9
minlength=200
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
primer_info=primerlist.txt 
primer_names=($(cut -f 1 "${primer_info}" | uniq)) 

for primer in "${primer_names[@]}" 
do
#16SMAM:CGGTTGGGGTGACCTCGGA GCTGTTATCCCTAGGGTAACT #97bp, here is about 120~
#12SV5:ACTGGGATTAGATACCCC  TAGAACAGGCTCCTCTAG -r 19:125 #80-110bp, here is about 105~
#12SU:GTGCCAGCNRCCGCGGTYANAC ATAGTRGGGTATCTAATCCYAGT -r 23:236 #207bp, here is about 211~
if [ $primer == "12sv" ]
then
mkdir -p 12sv 
echo "triming sequences by primer 12sv5"
seqkit amplicon -F ACTGGGATTAGATACCCC -R TAGAACAGGCTCCTCTAG  -r 19:125 -w 0 -j 8  -m 4  filtered_chopper_${Q}_${minlength}.fastq | seqkit seq -m 80 > $primer/filtered_seqkit_12sv.fastq 
seqkit stat $primer/filtered_seqkit_$primer.fastq
seqtk seq -a $primer/filtered_seqkit_$primer.fastq > $primer/filtered_seqkit_$primer.fasta
Rscript ${script}rename_single $primer/filtered_seqkit_$primer.fasta $primer/filtered_seqkit_rename_$primer.fasta
elif [ $primer == "12su" ] 
then
mkdir -p 12su 
echo "triming sequences by primer 12su"
seqkit amplicon -F GTGCCAGCNRCCGCGGTYANAC -R ATAGTRGGGTATCTAATCCYAGT  -r 23:236 -w 0 -j 8  -m 4  filtered_chopper_${Q}_${minlength}.fastq | seqkit seq -m 150 > $primer/filtered_seqkit_12su.fastq 
seqkit stat $primer/filtered_seqkit_$primer.fastq
seqtk seq -a $primer/filtered_seqkit_$primer.fastq > $primer/filtered_seqkit_$primer.fasta
Rscript ${script}rename_single $primer/filtered_seqkit_$primer.fasta $primer/filtered_seqkit_rename_$primer.fasta
else [ $primer == "16s" ] 
mkdir -p 16s 
echo "triming sequences by primer 16smam2"
seqkit amplicon -F CGGTTGGGGTGACCTCGGA -R GCTGTTATCCCTAGGGTAACT  -r 20:150 -w 0 -j 8  -m 4  filtered_chopper_${Q}_${minlength}.fastq | seqkit seq -m 80 > $primer/filtered_seqkit_16s.fastq 
seqkit stat $primer/filtered_seqkit_$primer.fastq
seqtk seq -a $primer/filtered_seqkit_$primer.fastq > $primer/filtered_seqkit_$primer.fasta
Rscript ${script}rename_single $primer/filtered_seqkit_$primer.fasta $primer/filtered_seqkit_rename_$primer.fasta
fi
##############################4.vsearch pipeline###############################
echo "vsearch pipeline to derep and de novo"--minuniquesize 2
vsearch --derep_fulllength  $primer/filtered_seqkit_rename_$primer.fasta  --fasta_width 0 --sizeout --output $primer/uniques_min2_$primer.fasta 
vsearch --sortbysize $primer/uniques_min2_$primer.fasta --output $primer/FilteredReads.forvsearch_sorted_$primer.fasta  
vsearch --uchime_denovo $primer/FilteredReads.forvsearch_sorted_$primer.fasta  --nonchimeras $primer/FilteredReads.forvsearch_sorted_nochimeras_$primer.fasta

#######clusting###############
gsed 's/;size/ count/' $primer/FilteredReads.forvsearch_sorted_nochimeras_$primer.fasta > $primer/FilteredReads.forcluster_forsumaclust_$primer.fasta
sumaclust -t .${OTUSIM} -e $primer/FilteredReads.forcluster_forsumaclust_$primer.fasta > $primer/OTUs_${OTUSIM}_sumaclust_$primer.fna
python ${script}tabulateSumaclust.py -i $primer/OTUs_${OTUSIM}_sumaclust_$primer.fna -o $primer/OTUs_${OTUSIM}_sumaclust_forblast_$primer.txt -blast
mv $primer/OTUs_${OTUSIM}_sumaclust_forblast_$primer.txt.blast.txt  $primer/OTUs_${OTUSIM}_sumaclust_forblast_$primer.fasta
seqkit stat $primer/OTUs_${OTUSIM}_sumaclust_forblast_$primer.fasta

########lulu pipeline###############
#vsearch --usearch_global $primer/OTUs_${OTUSIM}_sumaclust_forblast_$primer.fasta --db $primer/OTUs_${OTUSIM}_sumaclust_forblast_$primer.fasta --strand plus --self --id .80 --iddef 1 --userout matchlist.txt --userfields query+target+id --maxaccepts 0 --query_cov .9 --maxhits 10
#Rscript  ${script}lulu $primer/OTUs_${OTUSIM}_sumaclust_forblast_$primer.txt ${luluSIM} $primer
#rm lulu.log_* matchlist.txt
#echo "the output fasta file is OTUs_sumaclust_lulu${luluSIM}_forblast_${primer}.fasta, and otutable is table_${OTUSIM}_lulu${luluSIM}_$primer.txt"
#seqkit stat OTUs_sumaclust_lulu${luluSIM}_forblast_$primer.fasta

done
########local blast###############
#echo "local blast against nt database"
#blastn -query OTUs_sumaclust_lulu${luluSIM}forblast.fasta -db /Users/caiwang/Documents/projects/MARTi/blastdb/nt -out out_min2_${OTUSIM}_sumaclust.xml -evalue 1e-5 -outfmt 5 -mt_mode 2 -num_threads 8
#subsample if need
#seqkit sample -s 100 -n 500 OTUs_${OTUSIM}_sumaclust_forblast_16s.fasta >OTUs_${OTUSIM}_sumaclust_forblast_16s_500.fasta
#seqkit stat otus.vsearch97_sampling_5000.fasta
