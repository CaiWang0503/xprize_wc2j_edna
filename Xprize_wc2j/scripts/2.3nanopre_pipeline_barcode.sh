#!/bin/bash
set -e
set -u
set -o
##########################################################################################################
##########################################################################################################
### script for nanopore metabarcoding: chopper to blastn pipeline ########################################
##########################################################################################################
##########################################################################################################
echo "starting analysis" 
#DAME="/usr/local/bin/DAMe/bin/"
#set up the working_dir where the MinKNOW raw data should be
working_dir="/Users/caiwang/MinKNOW/data/zoo_3primers"
HOMEFOLDER="/Users/caiwang/Documents/projects/MinKNOW_output"
cd ${HOMEFOLDER}
mkdir zoo_3primers
HOMEFOLDER="/Users/caiwang/Documents/projects/MinKNOW_output/zoo_3primers"
cd ${HOMEFOLDER}
script="/Users/caiwang/Documents/projects/MinKNOW_output/scripts/"

echo "Home folder is" ${HOMEFOLDER}
Q=9
minlength=150
OTUSIM=97
luluSIM=97
########## merge data ##########
echo "merge raw data..."
#merge_seuqnceses, the samplelist.txt should be saved on the homefolder, has barcode info.
sample_info=samplelist.txt 
sample_names=($(cut -f 1 "${sample_info}" | uniq)) 

for sample in "${sample_names[@]}" 
do
gunzip ${working_dir}/fastq_pass/$sample/*.gz
cat ${working_dir}/fastq_pass/$sample/*.fastq > all_$sample.fastq
#rm fastq_pass/$sample/FAY*.fastq 
done
cat all_*.fastq > merge_all.fastq
rm all_*.fastq
#1.filer qulity >Q9, mini length is 100,  "-c, --contam <CONTAM>, Filter contaminants against a fasta"
echo "filering sequences with QC lower than" ${Q} 
echo "filering sequences shorter than" ${minlength} 
chopper -q ${Q} -l ${minlength} -t 10  -i merge_all.fastq > filtered_chopper_${Q}_${minlength}.fastq 

########## trim by primer ##########
#16SMAM:CGGTTGGGGTGACCTCGGA GCTGTTATCCCTAGGGTAACT # 97bp
#12SV5:ACTGGGATTAGATACCCC  TAGAACAGGCTCCTCTAG  #80 and 104 bp
#12SU:GTGCCAGCNRCCGCGGTYANAC ATAGTRGGGTATCTAATCCYAGT #207bp, here is about 211~
primer=12sv
mkdir -p ${primer}
echo "triming sequences by primer 12sv"
seqkit amplicon -F ACTGGGATTAGATACCCC -R TAGAACAGGCTCCTCTAG -r 19:-19 -w 0 -j 8  -m 4  filtered_chopper_${Q}_${minlength}.fastq | seqkit seq -m 80 -M 110 > ${primer}/filtered_seqkit_${primer}.fastq

seqkit stat ${primer}/filtered_seqkit_${primer}.fastq
seqtk seq -a ${primer}/filtered_seqkit_${primer}.fastq > ${primer}/filtered_seqkit_${primer}.fasta
Rscript ${script}rename_single ${primer}/filtered_seqkit_${primer}.fasta ${primer}/filtered_seqkit_rename_${primer}.fasta

primer=12su 
mkdir -p ${primer}
echo "triming sequences by primer 12su"
seqkit amplicon -F GTGCCAGCNRCCGCGGTYANAC -R ATAGTRGGGTATCTAATCCYAGT  -r 23:-24 -w 0 -j 8  -m 4  filtered_chopper_${Q}_${minlength}.fastq | seqkit seq -m 200 -M 211 > ${primer}/filtered_seqkit_${primer}.fastq 
seqkit stat ${primer}/filtered_seqkit_${primer}.fastq
seqtk seq -a ${primer}/filtered_seqkit_${primer}.fastq > ${primer}/filtered_seqkit_${primer}.fasta
Rscript ${script}rename_single ${primer}/filtered_seqkit_${primer}.fasta ${primer}/filtered_seqkit_rename_${primer}.fasta

primer=16s 
mkdir -p ${primer}
echo "triming sequences by primer 16smam2"
seqkit amplicon -F CGGTTGGGGTGACCTCGGA -R GCTGTTATCCCTAGGGTAACT  -r 20:-22 -w 0 -j 8  -m 4  filtered_chopper_${Q}_${minlength}.fastq | seqkit seq -m 90 -M 110> $primer/filtered_seqkit_${primer}.fastq 
seqkit stat ${primer}/filtered_seqkit_${primer}.fastq
seqtk seq -a ${primer}/filtered_seqkit_${primer}.fastq > ${primer}/filtered_seqkit_${primer}.fasta
Rscript ${script}rename_single ${primer}/filtered_seqkit_${primer}.fasta ${primer}/filtered_seqkit_rename_${primer}.fasta
rm ${primer}/filtered_seqkit_${primer}.fasta

########## vsearch pipeline ##############
input_vsearch=filtered_seqkit_rename_${primer}.fasta
cd ${primer}
echo "vsearch pipeline" 
# --minuniquesize 2 only run for big dataset
vsearch --threads 10 --derep_fulllength  ${input_vsearch} --fasta_width 0 --sizeout --output uniques_min2.fasta 
vsearch --sortbysize uniques_min2.fasta --output FilteredReads.forvsearch_sorted.fasta 
vsearch  --threads 10 --uchime_denovo FilteredReads.forvsearch_sorted.fasta --nonchimeras FilteredReads.forvsearch_sorted_nochimeras.fasta 
rm uniques_min2.fasta FilteredReads.forvsearch_sorted.fasta 

########## clusting ##############
echo "clusting by sumaclust with" ${OTUSIM}
gsed 's/;size/ count/' FilteredReads.forvsearch_sorted_nochimeras.fasta > FilteredReads.forcluster_forsumaclust.fasta
sumaclust -t .${OTUSIM} -e FilteredReads.forcluster_forsumaclust.fasta > OTUs_${OTUSIM}_sumaclust.fna 
python ${script}tabulateSumaclust.py -i OTUs_${OTUSIM}_sumaclust.fna -o OTUs_${OTUSIM}_sumaclust_forblast_${primer}.txt -blast
mv OTUs_${OTUSIM}_sumaclust_forblast_${primer}.txt.blast.txt  OTUs_${OTUSIM}_sumaclust_forblast_${primer}.fasta
seqkit stat OTUs_${OTUSIM}_sumaclust_forblast_${primer}.fasta
cp OTUs_${OTUSIM}_sumaclust_forblast_${primer}.txt /Users/caiwang/Documents/projects/xprize/R_code/data/OTUs_table_${primer}.txt

rm FilteredReads.forcluster_forsumaclust.fasta
########## lulu pipeline ##############
#lulu DOESN'T WORK FOR ONLY ONE SAMPLE!! skip lulu step.
vsearch --usearch_global OTUs_${OTUSIM}_sumaclust_forblast_${primer}.fasta --db OTUs_${OTUSIM}_sumaclust_forblast_${primer}.fasta --strand plus --self --id .80 --iddef 1 --userout matchlist.txt --userfields query+target+id --maxaccepts 0 --query_cov .9 --maxhits 10
Rscript  ${script}lulu OTUs_${OTUSIM}_sumaclust_forblast_${primer}.txt ${luluSIM} ${primer}
rm lulu.log_* matchlist.txt
echo "the output fasta file is OTUs_sumaclust_lulu${luluSIM}_forblast_${primer}.fasta, and otutable is table_${OTUSIM}_lulu${luluSIM}_${primer}.txt"
seqkit stat OTUs_sumaclust_lulu${luluSIM}_forblast_${primer}.fasta
cp  table_${OTUSIM}_lulu${luluSIM}_${primer}.txt /Users/caiwang/Documents/projects/xprize/R_code/data/OTUs_table_${primer}.txt

mkdir middle_files
mv OTUs_${OTUSIM}_sumaclust_forblast_${primer}.fasta middle_files/OTUs_${OTUSIM}_sumaclust_forblast_${primer}.fasta
mv OTUs_${OTUSIM}_sumaclust_forblast_${primer}.txt middle_files/OTUs_${OTUSIM}_sumaclust_forblast_${primer}.txt

########################### species identification #################################
########################### local blast ###############
echo "local blast against nt database"
date; blastn -query OTUs_sumaclust_lulu97.1_forblast_12sv_2.fasta -db /Users/caiwang/src/BLCA/db/nt -out out_min2_${OTUSIM}_sumaclust.xml -evalue 1e-5 -outfmt 5 -mt_mode 2 -num_threads 12 ;date

#14:30.23, 5,653 seqs, 12 threads
########################### BLCA ###############
mkdir BLCA_output
cd BLCA_output
mkdir defaul
cd defaul
date; python /Users/caiwang/src/BLCA/2.blca_main.py -a muscle -i ${HOMEFOLDER}/${primer}/OTUs_sumaclust_lulu${luluSIM}_forblast_${primer}.fasta -r /Users/caiwang/src/BLCA/db_${primer}/${primer}.ACC.taxonomy.txt -q /Users/caiwang/src/BLCA/db_${primer}/${primer}.fasta -p 12 -j 20  -c 0.3 -o OTUs_sumaclust_lulu97_defaul.txt ;date
#for single species
date; python /Users/caiwang/src/BLCA/2.blca_main.py -a muscle -i ${HOMEFOLDER}/${primer}/OTUs_97_sumaclust_forblast_${primer}.fasta -r /Users/caiwang/src/BLCA/db_${primer}/${primer}.ACC.taxonomy.txt -q /Users/caiwang/src/BLCA/db_${primer}/${primer}.fasta -p 12 -j 20  -c 0.3 -o OTUs_sumaclust_lulu97_defaul.txt ;date
#replace by hand
gsed 's/;/\t/g'  OTUs_sumaclust_lulu97_defaul.txt > OTUs_sumaclust_lulu97_defaul_forR.txt
gsed 's/:/\t/g'  OTUs_sumaclust_lulu97_defaul_forR.txt > OTUs_sumaclust_lulu97_defaul_forR2.txt
mv OTUs_sumaclust_lulu97_defaul_forR2.txt /Users/caiwang/Documents/projects/xprize/R_code/data/OTUs_sumaclust_lulu97_defaul_forR.txt

HOMEFOLDER="/Users/caiwang/Documents/projects/MinKNOW_output/zoo_3primers"
primer=16s
luluSIM=97
BLCA_output=BLCA_output
cd ${HOMEFOLDER}/${primer}/${BLCA_output}
mkdir 97
cd 97
date; python /Users/caiwang/src/BLCA/2.blca_main.py -a muscle -i ${HOMEFOLDER}/${primer}/OTUs_sumaclust_lulu${luluSIM}_forblast_${primer}.fasta -r /Users/caiwang/src/BLCA/db_${primer}/${primer}.ACC.taxonomy.txt -q /Users/caiwang/src/BLCA/db_${primer}/${primer}.fasta -p 12 -j 20  -c 0.3 --iset 97 -o OTUs_sumaclust_lulu97_97.txt ;date
#for single species
date; python /Users/caiwang/src/BLCA/2.blca_main.py -a muscle -i ${HOMEFOLDER}/${primer}/OTUs_97_sumaclust_forblast_${primer}.fasta -r /Users/caiwang/src/BLCA/db_${primer}/${primer}.ACC.taxonomy.txt -q /Users/caiwang/src/BLCA/db_${primer}/${primer}.fasta -p 12 -j 20  -c 0.3 --iset 97 -o OTUs_sumaclust_lulu97_97.txt ;date
gsed 's/;/\t/g' OTUs_sumaclust_lulu97_97.txt > OTUs_sumaclust_lulu97_97_forR.txt
gsed 's/:/\t/g' OTUs_sumaclust_lulu97_97_forR.txt > OTUs_sumaclust_lulu97_97_forR2.txt
mv OTUs_sumaclust_lulu97_97_forR2.txt /Users/caiwang/Documents/projects/xprize/R_code/data/OTUs_sumaclust_lulu97_97_forR.txt

HOMEFOLDER="/Users/caiwang/Documents/projects/MinKNOW_output/zoo_3primers"
primer=16s
luluSIM=97
BLCA_output=BLCA_output
cd ${HOMEFOLDER}/${primer}/${BLCA_output}
mkdir 99
cd 99
date; python /Users/caiwang/src/BLCA/2.blca_main.py -a muscle -i ${HOMEFOLDER}/${primer}/OTUs_sumaclust_lulu${luluSIM}_forblast_${primer}.fasta -r /Users/caiwang/src/BLCA/db_${primer}/${primer}.ACC.taxonomy.txt -q /Users/caiwang/src/BLCA/db_${primer}/${primer}.fasta  -p 12 -j 20  -c 0.3 --iset 99 -o OTUs_sumaclust_lulu97_99.txt ;date
#for single species
date; python /Users/caiwang/src/BLCA/2.blca_main.py -a muscle -i ${HOMEFOLDER}/${primer}/OTUs_97_sumaclust_forblast_${primer}.fasta -r /Users/caiwang/src/BLCA/db_${primer}/${primer}.ACC.taxonomy.txt -q /Users/caiwang/src/BLCA/db_${primer}/${primer}.fasta -p 12 -j 20  -c 0.3 --iset 99 -o OTUs_sumaclust_lulu97_99.txt ;date

gsed 's/;/\t/g' OTUs_sumaclust_lulu97_99.txt > OTUs_sumaclust_lulu97_99_forR.txt
gsed 's/:/\t/g' OTUs_sumaclust_lulu97_99_forR.txt > OTUs_sumaclust_lulu97_99_forR2.txt
mv OTUs_sumaclust_lulu97_99_forR2.txt /Users/caiwang/Documents/projects/xprize/R_code/data/OTUs_sumaclust_lulu97_99_forR.txt

#running merge three result, it's 1.merge_BLCAoutputs.Rmd, need to move three blca output to Rcode/data folder
#running collapses otus, it's 2.collases.R, need to set start_samples and end_samples manually
#make metacoder figure it's 3.figure.R

################### end ######################
#Rscript $Rscriptpath/owi_collapse -s 16 -e 104 -t 0.50 -i spongetank_SWARM3_OTU_LULU_20220220.csv
