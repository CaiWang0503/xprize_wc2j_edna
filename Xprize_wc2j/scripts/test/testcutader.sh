primer=12sv
mkdir -p ${primer}
cutadapt filtered_chopper_${Q}_${minlength}.fastq  --rc --match-read-wildcards -j 0 -m 80 --maximum-length 110 -e 0.2 --max-n 0 --discard-untrimmed -a ACTGGGATTAGATACCCC...TAGAACAGGCTCCTCTAG  -n 1 -o ${primer}/filtered_cutadapt_${primer}.fastq
seqkit stat ${primer}/filtered_cutadapt_${primer}.fastq
seqtk seq -a ${primer}/filtered_cutadapt_${primer}.fastq > ${primer}/filtered_cutadapt_${primer}.fasta
Rscript ${script}rename ${primer}/filtered_cutadapt_${primer}.fasta ${primer}/filtered_seqkit_rename_${primer}.fasta

primer=12sv
input_vsearch=filtered_seqkit_rename_${primer}.fasta
cd ${primer}
echo "vsearch pipeline" 
vsearch --derep_fulllength  ${input_vsearch}  --fasta_width 0 --sizeout --output uniques_min2.fasta 
#vsearch --derep_fulllength  filtered_seqkit_rename.fasta  --fasta_width 0 --sizeout --output uniques_all.fasta 
vsearch --sortbysize uniques_min2.fasta --output FilteredReads.forvsearch_sorted.fasta 
vsearch --uchime_denovo FilteredReads.forvsearch_sorted.fasta --nonchimeras FilteredReads.forvsearch_sorted_nochimeras.fasta
echo "clusting by sumaclust with" ${OTUSIM}
gsed 's/;size/ count/' FilteredReads.forvsearch_sorted_nochimeras.fasta > FilteredReads.forcluster_forsumaclust.fasta
sumaclust -t .${OTUSIM} -e FilteredReads.forcluster_forsumaclust.fasta > OTUs_${OTUSIM}_sumaclust.fna
python ${script}tabulateSumaclust.py -i OTUs_${OTUSIM}_sumaclust.fna -o OTUs_${OTUSIM}_sumaclust_forblast.txt -blast
mv OTUs_${OTUSIM}_sumaclust_forblast.txt.blast.txt  OTUs_${OTUSIM}_sumaclust_forblast.fasta
seqkit stat OTUs_${OTUSIM}_sumaclust_forblast.fasta
