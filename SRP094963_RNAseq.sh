#!/usr/bin/bash


# 2021-1-18
# ###### Wed Jan 20 17:24:56 CST 2021
# use mm10
# Run RNA seq for  
# START call RNA-seq count 
# ref RNAseq pipeline 
# Run in 202.195.187.3
#======================================================================
module load STAR-2.6.1b
module load samtools-1.9
module load HTSeq

THREAD=10
REF='/public/workspace/wulx/BIT/references/REF_mm10/refdata-cellranger-mm10-3.0.0/star'
JUNK=100
THREAD=10
FA='/public/workspace/wulx/BIT/references/REF_mm10/refdata-cellranger-mm10-3.0.0/fasta/genome.fa'


#ls /public/workspace/wangzy/work/Lung/Data/RNA/02_CleanData/

for i in `ls /public/workspace/lily/ncbi/public/sra/ | cut -d "_" -f 1 | sort -u`
do

read1="/public/workspace/lily/ncbi/public/sra/${i}_1.fastq.gz"
read2="/public/workspace/lily/ncbi/public/sra/${i}_2.fastq.gz"
outDir='/public/workspace/lily/metastasis/data/verify/KRAS_KD'
PREFIX="${i}"
mkdir -p /public/workspace/lily/metastasis/data/verify/KRAS_KD/S${i}

STAR 	--runThreadN $THREAD \
     	--genomeDir $REF \
     	--readFilesIn $read1 $read2 \
     	--outFileNamePrefix $outDir/$PREFIX \
     	--readFilesCommand zcat

mkdir	$outDir/${PREFIX}_tmpREF

STAR    --runMode genomeGenerate \
    --genomeDir $outDir/${PREFIX}_tmpREF \
    --genomeFastaFiles $FA \
    --sjdbFileChrStartEnd $outDir/${PREFIX}SJ.out.tab \
    --sjdbOverhang $JUNK \
    --runThreadN $THREAD

rm $outDir/${PREFIX}Aligned.out.sam
rm $outDir/${PREFIX}Log.final.out
rm $outDir/${PREFIX}Log.out
rm $outDir/${PREFIX}Log.progress.out
rm $outDir/${PREFIX}SJ.out.tab
rm -rf $outDir/${PREFIX}_STARtmp


STAR    --runThreadN $THREAD \
    --genomeDir $outDir/${PREFIX}_tmpREF \
    --readFilesIn $read1 $read2 \
    --outFileNamePrefix $outDir/$PREFIX \
    --readFilesCommand zcat


#==================================================
# HTSeq 

STRAND='no'
MINAQUAL=10
FeatureID='gene_id'
Mode='union'
GTF='/public/workspace/wulx/BIT/references/REF_mm10/refdata-cellranger-mm10-3.0.0/genes/genes.gtf'

samtools view -bS -F 4 $outDir/$PREFIX'Aligned.out.sam' > $outDir/$PREFIX.mapped.bam
samtools sort -n -o $outDir/$PREFIX.mapped.sort.bam $outDir/$PREFIX.mapped.bam
htseq-count -f bam -s $STRAND -a $MINAQUAL -i $FeatureID  -m $Mode $outDir/$PREFIX.mapped.sort.bam $GTF > $outDir/$PREFIX.count

done
























































