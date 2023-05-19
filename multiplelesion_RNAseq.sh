#!/usr/bin/bash


# 2020-10-9 
# START call RNA-seq count 
# ref RNAseq pipeline 
#======================================================================
bytlib load STAR-2.6.1b
bytlib load samtools-1.9
bytlib load HTSeq

THREAD=8
REF='/public/workspace/wulx/references/refdata-cellranger-hg19-1.2.0/star/'
JUNK=100
THREAD=1
FA='/public/workspace/wulx/references/refdata-cellranger-hg19-1.2.0/fasta/genome.fa'


#ls /public/workspace/wangzy/work/Lung/Data/RNA/02_CleanData/

for i in `ls /public/workspace/wangzy/work/Lung/Data/RNA/02_CleanData/ | cut -d "_" -f 1 | sort -u`
do

read1="/public/workspace/wangzy/work/Lung/Data/RNA/02_CleanData/${i}_clean_R1.fq.gz"
read2="/public/workspace/wangzy/work/Lung/Data/RNA/02_CleanData/${i}_clean_R2.fq.gz"
outDir='/public/workspace/lily/Lung2Brain/multiple_lesion'
PREFIX="S${i}"
mkdir -p /public/workspace/lily/Lung2Brain/multiple_lesion/S${i}

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
GTF='/public/workspace/wulx/references/refdata-cellranger-hg19-1.2.0/genes/genes.gtf'

samtools view -bS -F 4 $outDir/$PREFIX'Aligned.out.sam' > $outDir/$PREFIX.mapped.bam
samtools sort -n -o $outDir/$PREFIX.mapped.sort.bam $outDir/$PREFIX.mapped.bam
htseq-count -f bam -s $STRAND -a $MINAQUAL -i $FeatureID  -m $Mode $outDir/$PREFIX.mapped.sort.bam $GTF > $outDir/$PREFIX.count

done




