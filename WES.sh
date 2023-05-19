
# use Sequenza to run WES CNV
bytlib load python-3.6.6
bytlib load R-3.6.0



pip3 install sequenza-utils

#=======================================
# run bwa 
# A20190305
/public/workspace/lily/Lung2Brain/WES/DNASeq_preproc.sh \
-1 /public/workspace/lily/Lung2Brain/WES/A20190305/2019-3-5-778-tumortissue-WXS/2019-3-5-778-tumortissue-WXS_R1.fq.gz \
-2 /public/workspace/lily/Lung2Brain/WES/A20190305/2019-3-5-778-tumortissue-WXS/2019-3-5-778-tumortissue-WXS_R2.fq.gz \
-p A20190305_T \
-t 4 \
-r /public/workspace/wangzy/ref/hg19/ucsc.hg19.fasta \
-o /public/workspace/lily/Lung2Brain/WES/res/A20190305 \
-d /public/workspace/wangzy/ref/hg19/1000G_phase1.indels.hg19.sites.vcf \
-n /public/workspace/wangzy/ref/hg19/dbsnp_138.hg19.vcf

# A20190305 ref 
/public/workspace/lily/Lung2Brain/WES/DNASeq_preproc.sh \
-1 /public/workspace/lily/Lung2Brain/WES/A20190305/2019-3-5-778-WBC-WXS/2019-3-5-778-WBC-WXS_R1.fq.gz \
-2 /public/workspace/lily/Lung2Brain/WES/A20190305/2019-3-5-778-WBC-WXS/2019-3-5-778-WBC-WXS_R2.fq.gz \
-p A20190305_C \
-t 4 \
-r /public/workspace/wangzy/ref/hg19/ucsc.hg19.fasta \
-o /public/workspace/lily/Lung2Brain/WES/res/A20190305 \
-d /public/workspace/wangzy/ref/hg19/1000G_phase1.indels.hg19.sites.vcf \
-n /public/workspace/wangzy/ref/hg19/dbsnp_138.hg19.vcf




#======================================
# 20190312
#
# A20190312
/public/workspace/lily/Lung2Brain/WES/DNASeq_preproc.sh \
-1 /public/workspace/lily/Lung2Brain/WES/A20190312/2019-3-12-783-tumortissue-WXS/2019-3-12-783-tumortissue-WXS_R1.fq.gz \
-2 /public/workspace/lily/Lung2Brain/WES/A20190312/2019-3-12-783-tumortissue-WXS/2019-3-12-783-tumortissue-WXS_R2.fq.gz \
-p A20190312_T \
-t 4 \
-r /public/workspace/wangzy/ref/hg19/ucsc.hg19.fasta \
-o /public/workspace/lily/Lung2Brain/WES/res/A20190312 \
-d /public/workspace/wangzy/ref/hg19/1000G_phase1.indels.hg19.sites.vcf \
-n /public/workspace/wangzy/ref/hg19/dbsnp_138.hg19.vcf

# A20190312 ref 
/public/workspace/lily/Lung2Brain/WES/DNASeq_preproc.sh \
-1 /public/workspace/lily/Lung2Brain/WES/A20190312/2019-3-12-783-WBC-WXS/2019-3-12-783-WBC-WXS_R1.fq.gz \
-2 /public/workspace/lily/Lung2Brain/WES/A20190312/2019-3-12-783-WBC-WXS/2019-3-12-783-WBC-WXS_R2.fq.gz \
-p A20190312_C \
-t 4 \
-r /public/workspace/wangzy/ref/hg19/ucsc.hg19.fasta \
-o /public/workspace/lily/Lung2Brain/WES/res/A20190312 \
-d /public/workspace/wangzy/ref/hg19/1000G_phase1.indels.hg19.sites.vcf \
-n /public/workspace/wangzy/ref/hg19/dbsnp_138.hg19.vcf






#======================================
# T_Bsc 
#
# T_Bsc
/public/workspace/lily/Lung2Brain/WES/DNASeq_preproc.sh \
-1 /public/workspace/lily/Lung2Brain/WES/T_Bsc/1017636T-WES/1017636T-WES_R1.fq.gz \
-2 /public/workspace/lily/Lung2Brain/WES/T_Bsc/1017636T-WES/1017636T-WES_R2.fq.gz \
-p T_Bsc_T \
-t 4 \
-r /public/workspace/wangzy/ref/hg19/ucsc.hg19.fasta \
-o /public/workspace/lily/Lung2Brain/WES/res/T_Bsc \
-d /public/workspace/wangzy/ref/hg19/1000G_phase1.indels.hg19.sites.vcf \
-n /public/workspace/wangzy/ref/hg19/dbsnp_138.hg19.vcf

# T_Bsc ref 
/public/workspace/lily/Lung2Brain/WES/DNASeq_preproc.sh \
-1 /public/workspace/lily/Lung2Brain/WES/T_Bsc/white-cell/white-cell_R1.fq.gz \
-2 /public/workspace/lily/Lung2Brain/WES/T_Bsc/white-cell/white-cell_R2.fq.gz \
-p T_Bsc_C \
-t 4 \
-r /public/workspace/wangzy/ref/hg19/ucsc.hg19.fasta \
-o /public/workspace/lily/Lung2Brain/WES/res/T_Bsc \
-d /public/workspace/wangzy/ref/hg19/1000G_phase1.indels.hg19.sites.vcf \
-n /public/workspace/wangzy/ref/hg19/dbsnp_138.hg19.vcf









#========================================================================
# prepare to run sequenza 
#
bytlib load python-3.6.6
bytlib load R-3.6.0
# pip3 install sequenza-utils
#===============================================================================================================================

#!/bin/bash
bytlib load python-3.6.6
bytlib load samtools-1.9
bytlib load tabix-0.2.6
bytlib load R-3.6.0

GENOME=/public/workspace/wangzy/ref/hg19/ucsc.hg19.fasta
#out_dir=/public/workspace/wumin/TAM_scRNA/data/WES/sequenza
out_dir=/public/workspace/lily/Lung2Brain/WES/Sequenza
#sequenza-utils gc_wiggle -w 50 --fasta $GENOME -o $out_dir/hg19.gc50Base.wig.gz
wig=$out_dir/hg19.gc50Base.wig.gz

tumor=/public/workspace/lily/Lung2Brain/WES/res/A20190312/A20190312_T.sorted.dedup.realign.recal.bam
norm=/public/workspace/lily/Lung2Brain/WES/res/A20190312/A20190312_C.sorted.dedup.realign.recal.bam
sample=A20190312
# mkdir $out_dir/$sample

sequenza-utils bam2seqz -n $norm -t $tumor --fasta $GENOME -gc $wig -o $out_dir/$sample/$sample.out.seqz.gz 1>$out_dir/$sample/$sample.sequenza.error.log 2>$out_dir/$sample/$sample.sequenza.log.log

sequenza-utils seqz_binning --seqz $out_dir/$sample/$sample.out.seqz.gz -w 50 -o $out_dir/$sample/$sample.small.out.seqz.gz 1>$out_dir/$sample/$sample.sequenza.error.log 2>$out_dir/$sample/$sample.sequenza.log.log

R CMD BATCH --no-save "--args data.file='$out_dir/$sample/$sample.small.out.seqz.gz' out.dir='$out_dir/$sample' sample.id='$sample'" /public/workspace/lily/Lung2Brain/WES/sequenza2.R /dev/null 1>$out_dir/$sample/$sample.sequenza.error.log 2>$out_dir/$sample/$sample.sequenza.log.log




#====================================================================================================
# T_Bsc sample 

#!/bin/bash
bytlib load python-3.6.6
bytlib load samtools-1.9
bytlib load tabix-0.2.6
bytlib load R-3.6.0

GENOME=/public/workspace/wangzy/ref/hg19/ucsc.hg19.fasta
#out_dir=/public/workspace/wumin/TAM_scRNA/data/WES/sequenza
out_dir=/public/workspace/lily/Lung2Brain/WES/Sequenza
#sequenza-utils gc_wiggle -w 50 --fasta $GENOME -o $out_dir/hg19.gc50Base.wig.gz
wig=$out_dir/hg19.gc50Base.wig.gz

tumor=/public/workspace/lily/Lung2Brain/WES/res/T_Bsc/T_Bsc_T.sorted.dedup.realign.recal.bam
norm=/public/workspace/lily/Lung2Brain/WES/res/T_Bsc/T_Bsc_C.sorted.dedup.realign.recal.bam
sample=T_Bsc
# mkdir $out_dir/$sample

sequenza-utils bam2seqz -n $norm -t $tumor --fasta $GENOME -gc $wig -o $out_dir/$sample/$sample.out.seqz.gz 1>$out_dir/$sample/$sample.sequenza.error.log 2>$out_dir/$sample/$sample.sequenza.log.log

sequenza-utils seqz_binning --seqz $out_dir/$sample/$sample.out.seqz.gz -w 50 -o $out_dir/$sample/$sample.small.out.seqz.gz 1>$out_dir/$sample/$sample.sequenza.error.log 2>$out_dir/$sample/$sample.sequenza.log.log

R CMD BATCH --no-save "--args data.file='$out_dir/$sample/$sample.small.out.seqz.gz' out.dir='$out_dir/$sample' sample.id='$sample'" /public/workspace/lily/Lung2Brain/WES/sequenza2.R /dev/null 1>$out_dir/$sample/$sample.sequenza.error.log 2>$out_dir/$sample/$sample.sequenza.log.log




#===================================================================================================
# A20190305


#!/bin/bash
bytlib load python-3.6.6
bytlib load samtools-1.9
bytlib load tabix-0.2.6
bytlib load R-3.6.0

GENOME=/public/workspace/wangzy/ref/hg19/ucsc.hg19.fasta
#out_dir=/public/workspace/wumin/TAM_scRNA/data/WES/sequenza
out_dir=/public/workspace/lily/Lung2Brain/WES/Sequenza
#sequenza-utils gc_wiggle -w 50 --fasta $GENOME -o $out_dir/hg19.gc50Base.wig.gz
wig=$out_dir/hg19.gc50Base.wig.gz

tumor=/public/workspace/lily/Lung2Brain/WES/res/A20190305/A20190305_T.sorted.dedup.realign.recal.bam
norm=/public/workspace/lily/Lung2Brain/WES/res/A20190305/A20190305_C.sorted.dedup.realign.recal.bam
sample=A20190305
# mkdir $out_dir/$sample

sequenza-utils bam2seqz -n $norm -t $tumor --fasta $GENOME -gc $wig -o $out_dir/$sample/$sample.out.seqz.gz 1>$out_dir/$sample/$sample.sequenza.error.log 2>$out_dir/$sample/$sample.sequenza.log.log

sequenza-utils seqz_binning --seqz $out_dir/$sample/$sample.out.seqz.gz -w 50 -o $out_dir/$sample/$sample.small.out.seqz.gz 1>$out_dir/$sample/$sample.sequenza.error.log 2>$out_dir/$sample/$sample.sequenza.log.log

R CMD BATCH --no-save "--args data.file='$out_dir/$sample/$sample.small.out.seqz.gz' out.dir='$out_dir/$sample' sample.id='$sample'" /public/workspace/lily/Lung2Brain/WES/sequenza2.R /dev/null 1>$out_dir/$sample/$sample.sequenza.error.log 2>$out_dir/$sample/$sample.sequenza.log.log
















#######################################################################################################################################################
# 2021-4-9
# sequenza result show sample 20190305 is all amplifcation ,not very good ,use varscan 
# Varscan 
#======================================================================================================================================================
#!/bin/sh

# sample A20190312
bytlib load python-3.6.6
bytlib load samtools-1.9
bytlib load VarScan-2.3.9

REF=/public/workspace/wangzy/ref/hg19/ucsc.hg19.fasta
CONTROL=/public/workspace/lily/Lung2Brain/WES/res/A20190312/A20190312_C.sorted.dedup.realign.recal.bam
CASE=/public/workspace/lily/Lung2Brain/WES/res/A20190312/A20190312_T.sorted.dedup.realign.recal.bam
outDir=/public/workspace/lily/Lung2Brain/WES/Varscan/A20190312
PREFIX="A20190312"
samtools mpileup -f $REF $CASE > $outDir/$PREFIX.case.pileup &
samtools mpileup -f $REF $CONTROL > $outDir/$PREFIX.ctrl.pileup &
wait

echo "[INFO] varscan processing ..."
VarScan copynumber \
               $outDir/$PREFIX.ctrl.pileup \
               $outDir/$PREFIX.case.pileup \
               $outDir/$PREFIX.cnv 

wait 
echo "[IFNO] sample A20190312 done"



######################################################################################################################################
# sample A20190305
bytlib load python-3.6.6
bytlib load samtools-1.9
bytlib load VarScan-2.3.9

REF=/public/workspace/wangzy/ref/hg19/ucsc.hg19.fasta
CONTROL=/public/workspace/lily/Lung2Brain/WES/res/A20190305/A20190305_C.sorted.dedup.realign.recal.bam
CASE=/public/workspace/lily/Lung2Brain/WES/res/A20190305/A20190305_T.sorted.dedup.realign.recal.bam
outDir=/public/workspace/lily/Lung2Brain/WES/Varscan/A20190305
PREFIX="A20190305"
samtools mpileup -f $REF $CASE > $outDir/$PREFIX.case.pileup &
samtools mpileup -f $REF $CONTROL > $outDir/$PREFIX.ctrl.pileup &
wait

echo "[INFO] varscan processing ..."
VarScan copynumber \
               $outDir/$PREFIX.ctrl.pileup \
               $outDir/$PREFIX.case.pileup \
               $outDir/$PREFIX.cnv 
wait 
echo "[IFNO] sample A20190305 done"


######################################################################################################################################
# sample T_Bsc
bytlib load python-3.6.6
bytlib load samtools-1.9
bytlib load VarScan-2.3.9

REF=/public/workspace/wangzy/ref/hg19/ucsc.hg19.fasta
CONTROL=/public/workspace/lily/Lung2Brain/WES/res/T_Bsc/T_Bsc_C.sorted.dedup.realign.recal.bam
CASE=/public/workspace/lily/Lung2Brain/WES/res/T_Bsc/T_Bsc_T.sorted.dedup.realign.recal.bam
outDir=/public/workspace/lily/Lung2Brain/WES/Varscan/T_Bsc
PREFIX="T_Bsc"
samtools mpileup -f $REF $CASE > $outDir/$PREFIX.case.pileup &
samtools mpileup -f $REF $CONTROL > $outDir/$PREFIX.ctrl.pileup &
wait

echo "[INFO] varscan processing ..."
VarScan copynumber \
               $outDir/$PREFIX.ctrl.pileup \
               $outDir/$PREFIX.case.pileup \
               $outDir/$PREFIX.cnv 
wait 
echo "[IFNO] sample T_Bsc done"

















######################################################################################################################################################
# CNVkit software 
# 2021-4-9
#=====================================================================================================================================================
#!/bin/sh


######################################################################################################################################################
# sample T_Bsc
bytlib load python-3.7.4
bytlib load samtools-1.9

out_dir=/public/workspace/lily/Lung2Brain/WES/CNVkit/T_Bsc
GENOME=/public/workspace/wangzy/ref/hg19/ucsc.hg19.fasta
bed=/public/workspace/wangzy/ref/Homo_sapiens.GRCh37.75.bed
CONTROL=/public/workspace/lily/Lung2Brain/WES/res/T_Bsc/T_Bsc_C.sorted.dedup.realign.recal.bam
CASE=/public/workspace/lily/Lung2Brain/WES/res/T_Bsc/T_Bsc_T.sorted.dedup.realign.recal.bam

# mkdir $out_dir

/public/workspace/wangzy/tools/cnvkit/cnvkit.py batch $CASE --normal $CONTROL --targets ${bed} --fasta $GENOME \
--drop-low-coverage --scatter --diagram --method amplicon --output-reference my_reference.cnn \
--output-dir $out_dir 1>$out_dir/CNVkit.err.log 2>$out_dir/CNVkit.log.log



######################################################################################################################################################
# sample A20190305
bytlib load python-3.7.4
bytlib load samtools-1.9

out_dir=/public/workspace/lily/Lung2Brain/WES/CNVkit/A20190305
GENOME=/public/workspace/wangzy/ref/hg19/ucsc.hg19.fasta
bed=/public/workspace/wangzy/ref/Homo_sapiens.GRCh37.75.bed
CONTROL=/public/workspace/lily/Lung2Brain/WES/res/A20190305/A20190305_C.sorted.dedup.realign.recal.bam
CASE=/public/workspace/lily/Lung2Brain/WES/res/A20190305/A20190305_T.sorted.dedup.realign.recal.bam

# mkdir $out_dir

/public/workspace/wangzy/tools/cnvkit/cnvkit.py batch $CASE --normal $CONTROL --targets ${bed} --fasta $GENOME \
--drop-low-coverage --scatter --diagram --method amplicon --output-reference my_reference.cnn \
--output-dir $out_dir 1>$out_dir/CNVkit.err.log 2>$out_dir/CNVkit.log.log





######################################################################################################################################################
# sample A20190312
bytlib load python-3.7.4
bytlib load samtools-1.9

out_dir=/public/workspace/lily/Lung2Brain/WES/CNVkit/A20190312
GENOME=/public/workspace/wangzy/ref/hg19/ucsc.hg19.fasta
bed=/public/workspace/wangzy/ref/Homo_sapiens.GRCh37.75.bed
CONTROL=/public/workspace/lily/Lung2Brain/WES/res/A20190312/A20190312_C.sorted.dedup.realign.recal.bam
CASE=/public/workspace/lily/Lung2Brain/WES/res/A20190312/A20190312_T.sorted.dedup.realign.recal.bam

# mkdir $out_dir

/public/workspace/wangzy/tools/cnvkit/cnvkit.py batch $CASE --normal $CONTROL --targets ${bed} --fasta $GENOME \
--drop-low-coverage --scatter --diagram --method amplicon --output-reference my_reference.cnn \
--output-dir $out_dir 1>$out_dir/CNVkit.err.log 2>$out_dir/CNVkit.log.log


































###########################################################################################################################################################
# BICSeq2
# 2021-4-9
#==========================================================================================================================================================

























































