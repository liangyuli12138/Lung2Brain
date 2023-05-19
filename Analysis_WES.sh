
# this program is used to analysis WES data 
# 2022-1-10
# 1. use BWA for alignment
#===========================================================================================================================================
#===========================================================================================================================================
#===========================================================================================================================================


#!/usr/bin/bash
#PBS -l nodes=1:ppn=6
#PBS -d /public/workspace/lily/
#PBS -N A20190305_T
#PBS -o /public/workspace/lily/Lung2Brain/WES/hg38/logs/A20190305_T.out
#PBS -e /public/workspace/lily/Lung2Brain/WES/hg38/logs/A20190305_T.err
#PBS -l walltime=2400:00:00
# A20190305

/public/workspace/lily/Lung2Brain/WES/hg38/DNASeq_preproc.sh \
-1 /public/workspace/lily/Lung2Brain/WES/Data/A20190305/2019-3-5-778-tumortissue-WXS/2019-3-5-778-tumortissue-WXS_R1.fq.gz \
-2 /public/workspace/lily/Lung2Brain/WES/Data/A20190305/2019-3-5-778-tumortissue-WXS/2019-3-5-778-tumortissue-WXS_R2.fq.gz \
-p A20190305_T \
-t 4 \
-r /public/workspace/wangzy/ref/hg38/gatk_source/Homo_sapiens_assembly38.fasta \
-o /public/workspace/lily/Lung2Brain/WES/hg38/BAM/A20190305 \
-d /public/workspace/wangzy/ref/hg38/gatk_source/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
-n /public/workspace/wangzy/ref/hg38/gatk_source/dbsnp_146.hg38.vcf.gz


#!/usr/bin/bash
#PBS -l nodes=1:ppn=6
#PBS -d /public/workspace/lily/
#PBS -N A20190305_C
#PBS -o /public/workspace/lily/Lung2Brain/WES/hg38/logs/A20190305_C.out
#PBS -e /public/workspace/lily/Lung2Brain/WES/hg38/logs/A20190305_C.err
#PBS -l walltime=2400:00:00
# A20190305 ref 
/public/workspace/lily/Lung2Brain/WES/hg38/DNASeq_preproc.sh \
-1 /public/workspace/lily/Lung2Brain/WES/Data/A20190305/2019-3-5-778-WBC-WXS/2019-3-5-778-WBC-WXS_R1.fq.gz \
-2 /public/workspace/lily/Lung2Brain/WES/Data/A20190305/2019-3-5-778-WBC-WXS/2019-3-5-778-WBC-WXS_R2.fq.gz \
-p A20190305_C \
-t 4 \
-r /public/workspace/wangzy/ref/hg38/gatk_source/Homo_sapiens_assembly38.fasta \
-o /public/workspace/lily/Lung2Brain/WES/hg38/BAM/A20190305 \
-d /public/workspace/wangzy/ref/hg38/gatk_source/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
-n /public/workspace/wangzy/ref/hg38/gatk_source/dbsnp_146.hg38.vcf.gz





###########################################################################################################################################
#!/usr/bin/bash
#PBS -l nodes=1:ppn=6
#PBS -d /public/workspace/lily/
#PBS -N A20190312_T
#PBS -o /public/workspace/lily/Lung2Brain/WES/hg38/logs/A20190312_T.out
#PBS -e /public/workspace/lily/Lung2Brain/WES/hg38/logs/A20190312_T.err
#PBS -l walltime=2400:00:00
# A20190312
/public/workspace/lily/Lung2Brain/WES/hg38/DNASeq_preproc.sh \
-1 /public/workspace/lily/Lung2Brain/WES/Data/A20190312/2019-3-12-783-tumortissue-WXS/2019-3-12-783-tumortissue-WXS_R1.fq.gz \
-2 /public/workspace/lily/Lung2Brain/WES/Data/A20190312/2019-3-12-783-tumortissue-WXS/2019-3-12-783-tumortissue-WXS_R2.fq.gz \
-p A20190312_T \
-t 4 \
-r /public/workspace/wangzy/ref/hg38/gatk_source/Homo_sapiens_assembly38.fasta \
-o /public/workspace/lily/Lung2Brain/WES/hg38/BAM/A20190312 \
-d /public/workspace/wangzy/ref/hg38/gatk_source/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
-n /public/workspace/wangzy/ref/hg38/gatk_source/dbsnp_146.hg38.vcf.gz


#!/usr/bin/bash
#PBS -l nodes=1:ppn=6
#PBS -d /public/workspace/lily/
#PBS -N A20190312_C
#PBS -o /public/workspace/lily/Lung2Brain/WES/hg38/logs/A20190312_C.out
#PBS -e /public/workspace/lily/Lung2Brain/WES/hg38/logs/A20190312_C.err
#PBS -l walltime=2400:00:00
# A20190312 ref 
/public/workspace/lily/Lung2Brain/WES/hg38/DNASeq_preproc.sh \
-1 /public/workspace/lily/Lung2Brain/WES/Data/A20190312/2019-3-12-783-WBC-WXS/2019-3-12-783-WBC-WXS_R1.fq.gz \
-2 /public/workspace/lily/Lung2Brain/WES/Data/A20190312/2019-3-12-783-WBC-WXS/2019-3-12-783-WBC-WXS_R2.fq.gz \
-p A20190312_C \
-t 4 \
-r /public/workspace/wangzy/ref/hg38/gatk_source/Homo_sapiens_assembly38.fasta \
-o /public/workspace/lily/Lung2Brain/WES/hg38/BAM/A20190312 \
-d /public/workspace/wangzy/ref/hg38/gatk_source/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
-n /public/workspace/wangzy/ref/hg38/gatk_source/dbsnp_146.hg38.vcf.gz



###########################################################################################################################################
# 
#!/usr/bin/bash
#PBS -l nodes=1:ppn=6
#PBS -d /public/workspace/lily/
#PBS -N T_Bsc_T
#PBS -o /public/workspace/lily/Lung2Brain/WES/hg38/logs/T_Bsc_T.out
#PBS -e /public/workspace/lily/Lung2Brain/WES/hg38/logs/T_Bsc_T.err
#PBS -l walltime=2400:00:00
# T_Bsc
/public/workspace/lily/Lung2Brain/WES/hg38/DNASeq_preproc.sh \
-1 /public/workspace/lily/Lung2Brain/WES/Data/T_Bsc/1017636T-WES/1017636T-WES_R1.fq.gz \
-2 /public/workspace/lily/Lung2Brain/WES/Data/T_Bsc/1017636T-WES/1017636T-WES_R2.fq.gz \
-p T_Bsc_T \
-t 4 \
-r /public/workspace/wangzy/ref/hg38/gatk_source/Homo_sapiens_assembly38.fasta \
-o /public/workspace/lily/Lung2Brain/WES/hg38/BAM/T_Bsc \
-d /public/workspace/wangzy/ref/hg38/gatk_source/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
-n /public/workspace/wangzy/ref/hg38/gatk_source/dbsnp_146.hg38.vcf.gz

#!/usr/bin/bash
#PBS -l nodes=1:ppn=6
#PBS -d /public/workspace/lily/
#PBS -N T_Bsc_C
#PBS -o /public/workspace/lily/Lung2Brain/WES/hg38/logs/T_Bsc_C.out
#PBS -e /public/workspace/lily/Lung2Brain/WES/hg38/logs/T_Bsc_C.err
#PBS -l walltime=2400:00:00
# T_Bsc ref 
/public/workspace/lily/Lung2Brain/WES/hg38/DNASeq_preproc.sh \
-1 /public/workspace/lily/Lung2Brain/WES/Data/T_Bsc/white-cell/white-cell_R1.fq.gz \
-2 /public/workspace/lily/Lung2Brain/WES/Data/T_Bsc/white-cell/white-cell_R2.fq.gz \
-p T_Bsc_C \
-t 4 \
-r /public/workspace/wangzy/ref/hg38/gatk_source/Homo_sapiens_assembly38.fasta \
-o /public/workspace/lily/Lung2Brain/WES/hg38/BAM/T_Bsc \
-d /public/workspace/wangzy/ref/hg38/gatk_source/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
-n /public/workspace/wangzy/ref/hg38/gatk_source/dbsnp_146.hg38.vcf.gz










#===========================================================================================================================================
#===========================================================================================================================================
#===========================================================================================================================================
# now run CNVkit prepare 
# installing by conda https://github.com/etal/cnvkit

awk '{print "chr"$0"\t0\t+"}' CCDS_exons.20180614.txt > CCDS_exons.20180614.bed
awk '{print $1"\t"$6"\t"$7"\t"$3"\t"$2}' CCDS_exons.20180614.bed > CCDS_exons.20180614.bed
awk '{print $0"\t0\t+"}' mm10.exon > mm10.exon.bed

for sampleid in A20190305 A20190312 T_Bsc ;
do 
cat>/public/workspace/lily/Lung2Brain/WES/hg38/code/CNVkit_${sampleid}_run.sh<<EOF

#!/bin/sh
#PBS -l nodes=1:ppn=4
#PBS -N CNVkit_${sampleid}_run
#PBS -d /public/workspace/lily/
#PBS -o CNVkit_${sampleid}_run.out
#PBS -e CNVkit_${sampleid}_run.err
#PBS -l walltime=2400:00:00

# module load python-3.7.4
module load samtools-1.9

out_dir=/public/workspace/lily/Lung2Brain/WES/hg38/CNVkit/${sampleid}
GENOME=/public/workspace/wangzy/ref/hg38/gatk_source/Homo_sapiens_assembly38.fasta
bed=/public/workspace/lily/Lung2Brain/WES/hg38/hg38.exon.sort.chr.bed
Tumor=/public/workspace/lily/Lung2Brain/WES/hg38/BAM/${sampleid}/${sampleid}_T.sorted.dedup.realign.recal.bam
Control=/public/workspace/lily/Lung2Brain/WES/hg38/BAM/${sampleid}/${sampleid}_C.sorted.dedup.realign.recal.bam


mkdir \$out_dir
source ~/.bashrc
conda activate cnvkit
~/.conda/envs/cnvkit/bin/cnvkit.py batch \$Tumor --normal \$Control --targets \${bed} --fasta \$GENOME \
--scatter --diagram --method amplicon --output-reference my_reference.cnn \
--output-dir \$out_dir 1>\$out_dir/CNVkit.err.log 2>\$out_dir/CNVkit.log.log
EOF
done




cnvkit.py export seg /public/workspace/lily/Lung2Brain/WES/hg38/CNVkit/A20190305/A20190305_T.sorted.dedup.realign.cns \
-o /public/workspace/lily/Lung2Brain/WES/hg38/CNVkit/A20190305/A20190305_T.seg

cnvkit.py export seg /public/workspace/lily/Lung2Brain/WES/hg38/CNVkit/T_Bsc/T_Bsc_T.sorted.dedup.realign.cns \
-o /public/workspace/lily/Lung2Brain/WES/hg38/CNVkit/T_Bsc/T_Bsc_T.seg

cnvkit.py export seg /public/workspace/lily/Lung2Brain/WES/hg38/CNVkit/A20190312/A20190312_T.sorted.dedup.realign.cns \
-o /public/workspace/lily/Lung2Brain/WES/hg38/CNVkit/A20190312/A20190312_T.seg









#===========================================================================================================================================
#===========================================================================================================================================
#===========================================================================================================================================
# now run Varscan  


# sample 
for sampleid in A20190305 A20190312 T_Bsc ;
do 
cat>/public/workspace/lily/Lung2Brain/WES/hg38/code/Varscan_${sampleid}_run.sh<<EOF

#!/bin/sh
#PBS -l nodes=1:ppn=4
#PBS -N Varscan_${sampleid}_run
#PBS -d /public/workspace/lily/
#PBS -o Varscan_${sampleid}_run.out
#PBS -e Varscan_${sampleid}_run.err
#PBS -l walltime=2400:00:00

module load python-3.6.6
module load samtools-1.9
module load VarScan.v2.3.9

mkdir -p /public/workspace/lily/Lung2Brain/WES/hg38/Varscan/${sampleid}
REF=/public/workspace/wangzy/ref/hg38/gatk_source/Homo_sapiens_assembly38.fasta 
CONTROL=/public/workspace/lily/Lung2Brain/WES/hg38/BAM/${sampleid}/${sampleid}_C.sorted.dedup.realign.recal.bam
CASE=/public/workspace/lily/Lung2Brain/WES/hg38/BAM/${sampleid}/${sampleid}_T.sorted.dedup.realign.recal.bam
outDir=/public/workspace/lily/Lung2Brain/WES/hg38/Varscan/${sampleid}
PREFIX=${sampleid}
samtools mpileup -f \$REF \$CASE > \$outDir/\$PREFIX.case.pileup &
samtools mpileup -f \$REF \$CONTROL > \$outDir/\$PREFIX.ctrl.pileup &
wait


echo "[INFO] varscan processing ..."
VarScan copynumber \
               \$outDir/\$PREFIX.ctrl.pileup \
               \$outDir/\$PREFIX.case.pileup \
               \$outDir/\$PREFIX.cnv 

wait 
echo "[IFNO]  done"

EOF
done
































