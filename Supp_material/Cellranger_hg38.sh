
# 2021-12-30 
# this data run by xiapeng and lily run paired lung cancer samples
# this program is code recore for HG38 
#============================================================================================================================================
# this xiapeng's code 

#!/bin/sh
#PBS -l nodes=1:ppn=16
#PBS -N cellranger_T_Bsc1
#PBS -d /public/workspace/xiapeng/Brain_Tumor_sc/0data/LCBM/inhouse/T_Bsc1/
#PBS -o out.txt
#PBS -e err.txt
#PBS -l walltime=2400:00:00

#here, command for analysis
echo "hello"

# cell ranger
cellranger=/public/workspace/xiapeng/CellRanger/cellranger-3.0.2/cellranger
refdata=/public/workspace/xiapeng/CellRanger/refdata-cellranger-GRCh38-3.0.0

cd /public/workspace/xiapeng/Brain_Tumor_sc/0data/LCBM/inhouse/

id=T_Bsc1
echo $id
samples=`ls /public/workspace/xiapeng/Brain_Tumor_sc/0data/LCBM/inhouse/${id}/1.fq/ | cut -d'_' -f1 | uniq | tr "\n" , | sed 's/,$//'`
echo $samples

# 2022-6-14 
# change storge /public/workspace/xiapeng/Brain_Tumor_sc/0raw_data/LCBM/inhouse/
# /public/workspace/xiapeng/Brain_Tumor_sc/0raw_data/LCBM/inhouse/WSY/3cellranger/


cd /public/workspace/xiapeng/Brain_Tumor_sc/0data/LCBM/inhouse/${id}/
pwd
#mkdir 3cellranger

$cellranger count --id=3cellranger \
    --transcriptome=$refdata  \
    --fastqs=1.fq \
    --sample=$samples \
	--localcores=16









#==========================================================================================================================================
# for Paired Lung cancer samples
# 2021-12-30


module load cellranger-3.0.2
cellranger count --id=R21136163_hg38 \
    --sample=R21136163 \
    --fastqs=/public/workspace/lily/LCBM_pair \
    --transcriptome=/public/workspace/xiapeng/CellRanger/refdata-cellranger-GRCh38-3.0.0 \
    --chemistry=auto \
    --localcores=16 --nosecondary














#==========================================================================================================================================
# for Paired2 LCBM samples
# 2022-2-11

mv R22009109-10X-20220117-0117-1_combined_R1.fastq.gz R22009109_S1_L001_R1_001.fastq.gz
mv R22009109-10X-20220117-0117-1_combined_R2.fastq.gz R22009109_S1_L001_R2_001.fastq.gz
mv R22009109-10X-20220117-0117-2_combined_R1.fastq.gz R22009109_S1_L002_R1_001.fastq.gz
mv R22009109-10X-20220117-0117-2_combined_R2.fastq.gz R22009109_S1_L002_R2_001.fastq.gz
mv R22009109-10X-20220117-0117-3_combined_R1.fastq.gz R22009109_S1_L003_R1_001.fastq.gz
mv R22009109-10X-20220117-0117-3_combined_R2.fastq.gz R22009109_S1_L003_R2_001.fastq.gz
mv R22009109-10X-20220117-0117-4_combined_R1.fastq.gz R22009109_S1_L004_R1_001.fastq.gz
mv R22009109-10X-20220117-0117-4_combined_R2.fastq.gz R22009109_S1_L004_R2_001.fastq.gz


module load cellranger-3.0.2
cellranger count --id=R22009109_hg38 \
    --sample=R22009109 \
    --fastqs=/public/workspace/lily/LCBM_pair \
    --transcriptome=/public/workspace/xiapeng/CellRanger/refdata-cellranger-GRCh38-3.0.0 \
    --chemistry=auto \
    --localcores=16 --nosecondary




