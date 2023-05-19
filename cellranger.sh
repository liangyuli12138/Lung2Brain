#!/usr/bin/bash

#=====================
# 2020-11-9
#====================
#================================
# run cellranger fo early stage LUAD samples
#==============================================
# run in pang ;use cellranger-3.1.0

for i in `ls /public/workspace/lily/Lung2Brain/Data/oncogene_data | grep gz | cut -d "_" -f 1 | sort -u `
do
#statements
bytlib load cellranger-3.1.0

        cellranger count --id=$i \
                --sample=${i} \
                --fastqs=/public/workspace/lily/Lung2Brain/Data/oncogene_data/ \
                --transcriptome=/public/workspace/lily/REF/INDEX-hg19/refdata-cellranger-hg19-1.2.0 \
                --chemistry=SC3Pv2 \
                --localcores=10 --localmem=60 --nosecondary

done








#================================================================
# lung brain multiple 
#================================================================
cellranger count --id=D0927 \
    --sample=D0927 \
    --fastqs=/public/workspace/lily/Mutiple_LB/D0927 \
    --transcriptome=/public/workspace/lily/REF/INDEX-hg19/refdata-cellranger-hg19-1.2.0 \
    --chemistry=SC3Pv3 \
    --localcores=10 --localmem=60 --nosecondary



cellranger count --id=E0927 \
    --sample=E0927 \
    --fastqs=/public/workspace/lily/Mutiple_LB/E0927 \
    --transcriptome=/public/workspace/lily/REF/INDEX-hg19/refdata-cellranger-hg19-1.2.0 \
    --chemistry=SC3Pv3 \
    --localcores=10 --localmem=60 --nosecondary





cellranger count --id=E20927 \
    --sample=E20927 \
    --fastqs=/public/workspace/lily/Mutiple_LB/E20927 \
    --transcriptome=/public/workspace/lily/REF/INDEX-hg19/refdata-cellranger-hg19-1.2.0 \
    --chemistry=SC3Pv3 \
    --localcores=10 --localmem=60 --nosecondary








# 2021-11-20 
# a new LCBM 
#==========================================================================================================================================
mv R21125541-10X-WSY-WSY-1_combined_R1.fastq.gz R21125541_S1_L001_R1_001.fastq.gz
mv R21125541-10X-WSY-WSY-1_combined_R2.fastq.gz R21125541_S1_L001_R2_001.fastq.gz
mv R21125541-10X-WSY-WSY-2_combined_R1.fastq.gz R21125541_S1_L002_R1_001.fastq.gz
mv R21125541-10X-WSY-WSY-2_combined_R2.fastq.gz R21125541_S1_L002_R2_001.fastq.gz
mv R21125541-10X-WSY-WSY-3_combined_R1.fastq.gz R21125541_S1_L003_R1_001.fastq.gz
mv R21125541-10X-WSY-WSY-3_combined_R2.fastq.gz R21125541_S1_L003_R2_001.fastq.gz
mv R21125541-10X-WSY-WSY-4_combined_R1.fastq.gz R21125541_S1_L004_R1_001.fastq.gz
mv R21125541-10X-WSY-WSY-4_combined_R2.fastq.gz R21125541_S1_L004_R2_001.fastq.gz



module load cellranger-3.0.2
cellranger count --id=R21125541 \
    --sample=R21125541 \
    --fastqs=/public/workspace/lily/LCBM_pair \
    --transcriptome=/public/workspace/lily/REF/refdata-cellranger-hg19-1.2.0 \
    --chemistry=auto \
    --localcores=10 --localmem=60 --nosecondary









mv R21136163-10X-WSY-lung-lung1_combined_R1.fastq.gz R21136163_S1_L001_R1_001.fastq.gz
mv R21136163-10X-WSY-lung-lung1_combined_R2.fastq.gz R21136163_S1_L001_R2_001.fastq.gz
mv R21136163-10X-WSY-lung-lung2_combined_R1.fastq.gz R21136163_S1_L002_R1_001.fastq.gz
mv R21136163-10X-WSY-lung-lung2_combined_R2.fastq.gz R21136163_S1_L002_R2_001.fastq.gz
mv R21136163-10X-WSY-lung-lung3_combined_R1.fastq.gz R21136163_S1_L003_R1_001.fastq.gz
mv R21136163-10X-WSY-lung-lung3_combined_R2.fastq.gz R21136163_S1_L003_R2_001.fastq.gz
mv R21136163-10X-WSY-lung-lung4_combined_R1.fastq.gz R21136163_S1_L004_R1_001.fastq.gz
mv R21136163-10X-WSY-lung-lung4_combined_R2.fastq.gz R21136163_S1_L004_R2_001.fastq.gz



module load cellranger-3.0.2
cellranger count --id=R21136163 \
    --sample=R21136163 \
    --fastqs=/public/workspace/lily/LCBM_pair \
    --transcriptome=/public/workspace/lily/REF/refdata-cellranger-hg19-1.2.0 \
    --chemistry=auto \
    --localcores=10 --localmem=60 --nosecondary






































