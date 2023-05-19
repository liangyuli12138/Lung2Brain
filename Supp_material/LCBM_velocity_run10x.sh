

# run velocyto for LCBM samples
# run in 202.195.187.3 


for sample in A20190305 A20190312 D0927 E0927 T_Bsc1
do
cat > /public/workspace/lily/Lung2Brain/velocity/velocyto_hg38/velocyto_${sample}_code.sh <<EOF
#!/bin/sh
#PBS -l nodes=1:ppn=10
#PBS -N velocyto_${sample}
#PBS -d /public/workspace/lily/log_qusb/
#PBS -o velocyto_${sample}_out.txt
#PBS -e velocyto_${sample}_err.txt
#PBS -l walltime=2400:00:00

echo $sample
source /bio-apps/rhel7/Miniconda3/etc/profile.d/conda.sh
conda activate velocity
module load  samtools-1.9

velocyto run10x -m /public/workspace/lily/REF/GRCh38_rmsk.gtf \
    -@ 10 --samtools-memory 2000 \
    /public/workspace/xiapeng/Brain_Tumor_sc/0raw_data/LCBM/inhouse/${sample}/3cellranger/ \
    /public/workspace/xiapeng/CellRanger/refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf

echo ${sample} "done!"

EOF

done



# get data
for sample in A20190305 A20190312 D0927 E0927 T_Bsc1
do 
    cp /public/workspace/xiapeng/Brain_Tumor_sc/0raw_data/LCBM/inhouse/${sample}/3cellranger/velocyto/3cellranger.loom \
    /public/workspace/lily/Lung2Brain/velocity/velocyto_hg38/loom/${sample}.loom
done




























