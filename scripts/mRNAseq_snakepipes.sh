#bin/bash
#########################################################################
#   Script to run mRNA-seq analysis with SnakePipes
#########################################################################

conda init
conda activate snakePipes

#mRNAseq -i /home/cvillazon/RNA-seq_Fastq -o /home/cvillazon/mRNA_SNAKEPipes_out\
#    --DAG --trim --fastqc --rMats\
#    --sampleSheet /home/cvillazon/samplesheet.tsv GRCm39


ncRNAseq -i /home/cvillazon/mRNA_SNAKEPipes_out/filtered_bam -o /home/cvillazon/ncRNA_SNAKEPipes_out\
    --fromBAM --DAG --trim \
    --sampleSheet /home/cvillazon/samplesheet.tsv GRCm39


conda init
conda deactivate