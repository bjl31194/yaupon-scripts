#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=pv-gatk_haptest
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=160G
#SBATCH --time=20:00:00
#SBATCH --output=GATKHC.out
#SBATCH --error=GATKHC_testerr2.err

cd /scratch/sp27971/Zoysia/ZoysiaLastPool/ZoysiaRedown/AnalysisRedown/Trimmed

ml GATK/4.2.5.0-GCCcore-8.3.0-Java-1.8

gatk HaplotypeCaller -R ZJN_r1.1.fa -I pqX.Gr.sorted.bam -stand-call-conf 30 -ERC GVCF -O pqXHap.g.vcf
