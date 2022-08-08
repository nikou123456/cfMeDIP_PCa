#!/bin/bash
#SBATCH --mem=10G
#SBATCH -J rm_multimapper
#SBATCH -c 1
#SBATCH -N 1
#SBATCH -o /cluster/projects/hansengroup/sujunc/methylation/5mC/log/%x-%j.out
module load samtools/1.3.1

samtools view -F 1804 -f 2 -b $1 | samtools sort -o $2
samtools flagstat $2 > $3

