#!/bin/bash
#SBATCH --mem=10G
#SBATCH -J sum_sig
#SBATCH -c 1
#SBATCH -N 1
#SBATCH -o /cluster/projects/hansengroup/sujunc/methylation/5hmC/log/%x-%j.out
module load bwtool/1.0
## 1st arg, input bed region file
## 2nd arg, input bw signal file
## 3rd arg, output txt file name
bwtool summary $1 $2 -header -with-sum $3

