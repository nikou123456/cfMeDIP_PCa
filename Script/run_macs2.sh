#!/bin/bash
#SBATCH --mem=20G
#SBATCH -J run_macs2
#SBATCH -c 1
#SBATCH -N 1
#SBATCH -o /cluster/projects/hansengroup/sujunc/methylation/5mC/log/%x-%j.out
module load MACS/2.2.5
module load R/3.5.0
module load ucsctools/378
module load bwtool/1.0
###takes in 3 arguments from commandline
### 1, input fq director
### 2, name stem
### 3, output master directory
macs2 predictd -i $1 -g hs -m 5 50 --rfile $2.R --outdir $3
cd $3
d=`grep 'altd ' $2.R|cut -f 2 -d '('|sed 's/)//g'|grep ,`
if [[ $d ]]
then
d=`echo $d | rev|cut -f 1 -d ' '|rev`
echo $d
else
d=`grep 'altd ' $2.R|cut -f 2 -d '('|sed 's/)//g'`
echo $d
fi

macs2 callpeak -t $1 -g hs -n $2 --bw 250 --mfold 5 50  --extsize $d --seed 1 --fix-bimodal --qvalue 0.05 --SPMR -B
bedSort $2"_treat_pileup.bdg" $2"_treat_pileup.bdg"
bedClip $2"_treat_pileup.bdg" /cluster/home/sujunc/chensj/usr/lib/sizes.ucsc.hg19.fasta.fai $2.bdg.c
bedGraphToBigWig $2.bdg.c /cluster/home/sujunc/chensj/usr/lib/sizes.ucsc.hg19.fasta.fai $2.bw
bedSort $2"_control_lambda.bdg" $2"_control_lambda.bdg"
bedClip $2"_control_lambda.bdg" /cluster/home/sujunc/chensj/usr/lib/sizes.ucsc.hg19.fasta.fai $2.bdg.c
bedGraphToBigWig $2.bdg.c /cluster/home/sujunc/chensj/usr/lib/sizes.ucsc.hg19.fasta.fai $2"_ctl.bw"
if [ -f $2.bw -a -f $2"_ctl.bw" ]
        then
        rm $2"_treat_pileup.bdg" $2"_control_lambda.bdg"
fi

