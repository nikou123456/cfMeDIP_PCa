#!/bin/bash
#SBATCH -t 5-00:00:00
#SBATCH --mem=30G
#SBATCH -J BWA_PE_ompi
#SBATCH -p ompi
#SBATCH -c 24
#SBATCH -N 1
#SBATCH -o ./log/%x-%j.out
#SBATCH -e ./log/%x-%j.error

#===========================================
#/cluster/projects/hansengroup/sujunc/methylation/5hmC/output/macs2/TJ9.R  #fragement for SE data
#fastq
echo "##############START FASTQC#######################"
module load fastqc/0.11.5
module load bedtools/2.27.1
module load deeptools/3.2.1
module load jdk7/1.7.0_79
module load java/8
module load picard/2.6.0
module load samtools/1.3.1
module load bwa/0.7.15
module load trim_galore/0.5.0
module load R/3.5.0
data_path=/cluster/projects/tcge/cell_free_epigenomics/raw_data/TCGE-CFMe-HNSC
input=/cluster/projects/hansengroup/sujunc/wye/Health_20220421
output=/cluster/projects/hansengroup/sujunc/wye/Health_20220421/Align

cd ${input}
if [ -d "${output}/fastq" ]
then
echo "exist fastq files"
else
echo "creat dir dirctly"
mkdir ${output}/fastq
fi


if [ -d "${output}/fastqc" ]
then
echo "exist fastq files"
else
echo "creat dir dirctly"
mkdir ${output}/fastqc
fi


if [ -d "${output}/Sam_files" ]
then
echo "exist Sam files"
else
echo "creat dir dirctly"
mkdir ${output}/Sam_files
fi


if [ -d "${output}/Bam_files" ]
then
echo "exist Bam files"
else
echo "creat dir dirctly"
mkdir ${output}/Bam_files
fi


if [ -d "${output}/Sort_Bamfiles" ]
then
echo "exist Sort_Bamfiles"
else
echo "creat dir dirctly"
mkdir ${output}/Sort_Bamfiles
fi


if [ -d "${output}/Insert_size" ]
then
echo "exist Insert_size"
else
echo "creat dir dirctly"
mkdir ${output}/Insert_size
fi


cd ${data_path}
keyFile=${1}


echo "##############################################################################"
echo "##############Start analysis file ${data_path}####################################"
echo "##############start trim_galore ${readsOne}###################################"
readsOne=${keyFile}$"_R1.fastq.gz"
readsTwo=${keyFile}$"_R2.fastq.gz"


trim_galore  \
 -q 20 \
 --phred33 --stringency 3 --length 20 -e 0.1 \
 --gzip -o ${output}/fastq \
 --paired  ${data_path}/${readsOne} ${data_path}/${readsTwo} \
  >  ${output}/fastq/${keyFile}.trimgalore.log


#===========================================
#fastq
cd ${output}/fastq
readsOne=${keyFile}$"_R1_val_1.fq.gz"
readsTwo=${keyFile}$"_R2_val_2.fq.gz"
fastqc ${readsOne}
fastqc ${readsTwo}
mv ${keyFile}*.zip   ${output}/fastqc
mv ${keyFile}*.html  ${output}/fastqc
#fastq
#===========================================



#===========================================
#Alignment
#cd /cluster/projects/hansengroup/wye/genome/
#bowtie2-build  /cluster/tools/data/genomes/human/hg19/iGenomes/Sequence/WholeGenomeFasta/genome.fa /cluster/projects/ha
#This module tranlate the .fastq to .sam
echo "#####################START BWA#################"
cd ${output}/fastq
echo "###begin BWA ${keyFile}"
bwa mem   -t 24 -M  /cluster/projects/hansengroup/wye/genome/bwa/genome/genome.fa   ${readsOne}   ${readsTwo}  > ${keyFile}$".bwa.sam"
rm -f ${readsOne}
rm -f  ${readsTwo}
mv ${output}/fastq/${keyFile}$".bwa.sam"   ${output}/Sam_files
echo "###end BWA ${keyFile}"
#Alignment
#===========================================



#===========================================
#samtools
# This module translate .sam to .bam
cd ${output}/Sam_files
echo "###begin SAM2BAM ${keyFile}"
samFile=${keyFile}
bamFile=${samFile}$".bam"
dedup=${samFile}$"_s.bam"
qcbamFile=${samFile}$"_qcs.bam"
rmdupFile=${samFile}$"_rmdup.bam"
rmdupSortFile=${samFile}$"_rmdup_sort.bam"

echo "####begin ${samFile} #####"
samtools view -S -b ${keyFile}$".bwa.sam" > ${bamFile}
echo "####${bamFile}" >> ${output}/${keyFile}_flagstat_report.txt
samtools flagstat ${bamFile} >> ${output}/${keyFile}_flagstat_report.txt

samtools rmdup -s ${bamFile}  ${dedup}
echo "####${dedup}" >> ${output}/${keyFile}_flagstat_report.txt
samtools flagstat ${dedup} >> ${output}/${keyFile}_flagstat_report.txt

samtools view  -h  -f 2 -F 512 -b  ${dedup} > ${qcbamFile}
echo "####${qcbamFile}" >> ${output}/${keyFile}_flagstat_report.txt
samtools flagstat ${qcbamFile} >> ${output}/${keyFile}_flagstat_report.txt
echo  "###remove sam and bem file ${samFile} ###"
rm -f ${keyFile}$".bwa.sam"
rm -f ${dedup}

samtools rmdup -s ${qcbamFile}  ${rmdupFile}
samtools sort -@ 24 -o ${rmdupSortFile}  ${rmdupFile}
samtools index  ${rmdupSortFile}
echo  "#### ${file}#####"  >> ${output_path}/remove_dup_flagstat_report_rest.txt
samtools flagstat ${rmdupSortFile} >> ${output_path}/remove_dup_flagstat_report_rest.txt
rm ${rmdupFile}
mv  ${rmdupSortFile}  ${output}/Sort_Bamfiles
mv  ${rmdupSortFile}.bai   ${output}/Sort_Bamfiles


cd ${output}/Sort_Bamfiles


#======================================================
#=======================================================
#insert size
cd ${output}/Sort_Bamfiles
echo "########STAR anaysis insert size PE ${keyFile}"
java -jar $picard_dir/picard.jar CollectInsertSizeMetrics \
      I=${output}/Sort_Bamfiles/${rmdupSortFile} \
      O=${output}/Insert_size/${keyFile}_insert_size_metrics.txt \
      H=${output}/Insert_size/${keyFile}_insert_size_histogram.pdf \
      M=0.5
echo "#########end picard ${keyFile}########"



if [ -d "${output}/Wiggle" ]
then
echo "exist Wiggle"
else
echo "creat dir dirctly"
mkdir ${output}/Wiggle
fi


if [ -d "${output}/QC_Rdata" ]
then
echo "exist QC_Rdata"
else
echo "creat dir dirctly"
mkdir ${output}/QC_Rdata
fi

#wiggle_path=/cluster/projects/hansengroup/wye/5mc/Paired_End/BWA/Sort_Bamfiles_redo/Wiggle
Rscript  /cluster/projects/hansengroup/wye/5mc/Paired_End/BWA/Sort_Bamfiles_redo/Main_MEDIPS_Count.R  ${rmdupSortFile} ${output}/Sort_Bamfiles ${output}/Wiggle
#qc_path=/cluster/projects/hansengroup/wye/5mc/Paired_End/BWA/Sort_Bamfiles_redo/QC_Rdata
Rscript /cluster/projects/hansengroup/wye/5mc/Paired_End/BWA/Sort_Bamfiles_redo/MEDIPS_QC_report.R   ${rmdupSortFile}   ${output}/Sort_Bamfiles   ${output}/QC_Rdata



