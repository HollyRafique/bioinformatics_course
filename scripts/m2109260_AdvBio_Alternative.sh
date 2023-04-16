#!/bin/bash
# ALTERNATIVE VERSION
#   uses bowtie2 as the aligner instead of bwa mem
#
# Installation of bowtie2 required using the following command:
# conda install -c bioconda bowtie2

# requires the reference genome to be stored in DATA_DIR/reference
# wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
# This script requires that your reference index is created in advance by running
# bwa index ~/ngs_course/dnaseq/data/reference/hg19.fa.gz
#

################## Step 2.1 INIT STEPS #################

echo "** Step 2.1 - init"

# ILLUMINACLIP Value is Machine Specific so need to set this before running the script
nextera_path="/home/ubuntu/anaconda3/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic-0.39-2/adapters/NexteraPE-PE.fa:2:30:10"

## Set up directory variables to make it easier to change the location in the future
FILE_NAME="NGS"

PROJECT_DIR="/home/ubuntu/bioinformatics_course"
RESULTS_DIR="${PROJECT_DIR}/results"
DATA_DIR="${PROJECT_DIR}/data"
ALIGN_DIR="${DATA_DIR}/aligned_data"
STATS_DIR="${DATA_DIR}/stats"


##check the folders exist - if they don't then create them
if [ ! -d "${ALIGN_DIR}" ]; then
  mkdir -p $ALIGN_DIR
fi

if [ ! -d "${STATS_DIR}" ]; then
  mkdir -p $STATS_DIR
fi

echo "download data"
cd $DATA_DIR

if [ ! -d "untrimmed_fastq" ]; then
  mkdir -p untrimmed_fastq
fi

if [ ! -d "trimmed_fastq" ]; then
  mkdir -p trimmed_fastq
fi

#setting the names of the files into variables as it's easy to make a mistake
R1=NGS0001.R1.fastq
R2=NGS0001.R2.fastq

##Download the data
#commented out as it is unnecessary to do this every time we run the script
#wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/annotation.bed

#cd untrimmed_fastq
#wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/NGS0001.R1.fastq.qz
#wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/NGS0001.R2.fastq.qz
#zcat NGS0001.R1.fastq.qz > $R1
#zcat NGS0001.R1.fastq.qz > $R2
cd


###################  Section 2.2 PRE-ALIGNMENT QC ################
echo "** Section 2.2 Pre-Alignment QC*"
echo -e "****fastqc 1\n\n"
## FASTQC
fastqc -t 4 ${DATA_DIR}/untrimmed_fastq/*.fastq.qz \
	> $STATS_DIR/untrimmed_fastqc.log



## Trimming
echo -e "****trimming\n\n"
#nextera_path set at the start of the script in section 2.1

trimmomatic PE -threads 4 -phred33 ${DATA_DIR}/untrimmed_fastq/${R1} ${DATA_DIR}/untrimmed_fastq/${R2} -baseout ${DATA_DIR}/trimmed_fastq/trimmed_data ILLUMINACLIP:$nextera_path TRAILING:25 MINLEN:50

chmod +x ${DATA_DIR}/trimmed_fastq/trimmed_data*



## FastQC
echo -e "****fastqc 2\n\n"
fastqc -t 4 ${DATA_DIR}/trimmed_fastq/trimmed_data_1P \
	${DATA_DIR}/trimmed_fastq/trimmed_data_2P \
	> $STATS_DIR/trimmed_fastqc.log


###only create the trimmed reads folder in the Results Dir
### if it does not already exist (required to avoid generating an error)

if [ ! -d "${RESULTS_DIR}/fastqc_trimmed_reads" ]; then
  mkdir -p ${RESULTS_DIR}/fastqc_trimmed_reads 
fi
## move the trimmed reads into the results folder
mv ${DATA_DIR}/trimmed_fastq/*fastqc* ${RESULTS_DIR}/fastqc_trimmed_reads/

echo -e "\n\n"



################## Section 2.3 ALIGNMENT ##################
#####      ALTERNATIVE    ######

echo "** Step 2.3 Alignment *"
## 0) Reference index must be created before running this script to avoid timing out

## 1) align with read group info

#ID 11V6WR1.111.D1375ACXX.1
#SM m2109260 (Sample Name: using my student ID) 
#PL ILLUMINA
#LB NGS-AdvBio (Library Name)
#DT 2017-01-01 (Date run was produced - making an assumption here as we don't have the data)
#PU D1375ACXX.1 (Platform Unit)

read_grp_info="@RG\tID:11V6WR1.111.D1375ACXX.1\tSM:m2109260\tPL:ILLUMINA\tLB:NGS-AdvBio\tDT:2017-01-01\tPU:D1375ACXX.1"


############### ORIGINAL HERE  ########################
#echo "**** bwa mem"
#bwa mem -t 4 -v 1 -R $read_grp_info -I 250,50 ${DATA_DIR}/reference/hg19.fa.gz \
#       ${DATA_DIR}/trimmed_fastq/trimmed_data_1P ${DATA_DIR}/trimmed_fastq/trimmed_data_2P \
#       > ${ALIGN_DIR}/${FILE_NAME}_aligned.sam

############## FINISHED ORIGINAL   ########################




###############  ALTERNATIVE HERE  ########################
#
## using bowtie2 as an alternative alignment tool

echo "**** ALTERNATIVE ALIGNER: bowtie2 "

echo "**** build bowtie2 index"
# first create the necessary bowtie index files for the reference genome:
# Creating the bowtie index takes over an hour to run

#bowtie2-build ${DATA_DIR}/reference/hg19.fa ${DATA_DIR}/reference/hg19

echo "**** run bowtie2 aligner"
# -x is the bowtie reference prefix as built in the previous step
# -1 and -2 are the trimmed, paired reads
# -S is the output file
bowtie2 -x ${DATA_DIR}/reference/hg19 \
        -1 ${DATA_DIR}/trimmed_fastq/trimmed_data_1P -2 ${DATA_DIR}/trimmed_fastq/trimmed_data_2P \
        -S ${ALIGN_DIR}/${FILE_NAME}_aligned.sam


############## FINISHED ALTERNATIVE ######################


echo "**** convert sam to bam"
# Convert SAM to BAM
samtools view -h -b ${ALIGN_DIR}/${FILE_NAME}_aligned.sam > ${ALIGN_DIR}/${FILE_NAME}_aligned.bam 


#cleanup - sam file no longer needed
rm ${ALIGN_DIR}/${FILE_NAME}_aligned.sam  


echo "**** sort and index bam"
# sort bam and create index (.bai file)
samtools sort ${ALIGN_DIR}/${FILE_NAME}_aligned.bam  > ${ALIGN_DIR}/${FILE_NAME}_sorted.bam 

samtools index ${ALIGN_DIR}/${FILE_NAME}_sorted.bam



## 2) Duplicate Marking
echo -e "**** duplicate marking\n\n"
picard MarkDuplicates I=${ALIGN_DIR}/${FILE_NAME}_sorted.bam O=${ALIGN_DIR}/${FILE_NAME}_sorted_marked.bam M=$STATS_DIR/marked_dup_metrics.txt > $STATS_DIR/mark_duplicates.log

samtools index ${ALIGN_DIR}/${FILE_NAME}_sorted_marked.bam
echo -e "\n\n"



## 3) Post-Alignment Read Filtering
echo "**** post-alignment read filtering"
samtools view -F 1796  -q 20 -o $ALIGN_DIR/${FILE_NAME}_sorted_filtered.bam $ALIGN_DIR/${FILE_NAME}_sorted_marked.bam

samtools index ${ALIGN_DIR}/${FILE_NAME}_sorted_filtered.bam




## 4) Generate Alignment Statistics
echo "**** Generate Alignment Statistics"


echo "***** 0 samtools stats"
samtools stats $ALIGN_DIR/${FILE_NAME}_sorted_filtered.bam > $STATS_DIR/filtered.stats


#flagstats - compare marked to filtered
#INPUT: output from picard mark duplicates then samtools view read filtering
echo "****** 1 - flagstats"
samtools flagstat $ALIGN_DIR/${FILE_NAME}_sorted_marked.bam > $STATS_DIR/marked.flagstat
samtools flagstat $ALIGN_DIR/${FILE_NAME}_sorted_filtered.bam > $STATS_DIR/marked_and_filtered.flagstat


echo "****** 2-idxstats"
#Generate alignment statistics per chromosome (samtools)
samtools idxstats $ALIGN_DIR/${FILE_NAME}_sorted_filtered.bam > $STATS_DIR/marked_and_filtered.idxstats

echo "****** 3 insert size"
echo -e "\n\n"
#Determine the distribution of insert sizes - picard
picard CollectInsertSizeMetrics I=$ALIGN_DIR/${FILE_NAME}_sorted_filtered.bam \
        O=$STATS_DIR/marked_and_filtered_insertsize.stats H=$STATS_DIR/marked_and_filtered_insertsize.pdf \
	> $STATS_DIR/insertsize.log
echo -e "\n\n"



echo "****** 4 Depth of Coverage"
# Depth of Coverage - bedtools

# need to use bedtools coverage -sorted option to reduce how much memory the process uses

# bam file is sorted by hg19 chromosome order (for some reason puts chrX between chr7 and chr8)
# but annotation.bed is sorted lexographically
# have to resort the annotation.bed to the same order as the bam files 
# in order to use the -sorted option of bedtools coverage 
bedtools sort -i $DATA_DIR/annotation.bed -faidx ${DATA_DIR}/reference/hg19_chrom.bed > $DATA_DIR/sorted_annotation.bed

bedtools coverage -sorted -g ${DATA_DIR}/reference/hg19_chrom.bed -a $DATA_DIR/sorted_annotation.bed -b $ALIGN_DIR/${FILE_NAME}_sorted_filtered.bam | sort -k5,5n > $STATS_DIR/DOC_sorted.bed


# DOC Summary Statistics

#need the values from column 5 as this is Depth of Coverage for each feature  in the latest version of bedtools coverage
#for each row 'count':
#	add the value in col 1 for that row into an array at position 'count' and increment count by 1
#	add the value to sum and the square of the value to sum
#	if the value is smaller than the min then set min to the value
#	if the value is greater than the max then set max to the value
#finally calc mean by dividing sum by the number of rows (could use NR instead of count for this)
#	calc standard deviation from the sum of squares and the mean
# 	calc median by finding the middle value 
#	print all the values into a new summary file

echo "DEPTH OF COVERAGE STATS: "
awk 'BEGIN{count=0; sum=0; sumsq=0; min=0; max=0}
     {a[count++] = $5; sum+=$5; sumsq+=$5^2; if(count==1){min=$5; max=$5}; if($5<min){min=$5}; if($5>max){max=$5}}
     END{ mean=sum/count; stdev=sqrt((sumsq/count)-(mean^2)); 
     print "Mean: " mean; print "Stdev: " stdev; 
     print "Median: " a[int(count/2)]; print "Min: " min; print "Max: " max}
     ' $STATS_DIR/DOC_sorted.bed | tee $STATS_DIR/DOC_summmary.stats 



#cleanup - essential to minimise storage requirements
rm $ALIGN_DIR/${FILE_NAME}_sorted.*
rm $ALIGN_DIR/${FILE_NAME}_aligned.*
rm $ALIGN_DIR/${FILE_NAME}_sorted_marked*




################### Section 2.4 VARIANT CALLING #################

echo "** Step 2.4 Variant Calling *"

## 1) Call Variants using Freebayes restricting the analysis to the regions in the bed file provided 
#index it with samtools faidx, call variants with Freebayes,
#compress the resulting variant file (VCF) and index the VCF with tabix:

#decompress the reference genome
zcat $DATA_DIR/reference/hg19.fa.gz > $DATA_DIR/reference/hg19.fa

echo "**** starting samtools faidx"
samtools faidx $DATA_DIR/reference/hg19.fa

echo "**** starting freebayes"
freebayes --bam $ALIGN_DIR/${FILE_NAME}_sorted_filtered.bam \
	--fasta-reference $DATA_DIR/reference/hg19.fa \
	--vcf $RESULTS_DIR/${FILE_NAME}.vcf

echo "**** starting bgzip"
#zip up the vcf to reduce storage
bgzip $RESULTS_DIR/${FILE_NAME}.vcf

echo "**** starting tabix"
#index the variant file
tabix -p vcf $RESULTS_DIR/${FILE_NAME}.vcf.gz



## 2) Quality Filter Variants using your choice of filters 

#We will apply the following filters to the Freebayes output:
#QUAL > 1: removes extremely low quality sites 
#QUAL / AO > 10 : additional contribution of each observation should be 10 log units (~ Q10 per read)
#SAF > 0 & SAR > 0 : reads on both forward and reverse strands 
#RPR > 1 & RPL > 1 : at least two reads “balanced” to each side of the site

echo "**** starting vcffilter"
vcffilter -f "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1"  \
       $RESULTS_DIR/${FILE_NAME}.vcf.gz > $RESULTS_DIR/${FILE_NAME}_filter.tmp.vcf


echo "bedtools intersect"
#select only the variants that exist in the annotation file 
#as they are what we are interested in
bedtools intersect -header -wa -a $RESULTS_DIR/${FILE_NAME}_filter.tmp.vcf \
       -b $DATA_DIR/annotation.bed > $RESULTS_DIR/${FILE_NAME}_filtered.vcf

bgzip $RESULTS_DIR/${FILE_NAME}_filtered.vcf
tabix -p vcf $RESULTS_DIR/${FILE_NAME}_filtered.vcf.gz

#cleanup
rm $RESULTS_DIR/${FILE_NAME}_filter.tmp.vcf
 


################## Section 2.5 VARIANT ANNOTATION & PRIORITISATION ################

echo "** Step 2.5 Annotation"

echo "A - Annovar"
## 1) A - Annotate variants using ANNOVAR 
## 2) A - Perform basic variant prioritization: filter to exonic variants not seen in dbSNP

## ***********  not finished

#cd
#wget *ANNOVAR LINK*
#unzip *ANNOVAR TAR*
ANN_DIR="/home/ubuntu/annovar"

#annotate_variation.pl -downdb -buildver hg19 -webfrom annovar refGene $ANN_DIR/humandb/
#annotate_variation.pl -downdb -buildver hg19 -webfrom annovar avsnp150 $ANN_DIR/humandb/

#vcf to annovar
#$ANN_DIR/convert2annovar.pl -format vcf4 $RESULTS_DIR/${FILE_NAME}_filtered.vcf.gz \
#      > $RESULTS_DIR/${FILE_NAME}_filtered.avinput





echo "B - snpEFF"
## 1) B - Annotate using snpEFF
## 2) B - Perform basic variant prioritization: filter to exonic variants not seen in dbSNP

## to install snpEff
#cd
#wget https://sourceforge.net/projects/snpeff/files/snpEff_v4_3t_core.zip
#unzip snpEff_v4_3t_core.zip

EFF_DIR="/home/ubuntu/snpEff/"

#download FULL dbSNP database
# these commands should only be performed once as the file is very large

#download dbSNP database
#cd $DATA_DIR/reference
#wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/00-common_all.vcf.gz
#gunzip -c 00-common_all.vcf.gz | bgzip -c > dbSNP-b151-small.vcf.bgz
#tabix -p vcf $DATA_DIR/reference/dbSNP-b151-small.vcf.bgz

#download the ref genome for hg19/GRCh37
#java -jar $EFF_DIR/snpeff.jar download hg19


##1) snpsift annotate to add dbsnp ids
echo "**** running snpsift annotate to add dbsnp ids"
java -jar $EFF_DIR/SnpSift.jar annotate -id -exists "DBSNP" $DATA_DIR/reference/dbSNP-b151-small.vcf.bgz \
        $RESULTS_DIR/${FILE_NAME}_filtered.vcf.gz  > $RESULTS_DIR/${FILE_NAME}_snpEff1_dbsnpids.vcf


##2) snpeff to annotate with functions against hg19/grch37 genome
#java parameter -xmx8g to increase the memory available to the java virtual machine to 4g.
echo "**** running snpeff"
java -Xmx8g -jar $EFF_DIR/snpEff.jar hg19 $RESULTS_DIR/${FILE_NAME}_snpEff1_dbsnpids.vcf \
	> $RESULTS_DIR/${FILE_NAME}_snpEff2_annotated.vcf


##3) snpsift filter to remove non-exonic and any that exist in dbsnp - retain selected and remove the rest
echo "**** running snpsift filter"
java -jar $EFF_DIR/SnpSift.jar filter \
        "ANN[*].EFFECT has 'missense_variant' & !exists DBSNP"\
         $RESULTS_DIR/${FILE_NAME}_snpEff2_annotated.vcf > $RESULTS_DIR/${FILE_NAME}_snpEff3_prioritised.vcf







### END



