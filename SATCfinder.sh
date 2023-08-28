#!/bin/bash
#########################################################################################
################# These variables you should change for your datasets. ##################
#########################################################################################
#Name this dataset
datasetPrefix=HEK-test

## Define your repeat type and the minimum number of tandem repeats. IUPAC bases OK
repeatType=CAG
minRepeatLength=3

#Set up some file paths
#Your raw reads, such as the small example dataset
rawRead1=./WT_HEK_1.fq.gz
rawRead1=./WT_HEK_2.fq.gz
#The output directory. This will be made if necessary
outputDir=/lab/jain_imaging/Rachel/STARoutput/test
#Path to tools
bbdukDir=~/tools/bbmap_38.86/bbmap
SATCfinderDir=~/SATCfinder/
#Your STAR genome directory
genomeDir=/lab/jain_imaging/genomes/hg38_with_rRNA_overhang_100/
#If you need a spike-in fasta file, like for mapping to a plasmid, add the path to the fasta here.
fasta=""

#How many cores to run with
numCores=4

#########################################################################################
##################### Here's where the SATCfinder pipeline begins #######################
#########################################################################################
#Make sure output dir exists
mkdir -p ${outputDir}

#Generate a repeat string for bbduk
repeatSeq=`for i in $(seq 1 ${minRepeatLength}); do printf "${repeatType}"; done`

#Select reads with repeats using bbduk
${bbdukDir}/bbduk.sh \
			in=${rawRead1} \
			in2=${rawRead2} \
			outm=${outputDir}/${datasetPrefix}_bbduk_1.fastq.gz \
			outm2=${outputDir}/${datasetPrefix}_bbduk_2.fastq.gz \
			k=${#repeatSeq} \
			literal=${repeatSeq} \
			rcomp=t copyundefined hdist=0 mm=f \
			stats=${outputDir}/${datasetPrefix}_bbduk.log

#remove sequencing adapters. This should work for illumina-type adapters
cutadapt -a AGATCGGAAGAG \
      -o ${outputDir}/${datasetPrefix}_cutadapt_1.fastq.gz \
      -p ${outputDir}/${datasetPrefix}_cutadapt_2.fastq.gz \
			--cores=${numCores} \
			--error-rate=0.1 \
			--times=1 \
			--overlap=5 \
			--minimum-length=20 \
			--quality-cutoff=20 \
			${outputDir}/${datasetPrefix}_bbduk_1.fastq.gz \
			${outputDir}/${datasetPrefix}_bbduk_2.fastq.gz

${SATCfinderDir}/trimRepeatsAndConvertFASTQToSAM.py --sra=${study} --disease=${ID} \
			--in1=${outputDir}/${datasetPrefix}_bbduk_cutadapt_1.fastq.gz \
			--in2=${outputDir}/${datasetPrefix}_bbduk_cutadapt_2.fastq.gz \
			--out=${outputDir}/${datasetPrefix}.sam.gz \
			--minRepeats=${minRepeatLength} --repeatSequence=${repeatType} \
			--log=${outputDir}/${datasetPrefix}_repeatselected.log

STAR --genomeDir ${genomeDir} \
			--runThreadN ${numCores} ${fasta}\
			--readFilesCommand zcat \
			--readFilesType SAM PE \
			--readFilesIn ${outputDir}/${datasetPrefix}.sam.gz \
			--outFileNamePrefix ${outputDir}/${datasetPrefix}_ \
			--outSAMtype BAM SortedByCoordinate --outSAMunmapped Within

#Cleanup
mv ${outputDir}/${datasetPrefix}_Aligned.sortedByCoord.out.bam ${outputDir}/${datasetPrefix}.bam
samtools index ${outputDir}/${datasetPrefix}.bam
#Select only reads that had repeats (ignore the mates)
samtools view -h -F 256 ${outputDir}/${datasetPrefix}.bam \
      | awk -F'\t' '{if((substr(\$1, 0, 1) == \"@\" || ((\$0 !~ /tL\:i\:0/) || (\$0 !~ /aL\:i\:0/)))){print \$0}}' \
      | samtools view -b >| ${outputDir}/${datasetPrefix}_selected_3x.bam
samtools index ${outputDir}/${datasetPrefix}_selected_3x.bam
samtools index ${outputDir}/${datasetPrefix}_selected_4x.bam
rm -I ${outputDir}/${datasetPrefix}_bbduk_1.fastq.gz
rm -I ${outputDir}/${datasetPrefix}_bbduk_2.fastq.gz
rm -I ${outputDir}/${datasetPrefix}_bbduk_cutadapt_1.fastq.gz
rm -I ${outputDir}/${datasetPrefix}_bbduk_cutadapt_2.fastq.gz
rm -I ${outputDir}/${datasetPrefix}.sam.gz
find ${outputDir}/${datasetPrefix}__STARtmp -type d -empty -delete