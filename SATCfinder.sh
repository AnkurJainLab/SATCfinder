#!/bin/bash
set -e                        # Exit immediately on errors
#########################################################################################
##### Configure these variables for your datasets.
#########################################################################################
datasetPrefix=HEK-test        # Name this dataset
repeatType="CAG"              # Define your repeat type. Degenerate IUPAC bases are OK
minRepeatLength=3             # How many repeats to search & trim?
numCores=4                    # How many cores to run with?
removeTemporaryFiles=false     # Clean up temporary files created by this script? true/false

#Set up some file paths
#Your raw reads, such as the small example dataset
rawRead1=~/SATCfinder/WT_HEK_1.fq.gz
rawRead2=~/SATCfinder/WT_HEK_2.fq.gz
#The output directory. This will be made if it does not exist
outputDir=/lab/jain_imaging/Rachel/STARoutput/testSATCfinder

#Your STAR genome directory
genomeDir=/lab/jain_imaging/genomes/hg38_with_rRNA_overhang_100/
#If you need a spike-in fasta file, e.g. for mapping to a plasmid, add the path to the fasta file here.
fastaPath=""
#Directory for bbtools (https://sourceforge.net/projects/bbmap/)
bbdukDir=~/tools/bbmap_38.86/bbmap
#Directory for SATCfinder scripts. Default: the same directory as this script, but can be changed
SATCfinderDir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
#Directory for STAR executable, if not in $PATH
STARDir=""
#Directory for cutadapt, if not in $PATH
cutadaptDir=""


#########################################################################################
#### Begin SATCfinder pipeline.
#########################################################################################
#Set up paths for executables
outputDir=${outputDir%%+(/)}  # Remove trailing slash if present
mkdir -p ${outputDir}         # Make sure output dir exists

# Fix paths to have one trailing slash if needed
bbdukDir=${bbdukDir%%+(/)}
cutadaptDir=${cutadaptDir%%+(/)}
if [[ ${#cutadaptDir} -gt 0 ]]; then cutadaptPath=${cutadaptDir}/cutadapt; else cutadaptPath="cutadapt"; fi
STARDir=${STARDir%%+(/)}
if [[ ${#STARDir} -gt 0 ]]; then STARPath=${STARDir}/STAR; else STARPath="STAR"; fi
# Prepend -gFF if needed for fasta
if [[ ${#fasta} -gt 0 ]]; then fasta="--genomeFastaFiles ${fastaPath}"; fi

#Generate a repeat string for bbduk
repeatSeq=$(for i in $(seq 1 ${minRepeatLength}); do printf "${repeatType}"; done)

# Select reads with repeats using bbduk
${bbdukDir}/bbduk.sh \
    in=${rawRead1} \
    in2=${rawRead2} \
    outm=${outputDir}/${datasetPrefix}_bbduk_1.fastq.gz \
    outm2=${outputDir}/${datasetPrefix}_bbduk_2.fastq.gz \
    k=${#repeatSeq} \
    literal=${repeatSeq} \
    rcomp=t copyundefined hdist=0 mm=f \
    stats=${outputDir}/${datasetPrefix}_bbduk.log

# Remove sequencing adapters. This sequence provided is for Illumina-type adapters
${cutadaptPath} -a AGATCGGAAGAG \
    -o ${outputDir}/${datasetPrefix}_bbduk_cutadapt_1.fastq.gz \
    -p ${outputDir}/${datasetPrefix}_bbduk_cutadapt_2.fastq.gz \
    --cores=${numCores} \
    --error-rate=0.1 \
    --times=1 \
    --overlap=5 \
    --minimum-length=20 \
    --quality-cutoff=20 \
    ${outputDir}/${datasetPrefix}_bbduk_1.fastq.gz \
    ${outputDir}/${datasetPrefix}_bbduk_2.fastq.gz

# Move repeats to headers
${SATCfinderDir}/SATCfinder.py trim \
    --inFASTQ1=${outputDir}/${datasetPrefix}_bbduk_cutadapt_1.fastq.gz \
    --inFASTQ2=${outputDir}/${datasetPrefix}_bbduk_cutadapt_2.fastq.gz \
    --outSAM=${outputDir}/${datasetPrefix}.sam.gz \
    --minRepeats=${minRepeatLength} \
    --repeatSequence=${repeatType} \
    --outLog=${outputDir}/${datasetPrefix}_SATCfinder.log

# Align non-repetitive portion of reads to genome
${STARPath} --genomeDir ${genomeDir} \
    --runThreadN ${numCores} ${fasta}\
    --readFilesCommand zcat \
    --readFilesType SAM PE \
    --readFilesIn ${outputDir}/${datasetPrefix}.sam.gz \
    --outFileNamePrefix ${outputDir}/${datasetPrefix}_ \
    --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within

# Cleanup
# Make the output BAM file name nicer
#mv ${outputDir}/${datasetPrefix}_Aligned.sortedByCoord.out.bam ${outputDir}/${datasetPrefix}.bam
# Index BAM file
samtools index ${outputDir}/${datasetPrefix}.bam
# Select only reads that had repeats (aL/tL field is not 0), and stick them in a new indexed BAM file
# Here we only select primary alignments (-F 256)
samtools view -h -F 256 ${outputDir}/${datasetPrefix}.bam \
    | awk -F'\t' '{if((substr($1, 0, 1) == "@" || (($0 !~ /tL:i:0/) || ($0 !~ /aL:i:0/)))){print $0}}' \
    | samtools view -b >| ${outputDir}/${datasetPrefix}_selected.bam
samtools index ${outputDir}/${datasetPrefix}_selected.bam

# Remove temporary files
if [ $removeTemporaryFiles == true ]; then
    echo "Removing temporary files"
    rm -I ${outputDir}/${datasetPrefix}_bbduk_1.fastq.gz
    rm -I ${outputDir}/${datasetPrefix}_bbduk_2.fastq.gz
    rm -I ${outputDir}/${datasetPrefix}_bbduk_cutadapt_1.fastq.gz
    rm -I ${outputDir}/${datasetPrefix}_bbduk_cutadapt_2.fastq.gz
    rm -I ${outputDir}/${datasetPrefix}.sam.gz
    find ${outputDir}/${datasetPrefix}__STARtmp -type d -empty -delete
fi
