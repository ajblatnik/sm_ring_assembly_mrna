#!/bin/bash

# this script was used to process Sm-RIP-Seq libraries from Takara SMARTer Stranded Total RNA-Seq Kit v3 (product # 634451)

# Set paths to the indexed genome and annotation file (GTF)
index='/filepath to indexed genome'
gtf='/filepath to GTF/filename.gtf'

# Define output directories for the pipeline
outdir='/filepath to output directory''
rawout=$outdir/rawout              # Directory for raw fastq files
fastqs=$outdir/fastqs              # Directory for trimmed fastq files
star=$outdir/star                  # Directory for STAR alignment output
counts=$outdir/feature_counts      # Directory for featureCounts output
logs=$outdir/logs                  # Directory for log files
runscript=$outdir/sbatch           # Directory for SLURM sbatch scripts

# Create the output directories if they don't already exist
[ ! -d $outdir ] && mkdir -p $outdir
[ ! -d $rawout ] && mkdir -p $rawout
[ ! -d $fastqs ] && mkdir -p $fastqs
[ ! -d $star ] && mkdir -p $star
[ ! -d $counts ] && mkdir -p $counts
[ ! -d $logs ] && mkdir -p $logs
[ ! -d $runscript ] && mkdir -p $runscript

# Define system resource parameters
cpu=18           # Number of CPU cores to use
mem="128gb"      # Amount of memory to allocate

# Loop through all fastq files in the rawout directory that match *_1.fq
for fq in $(ls $rawout/*_1.fq)
do 
    # Extract the base name (without _1.fq) for the current sample
    name=$(basename $fq _1.fq)
    
    # Define paths to the paired fastq files (R1 and R2)
    R1=$fastqs/$name\_1.fq
    R2=$fastqs/$name\_2.fq
    
    # Define paths to the output BAM and feature count files
    outbam=$star/$name\_Aligned.sortedByCoord.out.bam
    outcount=$counts/$name.featureCounts.tsv

    # Print some information about the current sample being processed
    echo "========================================="
    echo $name
    echo $R1
    echo $R2
    echo $outbam 
    echo $outcount
    echo "========================================="
    echo
    
    # Generate the SLURM sbatch script for this sample
    echo """#!/bin/bash

#SBATCH -t 04:00:00              # Set time limit for the job
#SBATCH -c $cpu                  # Request CPU cores
#SBATCH --mem $mem               # Request memory
#SBATCH -J $name                 # Job name
#SBATCH -o $logs/$name.log       # Output log file
#SBATCH --account=PAS0631        # Account for job submission

# Load necessary modules (e.g., miniconda3, fastqc, star)
module load star

# Change to the rawout directory and activate the conda environment
cd rawout
source activate 'conda environment name'

# Trim the paired fastq files using cutadapt
cutadapt \\
    -m 20 \\
    --nextseq-trim=20 \\
    -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG \\
    -A NNNNNNNNNNNNNNAATGATACGGCGACCACCGAGATCTACACNNNNNNNNACACTCTTTCCCTACACGACGCTCTTCCGATCT \\
    -o $fastqs $name\_trimmed_1.fq -p $name\_trimmed_2.fq $rawout $name\_1.fq $name\_2.fq

source deactivate
cd ..

# Perform alignment using STAR
STAR \\
    --runThreadN $cpu \\
    --genomeDir $index \\
    --outFileNamePrefix $star/$name\_ \\
    --outFilterType BySJout \\
    --outFilterMismatchNoverLmax 0.6 \\
    --alignIntronMin 20 \\
    --alignIntronMax 10000 \\
    --outFilterMultimapNmax 999 \\
    --outSAMtype BAM SortedByCoordinate \\
    --readFilesIn $R1 $R2 \\
    --readFileCommand zcat

# Index the aligned BAM file using samtools
samtools index -@ $cpu $outbam

# Perform feature counting using featureCounts
source activate 'conda environment name'
featureCounts \\
    -T $cpu \\
    -t exon \\
    -g gene_id \\
    -a $gtf \\
    -s 1 \\
    -p \\
    -M \\
    --fraction \\
    -o $outcount.tmp \\
    $outbam 2> $outcount.log

# Process the output from featureCounts to generate the final count file
cat $outcount.tmp | grep -v \"^#\" | cut -f 1,7 | sed -e 's/_Aligned_SortedByCoord.bam//g' | awk '(NR>1)' > $outcount

# Remove the temporary count file
rm $outcount.tmp

# Compile QC results using MultiQC
source activate 'conda environment name'
multiqc fastqout
multiqc fastqs
multiqc star
multiqc counts
source deactivate

""" >$runscript/$name.sbatch
    
    # Submit the sbatch script to SLURM
    sbatch $runscript/$name.sbatch

done
