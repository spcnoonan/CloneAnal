#!/bin/bash

#SBATCH --job-name=Freebayes
#SBATCH --output=/storage/users/cnoonan22/experiments/EXP0001/scratch/Freebayes/freebayes_anal.log
#SBATCH --mem=150000
#SBATCH --cpus-per-task=64

###############################################################################
#             
#
#                 *** Pipeline for genotyping clones ***
#
#   This script was designed to analyse sequencing data obtained
#   from E. tenella recombinant clones selected for arprinocid resistance and precocity,
#   isolated during the make of the E. tenella linkage map.
#
#   The workflow is as follows:
#       1.  Reads are quality trimmed using BBTools.
#       2.  Reads are aligned to the Houghton E. tenella reference genome.
#       3.  Resulting BAM file is sorted, indexed, duplicates are marked
#           and read groups are added.
#       4.  Resulting BAM file is used as input into a freebayes workflow for
#           short variant discovery.
#       5.  Variants are filtered to remove multi-allelic variants, variants with a 
#           missing genotype rate of ~15% or more and variants which are present in both strains.
#       6.  Variants are annotated using SnpEff.
#
#   NOTE 1: Many of the processes use Java, and memory allocation for these
#   processes are assuming that they are free to use a maximum of 64G
#   memory. Change if less/more memory will be available to you.
#
#   NOTE 2: The value for MFH can be assessed using the ulimit command.
#
#   NOTE 3: Stringent read alignment parameters are used to minimise subsequent false positive variant calls.
#           Additional stringent parameters are used in the freebayes workflow and subsequent filtering steps 
#           to minimise false positive variant calls.
#
###############################################################################

set -e

#.......................................................
#					SET PATHS
#  These are paths that you will be required to change to suit your particular paths
#.......................................................
# Path to folder with raw reads
READS=/storage/users/cnoonan22/experiments/EXP0001/data
# Path to working directory
OUTPATH=/storage/users/cnoonan22/experiments/EXP0001/scratch/Freebayes
# Path to INDEXED E. tenella Houghton reference
REF=/storage/users/cnoonan22/experiments/EXP0001/scratch/WETis_VC/Ref/GCA_905310635.1_pEimTen1.1_genomic.fna
# Path to scripts folder
#SCRIPTS=
# Path to folder containing bbduk tool
BBDUK=/storage/users/cnoonan22/experiments/EXP0001/bin/bbmap/bbduk.sh
# Path to adaptor file
ADPTR=/storage/users/cnoonan22/experiments/EXP0001/bin/bbmap/resources/adapters.fa
# Path to Picard JAR file
PICARD=/storage/users/cnoonan22/experiments/EXP0014/bin/picard.jar
# Value for MAX_FILE_HANDLES parameter in MarkDuplicates
MFH="1000"
# Number of threads - should be equal to number of cpus requested
THREADS="64"


#........................................................
#			CHECKS
#........................................................

echo "activating sample_vc micromamba environment..."
unset PYTHONPATH
source /storage/users/cnoonan22/micromamba/etc/profile.d/mamba.sh
micromamba activate /storage/users/cnoonan22/micromamba/envs/sample_vc
echo

#........................................................
#			BEGINNING 
#........................................................

cd $OUTPATH
echo "Your working directory is ${OUTPATH}"

# Making folder for storing per-sample and combined data
PERSAMPLE=$OUTPATH/PerSample
COMBO=$OUTPATH/Combined
mkdir -p $PERSAMPLE
mkdir -p $COMBO
cd $PERSAMPLE

# Making folders to store logs and stats
STATS=$OUTPATH/Stats
LOGS=$OUTPATH/Logs
mkdir -p $STATS
mkdir -p $LOGS

#........................................................
#			PER-CLONE PROCESSING 
#........................................................

# Iterate over the files in the directory
for file_R1 in "$READS"/*_R1_001.fastq.gz; do
    # Extract the base name of the file (without the directory path)
    base_R1=$(basename "$file_R1")
    
    # Extract fields A, B, and C using parameter substitution and pattern matching
    ID=$(echo "$base_R1" | cut -d'_' -f1)
    PAIRNUM=$(echo "$base_R1" | cut -d'_' -f2)
    
    # Construct the corresponding R2 filename by replacing R1 with R2
    base_R2="${ID}_${PAIRNUM}_L001_R2_001.fastq.gz"
    file_R2="$READS/$base_R2"
    
    # Check if the R2 file exists
    if [[ -f "$file_R2" ]]; then

        # Print the matched pair and the value of A
        echo "Matched pair found:"
        echo "R1 file: $file_R1"
        echo "R2 file: $file_R2"
        echo "ID: $ID"

        # Making directory to store output files
        SAMPLE_DIR=$PERSAMPLE/${ID}
        mkdir -p $SAMPLE_DIR
        cd $SAMPLE_DIR

        # Making directory to store logs and stats for each clone
        SAMPLE_STATS=$SAMPLE_DIR/${ID}_Stats
        SAMPLE_LOGS=$SAMPLE_DIR/${ID}_Logs
        mkdir -p $SAMPLE_STATS
        mkdir -p $SAMPLE_LOGS

        # Trimming files            
        TRIMMING=$SAMPLE_DIR/Trimming
        mkdir -p $TRIMMING
        cd $TRIMMING
        echo "Running trimming step..."
        trim="$BBDUK in1=$file_R1 in2=$file_R2 ref=$ADPTR k=23 ftm=5 hdist=1 out1=${ID}_clean1.fastq out2=${ID}_clean2.fastq qtrim=rl trimq=10 minlen=20 2>&1 | tee $SAMPLE_LOGS/bbduk.txt"
        echo $trim
        eval $trim
        echo "Trimming complete"

        # Read alignment
        ALIGNMENT=$SAMPLE_DIR/Alignments
        mkdir -p $ALIGNMENT
        cd $ALIGNMENT
        echo "Running alignment step..."
        algn="bwa mem -t $THREADS -T 135 $REF $TRIMMING/${ID}_clean1.fastq $TRIMMING/${ID}_clean2.fastq | samtools view -@ $THREADS -b - | samtools sort -@ $THREADS -o ${ID}VsRef.bam"
        echo $algn
        eval $algn
        samtools index ${ID}VsRef.bam
        echo "unmapped extraction done"

        # Check mapping quality
        echo "Checking the mapping quality with Qualimap ..."
        qual="qualimap bamqc -nt $THREADS -bam ${ID}VsRef.bam -outdir $SAMPLE_STATS/Qualimap"
        echo $qual
        eval $qual

        # Picard workflow
        readg="java -Xmx64g -jar $PICARD AddOrReplaceReadGroups I=${ID}VsRef.bam O=${ID}VsRef_RG.bam RGLB=LIB1 RGPL=ILLUMINA RGPU=barcode RGSM=${ID} RGID=${ID} 2>&1 | tee $SAMPLE_LOGS/AddOrReplaceReadGroups.txt"
        echo $readg
        eval $readg
        markdup="java -Xmx64g -jar $PICARD MarkDuplicates I=${ID}VsRef_RG.bam O=${ID}VsRef_MD.bam METRICS_FILE=${ID}VsRef_MD.metrics MAX_FILE_HANDLES=$MFH 2>&1 | tee $SAMPLE_LOGS/MarkDuplicates.txt"
        echo $markdup
        eval $markdup
        samtools index ${ID}VsRef_MD.bam
        realpath ${ID}VsRef_MD.bam >> $COMBO/bamlist.txt

    else
        echo "Pair file not found for $file"
    fi
done

echo "Per-sample processes complete. Results will be stored in a folder for each clone, e.g. 'Clone_1'."
echo
echo "Moving to directory ${COMBO}"
cd $COMBO
echo

# Freebayes workflow
VARIANT=$COMBO/VariantCalling
mkdir -p $VARIANT
cd $VARIANT
varcall="freebayes -f $REF -p 1 --min-mapping-quality 20 --min-coverage 3 --min-base-quality 20 -L $COMBO/bamlist.txt -v combined.vcf 2>&1 | tee $LOGS/freebayes.txt"
echo $varcall
eval $varcall

# Variant Filtering
UNFILT=$COMBO/UnfilteredVars
FILT=$COMBO/FilteredVars
mkdir -p $UNFILT
mkdir -p $FILT

cd $COMBO/FilteredVars
var_filt="bcftools view --max-alleles 2 ${VARIANT}/combined.vcf > combined_nomulti.vcf"
echo $var_filt
eval $var_filt
echo

# Removing variants which have a missing genotype rate of ~15% (i.e. 4 samples or less missing)
var_filt='bcftools filter -e "AN < 17" combined_nomulti.vcf | bcftools filter -i "AF < 0" > combined_nomiss.vcf'
echo $var_filt
eval $var_filt
echo

# Filtering out variants which are present in both strains and splitting the output into two separate VCF files for each strain
var_py_filt="python3 ${SCRIPTS}/var_filter_freebayes.py combined_nomiss.vcf combined_pyfilt.vcf"
echo $var_py_filt
eval $var_py_filt
echo

#........................................................
#			VARIANT ANNOTATION 
#........................................................

ANNOT=$COMBO/AnnotatedVars
mkdir -p $ANNOT

cd $ANNOT
echo "Annotating variants in combined_pyfilt.vcf ..."

java -Xmx64g -jar $SNPEFF -v -o vcf -no-downstream -no-intergenic -no-intron -no-upstream -no-utr -stats $STATS/SNP_stats_coding.html $SNPEFFDATABASE $COMBO/FilteredVars/combined_pyfilt.vcf > combined_coding.ann.vcf
java -Xmx64g -jar $SNPEFF -v -o vcf -no-downstream -no-intergenic -no-upstream -stats $STATS/SNP_stats_genic.html $SNPEFFDATABASE $COMBO/FilteredVars/combined_pyfilt.vcf > combined_genic.ann.vcf
java -Xmx64g -jar $SNPEFF -v -o vcf -stats $STATS/SNP_stats_all.html $SNPEFFDATABASE $COMBO/FilteredVars/combined_pyfilt.vcf > combined_all.ann.vcf

$GATK IndexFeatureFile -I combined_coding.ann.vcf
$GATK IndexFeatureFile -I combined_genic.ann.vcf
$GATK IndexFeatureFile -I combined_all.ann.vcf

echo
echo "Annotation complete"
echo

#........................................................
#			SPLITTING VARIANTS
#........................................................

# Splitting files
SPLIT=$COMBO/SplitVariants
mkdir -p $SPLIT
cd $SPLIT

WEYVARS=$SPLIT/WeybridgeVars
WISVARS=$SPLIT/WisconsinVars
mkdir -p $WEYVARS $WISVARS

echo
echo "Splitting variants into Weybridge-derived and Wisconsin-derived ..."
echo

var_split1="python3 ${SCRIPTS}/var_split_freebayes.py ${ANNOT}/combined_coding.ann.vcf"
echo $var_split1
eval $var_split1
mv WisVars.vcf $SPLIT/WisconsinVars/WisCodingAnnot.vcf
mv WeyVars.vcf $SPLIT/WeybridgeVars/WeyCodingAnnot.vcf
echo

var_split2="python3 ${SCRIPTS}/var_split_freebayes.py ${ANNOT}/combined_genic.ann.vcf"
echo $var_split2
eval $var_split2
mv WisVars.vcf $SPLIT/WisconsinVars/WisGenicAnnot.vcf
mv WeyVars.vcf $SPLIT/WeybridgeVars/WeyGenicAnnot.vcf
echo

var_split3="python3 ${SCRIPTS}/var_split_freebayes.py ${ANNOT}/combined_all.ann.vcf"
echo $var_split3
eval $var_split3
mv WisVars.vcf $SPLIT/WisconsinVars/WisAllAnnot.vcf
mv WeyVars.vcf $SPLIT/WeybridgeVars/WeyAllAnnot.vcf
echo

echo "Pipeline Complete."