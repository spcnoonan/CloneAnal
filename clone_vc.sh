#!/bin/bash

#SBATCH --job-name=CloneAnalysis
#SBATCH --output=/storage/users/cnoonan22/experiments/EXP0001/CloneAnal/clone_analysis.log
#SBATCH --mem=64000
#SBATCH --cpus-per-task=64

###############################################################################
#             
#
#                 *** Pipeline for genotyping clones ***
#
#   This script was designed to analyse sequencing data obtained
#   from cross-bred Eimeria tenella clones. These data were produced by 
#   Justin, who cross-bred the Wisconsin and Weybridge strains of 
#   E. tenella and sequenced the offspring using Illumina short
#   read sequencing. 
#
#   The workflow is as follows:
#       1.  Reads are quality trimmed using BBTools.
#       2.  Reads are aligned to the Houghton E. tenella reference genome.
#       3.  Resulting BAM file is sorted, indexed, duplicates are marked
#           and read groups are added.
#       4.  Resulting BAM file is used as input into a GATK workflow for
#           short variant discovery
#
#   NOTE 1: Many of the processes use Java, and memory allocation for these
#   processes are assuming that they are free to use a maximum of 64G
#   memory. Change if less/more memory will be available to you.
#
#   NOTE 2: The value for MFH can be assessed using the ulimit command.
#
###############################################################################

set -e

#.......................................................
#					SET PATHS
#  These are paths that you will be required to change to suit your particular paths
#.......................................................
# Path to folder with raw reads
READS=/storage/users/cnoonan22/experiments/EXP0001/Data
# Path to working directory
OUTPATH=/storage/users/cnoonan22/experiments/EXP0001/CloneAnal
# Path to INDEXED E. tenella Houghton reference
REF=/storage/users/cnoonan22/experiments/EXP0001/Resources/References/GCA_905310635.1_pEimTen1.1_genomic.fna
# Path to scripts folder
SCRIPTS=/storage/users/cnoonan22/experiments/EXP0001/Scripts
# Path to folder containing bbduk tool
BBDUK=/storage/users/cnoonan22/experiments/EXP0001/Bin/bbmap/bbduk.sh
# Path to adaptor file
ADPTR=/storage/users/cnoonan22/experiments/EXP0001/Bin/bbmap/resources/adapters.fa
# Path to Picard JAR file
PICARD=/storage/users/cnoonan22/experiments/EXP0001/Bin/picard-3.2.0/picard.jar
# Value for MAX_FILE_HANDLES parameter in MarkDuplicates
MFH="1000"
# Path to GATK JAR file
GATK=/storage/users/cnoonan22/experiments/EXP0001/Bin/gatk-4.6.0.0/gatk
# Number of threads - should be equal to number of cpus requested
THREADS="64"
# Path to sample map for GenomicsDBImport
SMAP=/storage/users/cnoonan22/experiments/EXP0001/Resources/GenomeDB/clone.sample_map
# Path to interval list for GenomicsDBImport
INTLIST=/storage/users/cnoonan22/experiments/EXP0001/Resources/GenomeDB/ETenRef.intervals
# Path to SnpEFF JAR file
SNPEFF=/storage/users/cnoonan22/experiments/EXP0001/Bin/snpEff/snpEff.jar
# Path to SnpSift JAR file
SNPSIFT=/storage/users/cnoonan22/experiments/EXP0001/Bin/snpEff/SnpSift.jar
# SnpEFF database name
SNPEFFDATABASE=pEimTen1.1


#........................................................
#			CHECKS
#........................................................

echo "activating assembly conda environment..."
module load apps/anaconda-4.7.12.tcl
unset PYTHONPATH
source /storage/apps/anaconda-201910/etc/profile.d/conda.sh
conda activate /storage/users/cnoonan22/bin/mambaforge/envs/clones
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

# Making folder to store the output of per-sample processing
GVCF=$PERSAMPLE/GVCF
mkdir -p $GVCF

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
        echo "A: $ID"
            
        # Making directory to store output files
        CLONE_DIR=$PERSAMPLE/${ID}
        mkdir -p $CLONE_DIR
        cd $CLONE_DIR

        # Making directory to store logs and stats for each clone
        CLONE_STATS=$CLONE_DIR/${ID}_Stats
        CLONE_LOGS=$CLONE_DIR/${ID}_Logs
        mkdir -p $CLONE_STATS
        mkdir -p $CLONE_LOGS

        # Trimming files            
        TRIMMING=$CLONE_DIR/Trimming
        mkdir $TRIMMING
        cd $TRIMMING
        echo "Running trimming step..."
        trim="$BBDUK in1=$file_R1 in2=$file_R2 ref=$ADPTR k=23 ftm=5 hdist=1 out1=${ID}_clean1.fastq out2=${ID}_clean2.fastq qtrim=rl trimq=10 minlen=20 2>&1 | tee $CLONE_LOGS/bbduk.txt"
        echo $trim
        eval $trim
        echo "Trimming complete"

        # Read alignment
        ALIGNMENT=$CLONE_DIR/Alignments
        mkdir -p $ALIGNMENT
        cd $ALIGNMENT
        echo "Running alignment step..."
        algn="bwa mem -t $THREADS $REF $TRIMMING/${ID}_clean1.fastq $TRIMMING/${ID}_clean2.fastq | samtools view -@ $THREADS -b - | samtools sort -@ $THREADS -o ${ID}VsRef.bam"
        echo $algn
        eval $algn
        samtools index ${ID}VsRef.bam
        echo "unmapped extraction done"

        # Check mapping quality
        echo "Checking the mapping quality with Qualimap ..."
        qual="qualimap bamqc -nt $THREADS -bam ${ID}VsRef.bam -outdir $CLONE_STATS/Qualimap"
        echo $qual
        eval $qual

        # Picard workflow
        readg="java -Xmx64g -jar $PICARD AddOrReplaceReadGroups I=${ID}VsRef.bam O=${ID}VsRef_RG.bam RGLB=library RGPL=illumina RGPU=barcode RGSM=sample 2>&1 | tee $CLONE_LOGS/AddOrReplaceReadGroups.txt"
        echo $readg
        eval $readg
        markdup="java -Xmx64g -jar $PICARD MarkDuplicates I=${ID}VsRef_RG.bam O=${ID}VsRef_MD.bam METRICS_FILE=${ID}VsRef_MD.metrics MAX_FILE_HANDLES=$MFH 2>&1 | tee $CLONE_LOGS/MarkDuplicates.txt"
        echo $markdup
        eval $markdup
        samtools index ${ID}VsRef_MD.bam

        # GATK workflow
        VARIANT=$CLONE_DIR/VariantCalling
        mkdir -p $VARIANT
        cd $VARIANT
        #varcall="$GATK --java-options "-Xmx64g" HaplotypeCaller -R $REF -I $ALIGNMENT/${ID}VsRef_MD.bam -O ${ID}VsRef.g.vcf -ploidy 1 -ERC GVCF -bamout ${ID}_HC_Haplotypes.bam 2>&1 | tee $CLONE_LOGS/HaplotypeCaller.txt"
        varcall="$GATK --java-options "-Xmx64g" HaplotypeCaller -R $REF -I $ALIGNMENT/${ID}VsRef_MD.bam -O ${ID}VsRef.g.vcf --min-pruning 1 --min-dangling-branch-length 1 -ploidy 1 -ERC GVCF -bamout ${ID}_HC_Haplotypes.bam 2>&1 | tee $CLONE_LOGS/HaplotypeCaller.txt"
        echo $varcall
        eval $varcall

        # Copying GVCF file to a pooled GVCF folder for multi-sample analysis
        cp ${ID}VsRef.g.vcf $GVCF/${ID}.g.vcf
        cp ${ID}VsRef.g.vcf.idx $GVCF/${ID}.g.vcf.idx

    else
        echo "Pair file not found for $file"
    fi
done

#........................................................
#			GENOTYPING 
#........................................................

echo "Per-sample processes complete. Results will be stored in a folder for each clone, e.g. 'Clone_1'."
echo
echo "Moving to directory ${COMBO}"
cd $COMBO
echo
echo "Generating a GenomicsDB database ..."

# Generate a directory containing a GenomicsDB database with combined multi-sample data
# NOTE - Make sure Sample Map contains correct path and files are .idx indexed
gendb="$GATK --java-options "-Xmx64g" GenomicsDBImport -R $REF --genomicsdb-workspace-path CloneDB --sample-name-map $SMAP -L $INTLIST --reader-threads $THREADS 2>&1 | tee $LOGS/GenomicsDBImport.txt"
echo $gendb
eval $gendb

# Joint genotyping on cohort of clones
geno="$GATK --java-options "-Xmx64g" GenotypeGVCFs -R $REF -V gendb://CloneDB -O AllClones.vcf"
echo $geno
eval $geno

#........................................................
#			VARIANT FILTERING 
#........................................................

UNFILT=$COMBO/UnfilteredVars
FILT=$COMBO/FilteredVars
mkdir -p $UNFILT
mkdir -p $FILT

# SNP quality tagging and filtering

sep_SNP="$GATK --java-options "-Xmx64g" SelectVariants -V AllClones.vcf -select-type SNP -O CombinedSNPs.vcf"
echo "Generating SNP file ..."
echo
echo $sep_SNP
eval $sep_SNP
echo "SNPs are stored in file 'GenotypeGVCF_Output_SNP.vcf'"
echo "Filtering SNPs and flagging low-confidence calls ..."
echo 

SNP_filt="$GATK VariantFiltration --java-options "-Xmx64g" -V CombinedSNPs.vcf -R $REF --filter-expression 'QD < 2.0' --filter-name 'QD2' --filter-expression 'SOR > 3.0' --filter-name 'SOR3' --filter-expression 'FS > 60.0' --filter-name 'FS60' --filter-expression 'MQ < 40.0' --filter-name 'MQ40' --output TaggedCombinedSNPs.vcf"
echo $SNP_filt
eval $SNP_filt
echo
SNP_select="$GATK SelectVariants -V TaggedCombinedSNPs.vcf --exclude-filtered -O FilteredCombinedSNPs.vcf"
echo $SNP_select
eval $SNP_select
echo

mv TaggedCombinedSNPs.vcf* CombinedSNPs.vcf* $UNFILT

echo "SNP filtering complete - results stored in 'FilteredCombinedSNPs.vcf'"
echo

# INDEL quality tagging and filtering

sep_INDEL="$GATK --java-options "-Xmx64g" SelectVariants -V AllClones.vcf -select-type INDEL -O CombinedINDELs.vcf"
echo "Generating INDEL file ..."
echo
echo $sep_INDEL
eval $sep_INDEL
echo
echo "INDELs are stored in file 'CombinedINDELs.vcf'"
echo "Filtering INDELs and flagging low-confidence calls ..."
echo 
INDEL_filt="$GATK VariantFiltration --java-options "-Xmx64g" -V CombinedINDELs.vcf -R $REF --filter-expression 'QD < 2.0' --filter-name 'QD2' --filter-expression 'FS > 200.0' --filter-name 'FS200' --filter-expression 'MQ < 40.0' --filter-name 'MQ40' --output TaggedCombinedINDELs.vcf"
echo $INDEL_filt
eval $INDEL_filt
echo
INDEL_select="$GATK SelectVariants -V TaggedCombinedINDELs.vcf --exclude-filtered -O FilteredCombinedINDELs.vcf"
echo $INDEL_select
eval $INDEL_select
echo

mv TaggedCombinedINDELs.vcf* CombinedINDELs.vcf* $UNFILT

echo "INDEL filtering complete - results stored in 'FilteredCombinedINDELs.vcf'"
echo

# Moving remaining files to their respective directories

mv FilteredCombined* $FILT
mv AllClones.vcf* $UNFILT

# Merging SNP and INDEL files into a single VCF and sorting filtered variants 

cd $FILT
var_merge="java -Xmx64g -jar $PICARD MergeVcfs I=FilteredCombinedSNPs.vcf I=FilteredCombinedINDELs.vcf O=SNP_INDEL_Merged.vcf"
echo $var_merge
eval $var_merge
echo
var_sort="java -jar $PICARD SortVcf I=SNP_INDEL_Merged.vcf O=SNP_INDEL_Sorted.vcf"
echo $var_sort
eval $var_sort
echo

# Removing variants for which there are more than 1 alternate alleles and are present in either both parental strains or neither

var_filt="python3 $SCRIPTS/var_filter.py SNP_INDEL_Sorted.vcf SNP_INDEL_Processed.vcf"
echo $var_filt
eval $var_filt
echo
var_index="$GATK IndexFeatureFile -I SNP_INDEL_Processed.vcf"
echo $var_index
eval $var_index
echo

#........................................................
#			VARIANT ANNOTATION 
#........................................................

ANNOT=$COMBO/AnnotatedVars
mkdir -p $ANNOT

cd $ANNOT
echo "Annotating SNP_INDEL_Processed.vcf ..."

java -Xmx64g -jar $SNPEFF -v -o vcf -no-downstream -no-intergenic -no-intron -no-upstream -no-utr -stats $STATS/SNP_stats.html $SNPEFFDATABASE $FILT/SNP_INDEL_Processed.vcf > SNP_INDEL_Processed_Coding.ann.vcf
java -Xmx64g -jar $SNPEFF -v -o vcf -no-downstream -no-intergenic -no-upstream -stats $STATS/SNP_stats.html $SNPEFFDATABASE $FILT/SNP_INDEL_Processed.vcf > SNP_INDEL_Processed_Genic.ann.vcf
java -Xmx64g -jar $SNPEFF -v -o vcf -stats $STATS/SNP_stats.html $SNPEFFDATABASE $FILT/SNP_INDEL_Processed.vcf > SNP_INDEL_Processed_All.ann.vcf

$GATK IndexFeatureFile -I SNP_INDEL_Processed_Coding.ann.vcf
$GATK IndexFeatureFile -I SNP_INDEL_Processed_Genic.ann.vcf
$GATK IndexFeatureFile -I SNP_INDEL_Processed_All.ann.vcf

echo
echo "Annotation complete"
echo

WEYVARS=$ANNOT/WeybridgeVars
WISVARS=$ANNOT/WisconsinVars
mkdir -p $WEYVARS
mkdir -p $WISVARS

echo "Splitting variants into Weybridge-derived and Wisconsin-derived ..."
$SCRIPTS/var_split.py SNP_INDEL_Processed_Coding.ann.vcf
mv WeyVars.vcf $WEYVARS/WeyVars_Coding.vcf
mv WisVars.vcf $WISVARS/WisVars_Coding.vcf
$GATK IndexFeatureFile -I $WEYVARS/WeyVars_Coding.vcf
$GATK IndexFeatureFile -I $WISVARS/WisVars_Coding.vcf

$SCRIPTS/var_split.py SNP_INDEL_Processed_Genic.ann.vcf
mv WeyVars.vcf $WEYVARS/WeyVars_Genic.vcf
mv WisVars.vcf $WISVARS/WisVars_Genic.vcf
$GATK IndexFeatureFile -I $WEYVARS/WeyVars_Genic.vcf
$GATK IndexFeatureFile -I $WISVARS/WisVars_Genic.vcf

$SCRIPTS/var_split.py SNP_INDEL_Processed_All.ann.vcf
mv WeyVars.vcf $WEYVARS/WeyVars_All.vcf
mv WisVars.vcf $WISVARS/WisVars_All.vcf
$GATK IndexFeatureFile -I $WEYVARS/WeyVars_All.vcf
$GATK IndexFeatureFile -I $WISVARS/WisVars_All.vcf

# Extracting only the protein-coding variants
cat $WEYVARS/WeyVars_Coding.vcf | java -jar $SNPSIFT filter "ANN[0].FEATURE has 'transcript'" > $WEYVARS/WeyPCVars_only.vcf
cat $WISVARS/WisVars_Coding.vcf | java -jar $SNPSIFT filter "ANN[0].FEATURE has 'transcript'" > $WISVARS/WisPCVars_only.vcf

