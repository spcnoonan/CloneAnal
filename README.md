# Genotyping *Eimeria tenella* Recombinant Clones

A variant-calling pipeline for genotyping *E. tenella* recombinant clones isolated during the construction of the *E. tenella* linkage map. The clones come from a cross between the **Weybridge** and **Wisconsin** strains and were selected for arprinocid resistance and precocity.

The pipeline takes paired-end Illumina reads for each clone and the two parental strains through trimming, alignment, joint haploid variant calling, filtering and annotation. It then assigns each variant to the parent it came from. A companion R script uses PLINK outputs from those variants to look for genomic regions where linkage disequilibrium (LD) and allele-frequency skew are both high, which are candidate loci under selection.

Strict alignment, calling and filtering parameters are used throughout to keep false-positive variant calls to a minimum.

---

## Repository contents

| File | Role |
|---|---|
| [`clone_freebayes.sh`](clone_freebayes.sh) | Main SLURM pipeline: reads → filtered, annotated, parent-split VCFs |
| [`var_filter_freebayes.py`](var_filter_freebayes.py) | Companion script: keeps only variants that differ between the two parents |
| [`var_split_freebayes.py`](var_split_freebayes.py) | Companion script: splits variants by parent of origin |
| [`linkage_af_analysis.R`](linkage_af_analysis.R) | Downstream analysis: binned LD vs allele frequency across the genome |

---

## Main pipeline: `clone_freebayes.sh`

### Workflow

**Per sample** (for every `*_R1_001.fastq.gz` / `*_R2_001.fastq.gz` pair in `READS`):

1. **Trimming** with BBDuk: adapter removal (`k=23`, `hdist=1`) and quality trimming of both ends (`trimq=10`, `minlen=20`).
2. **Alignment** to the Houghton reference genome (`GCA_905310635.1`) with `bwa mem -T 135`, then sorting and indexing with samtools.
3. **Mapping QC** with Qualimap `bamqc`.
4. **BAM processing** with Picard: `AddOrReplaceReadGroups`, then `MarkDuplicates`.

**Combined** (all samples together):

5. **Joint variant calling** with FreeBayes in haploid mode (`-p 1`) with minimum mapping quality 20, base quality 20 and coverage 3.
6. **Filtering**:
   - bcftools: remove multi-allelic sites (`--max-alleles 2`).
   - bcftools: remove sites with a missing-genotype rate of about 15% or more (`AN < 17`).
   - `var_filter_freebayes.py`: remove sites that do not distinguish the two parents.
7. **Annotation** with SnpEff at three levels of stringency (coding only, genic, all), each indexed with GATK `IndexFeatureFile`.
8. **Parent-of-origin split**: each annotated VCF goes through `var_split_freebayes.py` to give Weybridge-derived and Wisconsin-derived variant sets.

### Requirements

- SLURM (the script requests 64 CPUs and about 150 GB RAM)
- A micromamba environment (`sample_vc`) containing `bwa`, `samtools`, `qualimap`, `freebayes`, `bcftools` and `python3`
- BBTools (`bbduk.sh` and its `adapters.fa`)
- Picard, SnpEff (with an *E. tenella* database) and GATK
- An indexed *E. tenella* Houghton reference FASTA

### Configuration

Edit the `SET PATHS` block at the top of the script before running:

| Variable | Description |
|---|---|
| `READS` | Folder of raw paired-end FASTQ files |
| `OUTPATH` | Working/output directory |
| `REF` | Indexed reference FASTA |
| `SCRIPTS` | Folder containing the companion Python scripts |
| `BBDUK`, `ADPTR` | BBDuk executable and adapter file |
| `PICARD` | Picard JAR |
| `SNPEFF`, `SNPEFFDATABASE` | SnpEff JAR and database name |
| `GATK` | GATK executable |
| `MFH` | `MAX_FILE_HANDLES` for MarkDuplicates (check with `ulimit -n`) |
| `THREADS` | Should match `--cpus-per-task` |

The Java steps assume up to 64 GB of heap (`-Xmx64g`). Change this to suit the memory you have available.

### Usage

```bash
sbatch clone_freebayes.sh
```

### Outputs

All outputs are written under `OUTPATH`:

```
OUTPATH/
├── PerSample/
│   └── <ID>/
│       ├── Trimming/            # <ID>_clean1.fastq, <ID>_clean2.fastq
│       ├── Alignments/          # <ID>VsRef.bam, <ID>VsRef_RG.bam,
│       │                        # <ID>VsRef_MD.bam (+ .bai), <ID>VsRef_MD.metrics
│       ├── <ID>_Stats/Qualimap/ # Qualimap BAM QC report
│       └── <ID>_Logs/           # BBDuk and Picard logs
├── Combined/
│   ├── bamlist.txt              # List of deduplicated BAMs passed to FreeBayes
│   ├── VariantCalling/
│   │   └── combined.vcf         # Raw joint FreeBayes calls
│   ├── FilteredVars/
│   │   ├── combined_nomulti.vcf # Biallelic sites only
│   │   ├── combined_nomiss.vcf  # + missingness filter
│   │   └── combined_pyfilt.vcf  # + parent-discriminating sites only
│   ├── AnnotatedVars/
│   │   ├── combined_coding.ann.vcf  # SnpEff, coding effects only
│   │   ├── combined_genic.ann.vcf   # SnpEff, genic (incl. introns/UTRs)
│   │   └── combined_all.ann.vcf     # SnpEff, all effects
│   └── SplitVariants/
│       ├── WeybridgeVars/       # WeyCodingAnnot.vcf, WeyGenicAnnot.vcf, WeyAllAnnot.vcf
│       └── WisconsinVars/       # WisCodingAnnot.vcf, WisGenicAnnot.vcf, WisAllAnnot.vcf
├── Stats/                       # SnpEff summary reports (SNP_stats_*.html)
└── Logs/                        # FreeBayes log
```

---

## Companion scripts

Both scripts assume the two parental strains are the **first two sample columns** of the VCF (0-based columns 9 and 10: Weybridge, then Wisconsin). The pipeline calls them automatically, but they can also be run on their own.

### `var_filter_freebayes.py`

```bash
python3 var_filter_freebayes.py <input.vcf> <output.vcf>
```

Keeps only variants that **distinguish the two parents**, meaning one parent carries REF and the other carries ALT. A site is discarded if:

- both parents share the same allele,
- either parent has a missing genotype (`.`), or
- either parent carries a higher-order alternate allele (2–5).

Header lines are copied through unchanged.

### `var_split_freebayes.py`

```bash
python3 var_split_freebayes.py <input.vcf>
```

Splits a VCF by parent of origin based on which parental column carries the ALT allele (`1`). It writes `WeyVars.vcf` and `WisVars.vcf` to the current directory, and both files keep the full header. The main pipeline renames and moves these into `SplitVariants/`.

---

## Downstream analysis: `linkage_af_analysis.R`

This script looks for genomic regions where **linkage disequilibrium and allele-frequency distortion occur together** among the recombinant clones. Such regions are candidate loci linked to the selected phenotypes. Because every retained site is biallelic between the two parents, the per-site frequency is the parental allele ratio across the haploid clones, not a population minor-allele frequency.

### Inputs

Generated from the filtered VCF with PLINK/VCFtools (not produced by `clone_freebayes.sh`):

| File | Source |
|---|---|
| `vcftools_test.ld.reformat` | Pairwise r² table (`CHR_A BP_A SNP_A CHR_B BP_B SNP_B R2`) |
| `clone_vars_freq.frq` | Per-site allele frequency (`plink --freq`) |
| `all_vars.bim` | PLINK `.bim` file, used to get SNP positions |

File names and the bin size are set in **STEP 0** of the script.

### Method

1. Assign every site to a genomic bin (`bin_size`, set to 20 kb by default).
2. **LD per bin**: each SNP pair counts towards the bin of *both* its SNPs, so long-range LD shows up in every bin it touches. The mean r² is then taken per bin.
3. **Allele frequency per bin**: mean A1 frequency per bin.
4. Merge the two with a full outer join, so bins with only one of the two metrics are kept (this makes gaps in marker coverage visible).
5. Compute a **candidate score**: each metric is normalised to its genome-wide maximum and the two are multiplied. Bins that are high in both rank highest. This score is a way to rank bins for follow-up, not a statistical test.

### Outputs

- **`ld_freq_binned.csv`**: per-chromosome, per-bin table with `mean_r2`, `n_pairs`, `mean_af`, `n_sites`, the normalised metrics and `candidate_score`. Sort by `candidate_score` to rank candidate regions.
- **Plots** (shown on the active graphics device; save them with `ggsave()` if needed):
  - **Chr5 and Chr12 close-up**: mean LD and mean allele frequency stacked against genomic position, with the candidate regions highlighted in red.
  - **Genome-wide LD**: mean r² per bin, one row per chromosome (Chr1–Chr15).
  - **Genome-wide allele frequency**: mean allele frequency per bin, one row per chromosome.
  - **Faceted overview**: all nuclear chromosomes with a 0.9 threshold line.

Apicoplast and mitochondrial contigs (`HG994976.1`, `HG994977.1`) are left out of the genome-wide plots, and GenBank accessions are relabelled as Chr1–Chr15.

### Requirements

R with `data.table`, `ggplot2` and `tidyverse`.

```bash
Rscript linkage_af_analysis.R
```
