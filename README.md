# Variant file annotation preparation pipeline (Step 1: Variant Annotation, Step 2: Individual-Level Carrier Analysis)

Two-step pipeline for gene-level variant work:
![Uploading annotation_flow.png…]()


- **Step 1** (`variant_annotation_pipeline.py`) - Two stages variant-level annotation:
  ANNOVAR + VEP (with CADD), optionally merged with a Franklin export. One
  row per variant.

Output example:

```Folder XXX```

- **Step 2** (`individual_level_pipeline.py`) - individual-level carrier
  analysis built on top of Step 1's output: who carries which variant,
  per-ancestry carrier frequencies, double-het flagging, monogenic report
  merging. One row per carrier (individual x variant), split by ancestry.

Output example:

```Folder XXX```


---

# Step 1 — Variant Annotation Pipeline

Implemented as two Python functions in `variant_annotation_pipeline.py`,
plus a CLI that chains them.

1. Optional region extraction (`chrom`/`start`/`end`).
2. Optional MAF filter (`maf_min`/`maf_max`) — computed on **all** samples
   in the input file.
3. Sample subsetting to `n_samples` (default 10) or a supplied
   `sample_list` — for annotation speed only.

## Step 1 CLI usage

```bash
python variant_annotation_pipeline.py \
  --input-file data/wgs/merged/BLTP3A/merged_release12_vwb \
  --file-type plink2 \
  --gene BLTP3A \
  --work-dir data/wgs/merged/BLTP3A \
  --chrom 6 --start 34700000 --end 34900000 \
  --maf-max 0.01 \
  --n-samples 10 \
  --subset-out-dir data/wgs/merged/BLTP3A \
  --subset-out-name BLTP3A_subset \
  --annovar-pl /home/jupyter/annovar/table_annovar.pl \
  --annovar-db /home/jupyter/annovar/humandb/ \
  --vep-cache-dir /home/jupyter/.vep/ \
  --vep-plugin-dir /home/jupyter/.vep/ \
  --cadd-snv-file /home/jupyter/.vep/CADD_v1.7_GRCh38/whole_genome_SNVs.tsv.gz \
  --cadd-indel-file /home/jupyter/.vep/CADD_v1.7_GRCh38/gnomad.genomes.r4.0.indel.tsv.gz \
  --franklin-file results/processed/wgs/annotation/BLTP3A_r12_franklin.txt \
  --out-dir results/processed/wgs/annotation \
  --out-name BLTP3A_final_annotation \
  --bcftools-path /opt/conda/envs/vep_env/bin/bcftools \
  --plink-path /opt/conda/envs/vep_env/bin/plink2 \
  --bgzip-path /opt/conda/envs/vep_env/bin/bgzip \
  --tabix-path /opt/conda/envs/vep_env/bin/tabix
```

```
usage: vars_annot.py [-h] --input-file INPUT_FILE --file-type
                     {plink1,plink2,vcf} --gene GENE --work-dir WORK_DIR
                     [--chrom CHROM] [--start START] [--end END]
                     [--maf-min MAF_MIN] [--maf-max MAF_MAX]
                     [--n-samples N_SAMPLES] [--sample-list SAMPLE_LIST]
                     [--subset-out-dir SUBSET_OUT_DIR]
                     [--subset-out-name SUBSET_OUT_NAME] --annovar-pl
                     ANNOVAR_PL --annovar-db ANNOVAR_DB --vep-cache-dir
                     VEP_CACHE_DIR --vep-plugin-dir VEP_PLUGIN_DIR
                     --cadd-snv-file CADD_SNV_FILE --cadd-indel-file
                     CADD_INDEL_FILE [--franklin-file FRANKLIN_FILE]
                     [--genome-build GENOME_BUILD]
                     [--vep-cache-version VEP_CACHE_VERSION]
                     [--out-dir OUT_DIR] [--out-name OUT_NAME]
                     [--plink-path PLINK_PATH] [--bcftools-path BCFTOOLS_PATH]
                     [--bgzip-path BGZIP_PATH] [--tabix-path TABIX_PATH]

Two-stage variant annotation pipeline.

options:
  -h, --help            show this help message and exit
  --input-file INPUT_FILE
  --file-type {plink1,plink2,vcf}
  --gene GENE
  --work-dir WORK_DIR   Where ANNOVAR/VEP intermediate files are written
  --chrom CHROM
  --start START
  --end END
  --maf-min MAF_MIN     Minimum MAF -- keeps variants >= this (drops rare),
                        computed on all samples before subsetting
  --maf-max MAF_MAX     Maximum MAF -- keeps variants <= this (drops common),
                        e.g. 0.01 to keep only rare variants; computed on all
                        samples before subsetting
  --n-samples N_SAMPLES
  --sample-list SAMPLE_LIST
  --subset-out-dir SUBSET_OUT_DIR
                        Directory for the Stage 1 subset VCF (default: work-
                        dir)
  --subset-out-name SUBSET_OUT_NAME
                        Base filename (no ext) for the Stage 1 subset VCF
                        (default: <gene>.subset)
  --annovar-pl ANNOVAR_PL
  --annovar-db ANNOVAR_DB
  --vep-cache-dir VEP_CACHE_DIR
  --vep-plugin-dir VEP_PLUGIN_DIR
  --cadd-snv-file CADD_SNV_FILE
  --cadd-indel-file CADD_INDEL_FILE
  --franklin-file FRANKLIN_FILE
  --genome-build GENOME_BUILD
  --vep-cache-version VEP_CACHE_VERSION
  --out-dir OUT_DIR     Directory for the FINAL annotation table (default:
                        work-dir)
  --out-name OUT_NAME   Base filename (no ext) for the FINAL annotation table
                        (default: <gene>_annotation_final)
  --plink-path PLINK_PATH
                        Path to the plink2 binary (default: 'plink2', resolved
                        via PATH)
  --bcftools-path BCFTOOLS_PATH
                        Path to the bcftools binary. If bcftools misbehaves
                        (e.g. 'invalid option' errors on modern flags like
                        -O/-o), PATH is probably resolving to a stale/wrong
                        binary -- pass the absolute path here instead, e.g.
                        output of `which bcftools` run inside your activated
                        env.
  --bgzip-path BGZIP_PATH
                        Path to the bgzip binary (default: 'bgzip', resolved
                        via PATH)
  --tabix-path TABIX_PATH
                        Path to the tabix binary (default: 'tabix', resolved
                        via PATH)
```

- Omit `--chrom/--start/--end` if the input file is already region-subset.
- Omit `--maf-min`/`--maf-max` to skip MAF filtering in that direction --
`--maf-min` drops rare variants below the cutoff, `--maf-max` drops common
variants above it (use `--maf-max` alone to keep only rare variants).
- `--subset-out-dir`/`--subset-out-name` and `--out-dir`/`--out-name` are all
optional — omit any of them to fall back to `--work-dir` and gene-based
default names.


### Franklin file option `--franklin-file`

This file can be obtained by running each variants using automate Franklin web crawl that we have developed (https://github.com/AAA-lab-genetic/frank_claw)



### Tool paths (`--bcftools-path` etc.)

All four tool-path flags (`--plink-path`, `--bcftools-path`,
`--bgzip-path`, `--tabix-path`) default to the bare command name resolved
via `PATH`, and are optional. 


```
[tools path default]
  plink: 'plink2' -> /PATH/TO/plink2
  bcftools: 'bcftools' -> /PATH/TO/bcftools
  bgzip: 'bgzip' -> /PATH/TO/bgzip
  tabix: 'tabix' -> /PATH/TO/tabix
```






---

# Step 2 — Individual-Level Carrier/Genotype Pipeline

Step 2 is dependent on Step 1, where
this step works at the **individual/genotype level**: who carries which
variant, per-ancestry carrier frequencies, double-het flagging, and
monogenic report merging.


## Step 2 pipeline stages

## Step 2 usage — CLI

```bash
python individual_level_pipeline.py \
  --bfile-prefix data/wgs/merged/BLTP3A/merged_release12_vwb_BLTP3A_updated \
  --gene BLTP3A \
  --work-dir data/wgs/merged/BLTP3A \
  --master-file data/wgs/master_key_release12_final_vwb.csv \
  --master-keep-cols "FID,GP2ID,wgs_label,diagnosis,GP2_phenotype_for_qc,biological_sex_for_qc,age_at_sample_collection,age_of_onset,age_at_diagnosis,age_at_death,age_at_last_follow_up,race_for_qc,family_history_for_qc,region_for_qc" \
  --master-id-rename "GP2ID:IID" \
  --master-filter-col wgs \
  --master-filter-value 1 \
  --label-col wgs_label \
  --labels "AAC,AFR,AJ,AMR,CAH,CAS,EAS,EUR,MDE,SAS" \
  --annotation-file results/processed/wgs/annotation/BLTP3A_final_annotation.txt \
  --plink19-path /home/jupyter/plink1.9 \
  --compute-ancestry-maf \
  --ancestry-bfile-template "data/wgs/{label}/BLTP3A/{label}_release12_vwb_{gene}_updated" \
  --compute-double-het \
  --add-monogenic \
  --monogenic-report-file /home/jupyter/workspace/gp2_tier2_eu_release12/variant_report/GP2_release12_variant_report.csv \
  --out-dir results/processed/wgs/annotation \
  --out-name BLTP3A_carrier_info_combined \
  --variant-out-name BLTP3A_variant_info
```

```
usage: individual_level_pipeline.py [-h] --bfile-prefix BFILE_PREFIX --gene
                                    GENE --work-dir WORK_DIR --master-file
                                    MASTER_FILE --master-keep-cols
                                    MASTER_KEEP_COLS
                                    [--master-id-rename MASTER_ID_RENAME]
                                    [--master-filter-col MASTER_FILTER_COL]
                                    [--master-filter-value MASTER_FILTER_VALUE]
                                    --label-col LABEL_COL --labels LABELS
                                    --annotation-file ANNOTATION_FILE
                                    [--plink19-path PLINK19_PATH]
                                    [--compute-ancestry-maf]
                                    [--ancestry-bfile-template ANCESTRY_BFILE_TEMPLATE]
                                    [--genes-for-maf GENES_FOR_MAF]
                                    [--compute-double-het]
                                    [--double-het-freq-col DOUBLE_HET_FREQ_COL]
                                    [--double-het-freq-threshold DOUBLE_HET_FREQ_THRESHOLD]
                                    [--double-het-pathogenicity-col DOUBLE_HET_PATHOGENICITY_COL]
                                    [--double-het-exclude-values DOUBLE_HET_EXCLUDE_VALUES]
                                    [--add-monogenic]
                                    [--monogenic-report-file MONOGENIC_REPORT_FILE]
                                    [--split-by-function]
                                    [--function-col FUNCTION_COL]
                                    [--function-split-out-dir FUNCTION_SPLIT_OUT_DIR]
                                    [--out-dir OUT_DIR] [--out-name OUT_NAME]
                                    [--variant-out-name VARIANT_OUT_NAME]
                                    [--final-column-order FINAL_COLUMN_ORDER]
                                    [--no-final-column-order]

Step 2: individual-level carrier/genotype pipeline.

options:
  -h, --help            show this help message and exit
  --bfile-prefix BFILE_PREFIX
                        plink1 fileset prefix (.bed/.bim/.fam, no extension)
  --gene GENE
  --work-dir WORK_DIR
  --master-file MASTER_FILE
  --master-keep-cols MASTER_KEEP_COLS
                        Comma-separated list of columns to keep from the
                        master file
  --master-id-rename MASTER_ID_RENAME
                        Comma-separated old:new pairs, e.g. 'GP2ID:IID'
  --master-filter-col MASTER_FILTER_COL
  --master-filter-value MASTER_FILTER_VALUE
  --label-col LABEL_COL
                        Ancestry/label column name in the master file
  --labels LABELS       Comma-separated list of labels to split by
  --annotation-file ANNOTATION_FILE
                        Step 1 variant annotation table
  --plink19-path PLINK19_PATH
  --compute-ancestry-maf
  --ancestry-bfile-template ANCESTRY_BFILE_TEMPLATE
                        e.g. 'data/wgs/{label}/{gene}/{label}_release12_vwb_{g
                        ene}_updated'
  --genes-for-maf GENES_FOR_MAF
                        Comma-separated gene list for the ancestry MAF step
                        (default: just --gene)
  --compute-double-het
  --double-het-freq-col DOUBLE_HET_FREQ_COL
  --double-het-freq-threshold DOUBLE_HET_FREQ_THRESHOLD
  --double-het-pathogenicity-col DOUBLE_HET_PATHOGENICITY_COL
  --double-het-exclude-values DOUBLE_HET_EXCLUDE_VALUES
  --add-monogenic
  --monogenic-report-file MONOGENIC_REPORT_FILE
  --split-by-function   Split each label's final carrier file into one file
                        per Func.refGene value (or --function-col)
  --function-col FUNCTION_COL
  --function-split-out-dir FUNCTION_SPLIT_OUT_DIR
                        Where split files are written. Default: <work-
                        dir>/by_variant_type
  --out-dir OUT_DIR     Directory for the FINAL combined carrier file (and
                        variant info file, if --compute-ancestry-maf).
                        Default: work-dir
  --out-name OUT_NAME   Base filename (no ext) for the FINAL combined carrier
                        file. Default: <gene>_carrier_info_combined
  --variant-out-name VARIANT_OUT_NAME
                        Base filename (no ext) for the FINAL variant info
                        file. Default: <gene>_variant_info
  --final-column-order FINAL_COLUMN_ORDER
                        Comma-separated explicit column list/order for the
                        FINAL combined carrier file, overriding the built-in
                        default layout (DEFAULT_FINAL_COLUMN_ORDER). Missing
                        columns (e.g. Mono_* if --add-monogenic wasn't used)
                        are added as all-null rather than erroring.
  --no-final-column-order
                        Skip column reordering entirely and leave columns in
                        whatever order the pipeline produced (default layout
                        is applied otherwise).
```

- Omit `--compute-ancestry-maf`/`--compute-double-het`/`--add-monogenic` (or
their Python equivalents) to skip those stages entirely.

- If you skip
`--compute-ancestry-maf`, the double-het/monogenic steps run against your
original `annotation_file` instead of the enriched variant-info file —
make sure it has whatever columns `find_double_het()` needs
(`gnomad41_genome_fafmax_faf95_max`, `Franklin` by default; both are
overridable via `--double-het-freq-col`/`--double-het-pathogenicity-col`).

## Input needed
### Per-gene input data
| File | Used by | Notes |
|---|---|---|
| `merged_release12_vwb_{gene}_updated.bed/.bim/.fam` | Step 1 (`--input-file`), Step 2 (`--bfile-prefix`) | full-cohort plink1 fileset for the gene region, all samples |
| `{label}_release12_vwb_{gene}_updated.bed/.bim/.fam` (one set per ancestry label) | Step 2 optional ancestry-MAF step (`--ancestry-bfile-template`) | only needed if using `--compute-ancestry-maf` |

### Reference/metadata files
| File | Used by | Notes |
|---|---|---|
| `master_key_release12_final_vwb.csv` | Step 2 (`--master-file`) | sample metadata + ancestry label column; column names are configurable via `--master-keep-cols`/`--master-id-rename` |
| `GP2_release12_variant_report.csv` | Step 2 optional (`--monogenic-report-file`) | only needed if using `--add-monogenic` |
| `{gene}_r12_franklin.txt` | Step 1 optional (`--franklin-file`) | only needed if merging Franklin annotations |
| `BLTP_fam_list.txt` (gene, chr, start, end) | not read directly by either script — your own driver loop over genes | needed if you're looping the pipeline across multiple genes |

## Directory structure

```
BLTP_fam/                                  # WORK_DIR
├── analysis/
│   ├── variant_annotation_pipeline.py
│   ├── individual_level_pipeline.py
│   └── export_tool_versions.py
│
├── data/
│   ├── BLTP_fam_list.txt                  # gene  chr  start  end  (your driver list, if looping genes)
│   ├── wgs/
│   │   ├── master_key_release12_final_vwb.csv
│   │   ├── merged/
│   │   │   └── {gene}/                    # e.g. BLTP3A/
│   │   │       ├── merged_release12_vwb_{gene}_updated.bed/.bim/.fam        # INPUT (yours)
│   │   │       ├── merged_release12_vwb_{gene}_updated_maf005.bed/.bim/.fam # INPUT (yours), alt. pre-MAF-filtered version
│   │   │       │
│   │   │       ├── {gene}_subset.vcf.gz(.tbi)            # Step 1 intermediate
│   │   │       ├── {gene}.annovar.*                      # Step 1 intermediate (ANNOVAR outputs)
│   │   │       ├── {gene}.VEP.vcf                         # Step 1 intermediate
│   │   │       │
│   │   │       ├── {gene}_updated.raw                     # Step 2 intermediate (plink1.9 recode)
│   │   │       ├── {gene}_updated_{label}.raw              # Step 2 intermediate, per ancestry
│   │   │       ├── {gene}_{label}_carrier_annot*.txt        # Step 2 intermediate, per ancestry
│   │   │       │
│   │   │       └── by_variant_type/                        # Step 2 optional split-by-function output
│   │   │           └── {gene}_{label}_{func}.txt
│   │   │
│   │   └── {label}/                        # per-ancestry plink filesets, e.g. EUR/, AFR/, ...
│   │       └── {gene}/
│   │           ├── {label}_release12_vwb_{gene}_updated.bed/.bim/.fam       # INPUT (yours), for ancestry MAF
│   │           ├── {label}_release12_vwb_{gene}_updated.case.raw            # Step 2 intermediate (case-only recode)
│   │           └── {label}_release12_vwb_{gene}_var_counts.txt             # Step 2 intermediate (per-ancestry counts)
│   │
│   └── ... (other genes' subfolders, same pattern)
│
└── results/
    └── processed/
        └── wgs/
            └── annotation/
                ├── {gene}_r12_franklin.txt                      # INPUT (yours), optional
                ├── {gene}_final_annotation.txt                   # Step 1 FINAL output
                ├── {gene}_variant_info.txt                       # Step 2 FINAL output (variant-level, ancestry MAF merged in)
                ├── {gene}_carrier_info_combined.txt              # Step 2 FINAL output (carrier-level, all ancestries)
                └── tool_versions.md / tool_versions.json         # reproducibility record

/home/jupyter/
├── plink1.9                                 # binary, not on PATH
├── annovar/
│   ├── table_annovar.pl
│   └── humandb/
├── .vep/
│   ├── (VEP cache, version 114)
│   └── CADD_v1.7_GRCh38/
│       ├── whole_genome_SNVs.tsv.gz
│       └── gnomad.genomes.r4.0.indel.tsv.gz
└── workspace/
    └── gp2_tier2_eu_release12/
        └── variant_report/
            └── GP2_release12_variant_report.csv   # INPUT (yours), optional (monogenic)
```


# Pipeline Tool Versions

| Tool | Version | 
|---|---|
| plink2 | v2.0.0-a.6.9LM (29 Jan 2025) | 
| plink1.9 | v1.90b6.9 (4 Mar 2019) |
| bcftools | 1.24 | 
| bgzip (htslib) | 1.23.1 |
| tabix (htslib) | 1.23.1 | 
| perl | 5.32.1 |
| VEP (ensembl-vep) | 116.0 |
| ANNOVAR (table_annovar.pl) | `$Date: 2022-08-02$` http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.latest.tar.gz | 
| python3 | 3.12.13 | 
| polars | 1.43.2 | 

| ANNOVAR annotation DB |
|---|
| refGene |
| cytoBand |
| exac03 |
| avsnp147|
| dbnsfp47a |
| clinvar_20250721 | 
| gnomad41_genome | 

| VEP annotation DB | link/source |
|---|--- |
| HGVS | ftp://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz |
| CADD | https://github.com/Ensembl/VEP_plugins/blob/release/116/CADD.pm |
