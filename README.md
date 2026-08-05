# Variant file annotation preparation pipeline (Step 1: Variant Annotation, Step 2: Individual-Level Carrier Analysis)

Two-step pipeline for gene-level variant work:

- **Step 1** (`variant_annotation_pipeline.py`) — Two stages variant-level annotation:
  ANNOVAR + VEP (with CADD), optionally merged with a Franklin export. One
  row per variant.

Output example:

```Folder XXX```

- **Step 2** (`individual_level_pipeline.py`) — individual-level carrier
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

Omit `--compute-ancestry-maf`/`--compute-double-het`/`--add-monogenic` (or
their Python equivalents) to skip those stages entirely. If you skip
`--compute-ancestry-maf`, the double-het/monogenic steps run against your
original `annotation_file` instead of the enriched variant-info file —
make sure it has whatever columns `find_double_het()` needs
(`gnomad41_genome_fafmax_faf95_max`, `Franklin` by default; both are
overridable via `--double-het-freq-col`/`--double-het-pathogenicity-col`).

## Step 2 output naming

- `work_dir` — where every per-ancestry **intermediate** file is written
  (the `.raw`, the per-label splits, the per-label carrier/doublehet
  files before the final combine).
- `out_dir` — where the two **final** outputs go. Defaults to `work_dir`
  if not given.
- `out_name` — base filename for the final combined carrier file.
  Defaults to `{gene}_carrier_info_combined`.
- `variant_out_name` — base filename for the final variant info file
  (only produced if `compute_ancestry_maf=True`). Defaults to
  `{gene}_variant_info`.
  
**Core path (always runs):**

1. `plink19_recode_raw()` — PLINK1.9 `--recode A` on the input plink1
   fileset → wide 0/1/2 (missing=NA) genotype `.raw`. Deliberately
   PLINK1.9, not plink2 — plink2's `--export A` output format differs.
2. `build_master_metadata()` — read a master/keyfile CSV, keep only the
   columns you specify, rename its ID column(s) to `FID`/`IID`,
   optionally filter rows (e.g. `wgs == 1`).
3. `merge_raw_with_master()` — left-join the wide `.raw` onto the master
   metadata on `IID`.
4. `split_by_label()` — split the merged wide table into one file per
   ancestry/label value.
5. `carrier_long_format()` — per label: unpivot to long format, keep only
   carriers (`Genotype > 0`), strip the trailing allele suffix off
   variant names, and join in the Step 1 variant annotation table.

**Optional stages** (each toggled independently — matches "it also has an
option to..." from the original workflow):

6. `run_ancestry_case_maf()` — PLINK1.9 `--filter-cases --recode A` on a
   per-ancestry plink fileset, parallelized across gene x label (mirrors
   the original bash `command & ... wait` fan-out).
7. `summarize_variant_counts()` — turn a `--recode A` `.raw` into
   per-variant `N_OBS`/`MISSING`/`HET`/`HOM_ALT`/`ALT_AF`/`MAF` for one
   ancestry.
8. `merge_annotation_with_counts()` — left-join the Step 1 annotation
   table with every ancestry's count file on `Variants`, producing **one**
   combined variant-level reference table (annotation + every ancestry's
   MAF, one row per variant). This is the file your original workflow
   referred to as `...multianno...counts.txt` — it's now something this
   pipeline produces itself, not an external file you have to build by hand.
9. `find_double_het()` — flag individuals carrying more than one
   qualifying (rare + not benign) variant.
10. `add_monogenic_data()` — left-join a GP2-style monogenic variant
    report onto the carrier table.
11. `combine_carrier_files()` — concatenate all per-ancestry carrier files
    into one final carrier-info file.

Stages 8 and 11 are the two final outputs described as "one carrier info
file and variant info file merged as one" — one variant-level table, one
carrier-level table.

## Master file flexibility

Different projects have different master/keyfile layouts, so nothing here
is hardcoded to one column set:

- `keep_cols` — pass whichever columns exist in *your* master file, using
  its own column names.
- `id_rename` — maps your ID column(s) to the standardized names the rest
  of the pipeline expects (e.g. `{"GP2ID": "IID"}`).
- `filter_col`/`filter_value` — optional row filter (e.g. `wgs == 1`).
  `filter_value` is coerced to match the column's actual dtype before
  comparing, so `1`, `"1"`, and `"1.0"` all work correctly against a
  numeric column — this matters because CLI arguments always arrive as
  strings, and master files often store flag columns as floats.

The metadata column list used downstream (`meta_cols`, needed to know
which columns are "not a variant" when unpivoting) is computed from
`merge_raw_with_master()`'s actual output columns — not a positional
slice/index — so it stays correct regardless of how many metadata columns
your master file has.


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
