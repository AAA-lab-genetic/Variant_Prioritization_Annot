# Variant Pipeline (Step 1: Variant Annotation, Step 2: Individual-Level Carrier Analysis)

Two-stage pipeline for gene-level variant work:

- **Step 1** (`variant_annotation_pipeline.py`) — variant-level annotation:
  ANNOVAR + VEP (with CADD), optionally merged with a Franklin export. One
  row per variant.
- **Step 2** (`individual_level_pipeline.py`) — individual-level carrier
  analysis built on top of Step 1's output: who carries which variant,
  per-ancestry carrier frequencies, double-het flagging, monogenic report
  merging. One row per carrier (individual x variant), split by ancestry.

Both use **polars** for table processing.

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

Omit `--chrom/--start/--end` if the input file is already region-subset.
Omit `--maf-min`/`--maf-max` to skip MAF filtering in that direction --
`--maf-min` drops rare variants below the cutoff, `--maf-max` drops common
variants above it (use `--maf-max` alone to keep only rare variants).
`--subset-out-dir`/`--subset-out-name` and `--out-dir`/`--out-name` are all
optional — omit any of them to fall back to `--work-dir` and gene-based
default names.

### Tool paths (`--bcftools-path` etc.)

All four tool-path flags (`--plink-path`, `--bcftools-path`,
`--bgzip-path`, `--tabix-path`) default to the bare command name resolved
via `PATH`, and are optional. Pin them to absolute paths if you hit tool
resolution issues — e.g. `bcftools` rejecting modern flags like `-O`/`-o`
with an `invalid option` error means `PATH` is resolving to a stale/legacy
binary rather than the one you get by typing `bcftools` interactively
(Python's `subprocess` doesn't see shell aliases/functions, only `PATH`).
Every run prints a `[tool resolution]` block up front showing exactly
which absolute path each tool name resolved to, so you can check this
without guessing:

```
[tool resolution]
  plink: 'plink2' -> /opt/conda/envs/vep_env/bin/plink2
  bcftools: 'bcftools' -> /usr/bin/bcftools   <- suspicious if this isn't in vep_env
  bgzip: 'bgzip' -> /opt/conda/envs/vep_env/bin/bgzip
  tabix: 'tabix' -> /opt/conda/envs/vep_env/bin/tabix
```

## Step 1 requirements

- `plink2`, `bcftools`, `bgzip`, `tabix` available (absolute paths
  recommended -- see above).
- `perl` + ANNOVAR (`table_annovar.pl`, `humandb/`) and `vep` (with the CADD
  plugin data files) available, for Stage 2.
- Python: `polars`.

## Step 1 known assumptions to revisit

- **Merge key**: uses the VCF ID string (`Variants`, e.g.
  `chr6:34792087:C:T`) rather than Chr/Start/Ref/Alt, since that's what your
  original ANNOVAR reduction and Franklin merge both used. This assumes the
  same VCF ID convention is present in both the ANNOVAR-run VCF and the
  VEP-run VCF — worth double-checking if you ever swap in a differently
  formatted input VCF.
- **Multi-transcript variants**: canonical-transcript selection is a
  simplification — flag if you need multi-transcript (long-format) output
  instead.
  
## Stage 1 — `prepare_variant_subset()`

Accepts **PLINK1** (`.bed/.bim/.fam`), **PLINK2** (`.pgen/.pvar/.psam`), or
**VCF** input. Outputs a bgzipped, indexed VCF.

```python
from variant_annotation_pipeline import prepare_variant_subset

vcf = prepare_variant_subset(
    input_file="data/wgs/merged/BLTP3A/merged_release12_vwb",  # plink prefix or vcf path
    file_type="plink2",             # "plink1" | "plink2" | "vcf"
    out_dir="data/wgs/merged/BLTP3A",
    out_name="BLTP3A_subset",       # -> data/wgs/merged/BLTP3A/BLTP3A_subset.vcf.gz
    chrom="6", start=34700000, end=34900000,   # omit all 3 if input is already region-subset
    maf_min=None,                   # keep variants with MAF >= this; None = no minimum
    maf_max=0.01,                   # keep variants with MAF <= this, i.e. rare variants; None = no maximum
    n_samples=10,                   # applied AFTER the MAF filter; or sample_list="my_samples.txt"
)
```

**Notes**
- `out_dir` is created if it doesn't exist; the final VCF is written to
  `f"{out_dir}/{out_name}.vcf.gz"` (indexed alongside it).
- If `chrom/start/end` are all `None`, the input is assumed to already be
  region-subset (e.g. you extracted the gene region upstream).
- `maf_min`/`maf_max` filter in **opposite directions** — don't mix them up:
  - `maf_min=0.01` keeps variants with MAF ≥ 0.01, i.e. drops rare variants.
  - `maf_max=0.01` keeps variants with MAF ≤ 0.01, i.e. keeps **only** rare
    variants below that cutoff (set `maf_min=None` if you just want this).
  - Set both for a MAF range (e.g. `maf_min=0.001, maf_max=0.01` for a
    rare-but-not-singleton band).
  - Both are computed across the full sample set before any subsetting.
  - PLINK inputs: implemented as `plink2 --export vcf-iid` (uncompressed
    VCF, IID-only sample names), then `bgzip -f` and `tabix -f -p vcf` --
    the same tools bcftools uses internally, and the confirmed-working
    path for this setup. (An earlier version relied on plink2's own
    internal `bgz` export + `bcftools index`, which crashed intermittently;
    this replaces that.)
  - VCF input: implemented as `bcftools view -q <maf_min>:minor` /
    `-Q <maf_max>:minor`, run on the region-extracted VCF before any `-S`
    sample subset, then indexed with `tabix -f -p vcf`.
- If `sample_list` is given, it overrides `n_samples`. One VCF sample name
  per line.
- Sample IDs for `n_samples` are read directly from the VCF's `#CHROM`
  header line via Python's stdlib `gzip` module -- no `bcftools query -l`
  call, so this step doesn't depend on bcftools at all.

## Stage 2 — `run_variant_annotation()`

Runs ANNOVAR (`table_annovar.pl`) and VEP (`--plugin CADD`) on the Stage-1
VCF, reduces each output to a fixed column set, merges them, and optionally
merges in a Franklin export.

```python
from variant_annotation_pipeline import run_variant_annotation

out_file = run_variant_annotation(
    vcf_file=vcf,
    gene="BLTP3A",
    work_dir="data/wgs/merged/BLTP3A",
    annovar_pl="/home/jupyter/annovar/table_annovar.pl",
    annovar_db="/home/jupyter/annovar/humandb/",
    vep_cache_dir="/home/jupyter/.vep/",
    vep_plugin_dir="/home/jupyter/.vep/",
    cadd_snv_file="/home/jupyter/.vep/CADD_v1.7_GRCh38/whole_genome_SNVs.tsv.gz",
    cadd_indel_file="/home/jupyter/.vep/CADD_v1.7_GRCh38/gnomad.genomes.r4.0.indel.tsv.gz",
    franklin_file="results/processed/wgs/annotation/BLTP3A_r12_franklin.txt",  # optional
    out_dir="results/processed/wgs/annotation",   # optional -- defaults to work_dir
    out_name="BLTP3A_final_annotation",           # optional -- defaults to <gene>_annotation_final
)
```

`work_dir` is where the ANNOVAR/VEP intermediate files land
(`{gene}.annovar.*`, `{gene}.VEP.vcf`). `out_dir`/`out_name` control only the
**final** merged table, so you can keep intermediates in a scratch directory
and send the finished table somewhere else — e.g. a shared results folder.
Both are optional; if omitted, the final table is written to
`{work_dir}/{gene}_annotation_final.txt` as before.

### What it runs

1. **ANNOVAR**: `table_annovar.pl` with
   `refGene,cytoBand,exac03,avsnp147,dbnsfp47a,clinvar_20250721,gnomad41_genome`.
2. **Column reduction (ANNOVAR)**: keeps the same column set as your original
   `awk` script (Chr/Start/End/Ref/Alt, `Otherinfo6` → renamed **`Variants`**,
   avsnp147, refGene fields, CLNSIG, CADD/EVE/AlphaMissense/SIFT/SIFT4G/
   Polyphen2/REVEL scores, gnomad41_genome AFs, GTEx eQTL/sQTL, eQTLGen).
3. **VEP**: same command you had, `--cache_version 114`, `--plugin CADD`,
   `--everything`.
4. **Column reduction (VEP)**: parses the `CSQ` INFO field using the field
   order read directly from the VCF header (so it never breaks if VEP
   version/flags change the column order). Keeps `Consequence`, `IMPACT`,
   `SIFT`, `PolyPhen`, `CADD_PHRED`, `CADD_RAW`.
   - **Transcript selection**: a variant can hit multiple transcripts. This
     keeps the block flagged `CANONICAL=YES`; if none is flagged canonical,
     it falls back to the first transcript block. If you'd rather keep
     every transcript per variant (long format), say so and I'll add a
     `canonical_only=False` mode.
5. **Merge**: ANNOVAR-reduced ⋈ VEP-reduced on `Variants` (the VCF ID,
   e.g. `chr6:34792087:C:T`) — a left join, so every ANNOVAR-annotated
   variant is kept even if VEP's ID string doesn't match for some reason.
6. **Franklin (optional)**: if `franklin_file` is given, its `Variants`
   column is normalized (`chr` prefix stripped and re-added, matching your
   original script) and merged in. Default join is **`left`** (keeps every
   annotated variant, blank Franklin score if unmatched) — pass
   `franklin_join="inner"` to only keep variants with a Franklin match
   (matches the original script's behavior).
7. **Output**: `{work_dir}/{gene}_annotation_final.txt`, tab-separated,
   with a header row.

Row counts at each step are printed to stdout so you can sanity-check drops
at each merge, same as your original `print(df.shape[0])` calls.



---

# Step 2 — Individual-Level Carrier/Genotype Pipeline

Second stage, built on top of Step 1's variant-level annotation table.
This stage works at the **individual/genotype level**: who carries which
variant, per-ancestry carrier frequencies, double-het flagging, and
monogenic report merging. Implemented as composable functions in
`individual_level_pipeline.py`, plus an orchestrator
(`run_individual_level_pipeline()`) and a CLI that chains them.

## Why this is separate from Step 1

Step 1 produces one row per **variant**. Step 2 produces one row per
**carrier** (an individual x variant pair), split out by ancestry, with
sample metadata, per-ancestry allele frequencies, double-het status, and
monogenic findings all joined in. Different unit of analysis, different
tooling (PLINK1.9 `--recode A` specifically, not plink2), so it's a
separate module.

## Step 2 pipeline stages

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

## Step 2 usage — Python

```python
from individual_level_pipeline import run_individual_level_pipeline

result = run_individual_level_pipeline(
    bfile_prefix="data/wgs/merged/BLTP3A/merged_release12_vwb_BLTP3A_updated",
    gene="BLTP3A",
    work_dir="data/wgs/merged/BLTP3A",
    master_file="data/wgs/master_key_release12_final_vwb.csv",
    master_keep_cols=[
        "FID", "GP2ID", "wgs_label", "diagnosis", "GP2_phenotype_for_qc",
        "biological_sex_for_qc", "age_at_sample_collection", "age_of_onset",
        "age_at_diagnosis", "age_at_death", "age_at_last_follow_up",
        "race_for_qc", "family_history_for_qc", "region_for_qc",
    ],
    master_id_rename={"GP2ID": "IID"},
    master_filter_col="wgs",
    master_filter_value=1,
    label_col="wgs_label",
    labels=["AAC", "AFR", "AJ", "AMR", "CAH", "CAS", "EAS", "EUR", "MDE", "SAS"],
    annotation_file="results/processed/wgs/annotation/BLTP3A_final_annotation.txt",  # Step 1's output
    plink19_path="/home/jupyter/plink1.9",   # pin the absolute path -- see Step 1's tool-path note

    # optional: ancestry case-only MAF
    compute_ancestry_maf=True,
    ancestry_bfile_template="data/wgs/{label}/BLTP3A/{label}_release12_vwb_{gene}_updated",

    # optional: double het
    compute_double_het=True,

    # optional: monogenic
    add_monogenic=True,
    monogenic_report_file="/home/jupyter/workspace/gp2_tier2_eu_release12/variant_report/GP2_release12_variant_report.csv",

    # final output naming/location
    out_dir="results/processed/wgs/annotation",
    out_name="BLTP3A_carrier_info_combined",
    variant_out_name="BLTP3A_variant_info",
)

print(result["combined_carrier_file"])   # final all-ancestry carrier file
print(result.get("variant_info_file"))   # final all-ancestry variant file (if compute_ancestry_maf=True)
```

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

## Step 2 known assumption to confirm: double-het join direction

`find_double_het()` replicates the original script's join order exactly:
the double-het-individuals table is joined as the **left** side onto the
full carrier table (`grouped.join(results, on="IID", how="left")`). This
means the **output only contains rows for double-het individuals**
(fanned out across each variant they carry) — not the full carrier table
with a `doublet_het_status` column added and defaulting to "No"/"." for
everyone else.

If you actually want the latter (every carrier row preserved, with
double-het status as an added column), the fix is to swap the join
direction: `results.join(grouped, on="IID", how="left")` instead. I kept
the original behavior since I wasn't sure which was intended — let me
know if you want it swapped.

## Step 2 requirements

- `plink1.9` available (absolute path recommended).
- Python: `polars`.
