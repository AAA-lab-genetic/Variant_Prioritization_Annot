"""
individual_level_pipeline.py

Step 2 of the variant pipeline: individual-level (carrier/genotype) analysis,
built on top of Step 1's variant-level annotation table
(variant_annotation_pipeline.py).

This is a set of composable functions, not a rigid fixed-column script --
master-file column names and the ancestry/label list are all passed in,
since those vary across projects. A CLI + orchestrator at the bottom chains
the common path, with the ancestry-MAF, double-het, and monogenic steps all
optional (matching "it also has an option to ..." in the original workflow).

Pipeline stages
---------------
1. plink19_recode_raw()       -- PLINK1.9 --recode A on the (same) input
                                  plink fileset -> wide 0/1/2 genotype .raw.
                                  Uses PLINK1.9 specifically, not plink2.
2. build_master_metadata()    -- read a master/keyfile CSV, keep only the
                                  columns you specify, rename its ID column(s)
                                  to FID/IID, optionally filter rows (e.g.
                                  wgs == 1). Column list is caller-supplied,
                                  so a different master file layout just
                                  means different keep_cols/id_rename args.
3. merge_raw_with_master()    -- left-join the wide .raw onto the master
                                  metadata on IID.
4. split_by_label()           -- split the merged wide table into one file
                                  per ancestry/label value.
5. carrier_long_format()      -- per label: unpivot to long format, keep
                                  only carriers (Genotype > 0), strip the
                                  trailing allele suffix off variant names,
                                  and join in the Step 1 variant annotation
                                  table.

Optional stages
----------------
6. run_ancestry_case_maf()    -- PLINK1.9 --filter-cases --recode A on a
                                  per-ancestry plink fileset (parallelized
                                  across gene x label, mirroring the
                                  original bash "&"+"wait" pattern).
7. summarize_variant_counts() -- turn a --recode A .raw into per-variant
                                  N_OBS/MISSING/HET/HOM_ALT/ALT_AF/MAF
                                  counts for that ancestry.
8. merge_annotation_with_counts() -- left-join the Step 1 annotation table
                                  with every ancestry's count file on
                                  Variants, producing ONE combined
                                  "variant info" reference table (this is
                                  the file your original snippets referred
                                  to as "...multianno...counts.txt" -- it's
                                  now something this pipeline produces,
                                  not an external mystery file).
9. find_double_het()          -- flag individuals carrying >1 qualifying
                                  (rare + not benign) variant.
10. add_monogenic_data()      -- left-join a GP2-style monogenic variant
                                  report onto the carrier table.
11. combine_carrier_files()   -- concatenate all per-ancestry carrier files
                                  into one final "carrier info" file.

Together, merge_annotation_with_counts() (variant-level) and
combine_carrier_files() (individual-level) are the two final outputs
described as "one carrier info file and variant info file merged as one".
"""

from __future__ import annotations

import argparse
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import polars as pl

# Fixed columns PLINK's --recode A always writes, in order, before the
# per-variant genotype columns (after FID is dropped in merge_raw_with_master).
PLINK_RAW_FIXED_COLS = ["IID", "PAT", "MAT", "SEX", "PHENOTYPE"]

# Default final column order/selection for the combined carrier file,
# applied automatically by run_individual_level_pipeline() /
# the CLI unless overridden (pass final_column_order=None, or
# --no-final-column-order on the CLI, to skip reordering entirely).
# Columns not present in a given run (e.g. Mono_* if add_monogenic=False,
# or MAF_PD_<label> if compute_ancestry_maf=False) are filled with nulls
# by reorder_columns() rather than raising.
DEFAULT_FINAL_COLUMN_ORDER = [
    "FID", "IID", "SEX", "diagnosis", "PHENOTYPE", "GP2_phenotype_for_qc",
    "biological_sex_for_qc", "age_at_sample_collection", "age_of_onset",
    "age_at_diagnosis", "age_at_death", "age_at_last_follow_up",
    "race_for_qc", "family_history_for_qc", "doublet_het_vars",
    "doublet_het_status", "wgs_label", "region_for_qc",
    "var", "Genotype", "Variants",
    "Mono_Variant_Name", "Mono_Zygosity", "Mono_Data_type",
    "Mono_Pathogenicity", "Mono_Variant_interpretation", "Mono_comment",
    "Chr", "Start", "End", "Ref", "Alt", "avsnp147",
    "Func.refGene", "Gene.refGene", "GeneDetail.refGene",
    "ExonicFunc.refGene", "AAChange.refGene", "Franklin", "CLNSIG",
    "Consequence", "IMPACT", "CADD_phred", "CADD_PHRED", "CADD_RAW",
    "EVE_score", "AlphaMissense_score", "AlphaMissense_pred",
    "SIFT_score", "SIFT_pred", "SIFT4G_score", "SIFT4G_pred",
    "Polyphen2_HVAR_score", "Polyphen2_HVAR_pred", "REVEL_score",
    "gnomad41_genome_fafmax_faf95_max",
    "gnomad41_genome_AF_afr", "gnomad41_genome_AF_ami",
    "gnomad41_genome_AF_amr", "gnomad41_genome_AF_asj",
    "gnomad41_genome_AF_eas", "gnomad41_genome_AF_fin",
    "gnomad41_genome_AF_mid", "gnomad41_genome_AF_nfe",
    "gnomad41_genome_AF_remaining", "gnomad41_genome_AF_sas",
    "GTEx_V8_eQTL_gene", "GTEx_V8_eQTL_tissue",
    "GTEx_V8_sQTL_gene", "GTEx_V8_sQTL_tissue", "eQTLGen_snp_id",
] + [
    f"{stat}_{label}"
    for label in ["AAC", "AFR", "AJ", "AMR", "CAH", "CAS", "EAS", "EUR", "MDE", "SAS"]
    for stat in ["N_OBS", "MISSING", "HET", "HOM_ALT", "ALT_AF", "MAF_PD"]
]


# ---------------------------------------------------------------------------
# Subprocess helper (shared style with variant_annotation_pipeline.py)
# ---------------------------------------------------------------------------

def _run(cmd: list[str], **kwargs) -> None:
    print(f"[run] {' '.join(str(c) for c in cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True, **kwargs)
    combined = result.stdout + result.stderr
    lines = combined.splitlines()
    max_lines = 300
    if len(lines) <= max_lines:
        if combined.strip():
            print(combined)
    else:
        print(f"[run] output truncated ({len(lines)} lines total) -- showing first 20 and last 20:")
        print("\n".join(lines[:20]))
        print(f"... [{len(lines) - 40} lines omitted] ...")
        print("\n".join(lines[-20:]))
    if result.returncode != 0:
        result.check_returncode()


# ---------------------------------------------------------------------------
# 1. PLINK1.9 recode to .raw
# ---------------------------------------------------------------------------

def plink19_recode_raw(
    bfile_prefix: str,
    out_prefix: str,
    plink19_path: str = "plink1.9",
    filter_cases: bool = False,
    extra_args: list[str] | None = None,
) -> str:
    """
    Run PLINK1.9 `--recode A` on a plink1 fileset (.bed/.bim/.fam prefix),
    producing a wide 0/1/2 (missing=NA) genotype .raw file. Deliberately
    PLINK1.9, not plink2 -- plink2's --export A output format differs.

    Parameters
    ----------
    filter_cases : if True, adds --filter-cases (subsets to PHENOTYPE==2
        i.e. cases only, per the .fam/pheno file) before recoding -- used
        for the case-only ancestry MAF calculation (stage 6).
    extra_args : any additional PLINK1.9 flags to pass through as-is.

    Returns
    -------
    str : path to the resulting f"{out_prefix}.raw"
    """
    cmd = [plink19_path, "--bfile", bfile_prefix]
    if filter_cases:
        cmd += ["--filter-cases"]
    if extra_args:
        cmd += extra_args
    cmd += ["--recode", "A", "--out", out_prefix]
    _run(cmd)
    return f"{out_prefix}.raw"


# ---------------------------------------------------------------------------
# 2. Master/keyfile metadata
# ---------------------------------------------------------------------------

def build_master_metadata(
    master_file: str,
    keep_cols: list[str],
    id_rename: dict[str, str] | None = None,
    filter_col: str | None = None,
    filter_value=None,
) -> pl.DataFrame:
    """
    Read a master/keyfile CSV and reduce it to a standardized metadata
    table. Column names and which columns to keep are entirely
    caller-supplied, so a differently-laid-out master file just means
    different keep_cols/id_rename/filter_col arguments -- nothing here is
    hardcoded to one project's column names.

    Parameters
    ----------
    keep_cols : columns to retain, using the master file's OWN column
        names (e.g. ["FID", "GP2ID", "wgs_label", "diagnosis", ...]).
    id_rename : optional {original_name: standardized_name} mapping,
        applied AFTER keep_cols selection (e.g. {"GP2ID": "IID"}). The
        renamed column becomes the join key used by merge_raw_with_master().
    filter_col, filter_value : optional row filter, e.g.
        filter_col="wgs", filter_value=1 to keep only wgs==1 rows. Leave
        both None to keep every row. filter_value is coerced to match the
        column's actual dtype before comparing (so filter_value=1,
        filter_value="1", and filter_value="1.0" all work against a
        numeric column, e.g. when passed in as a string from the CLI).

    Returns
    -------
    pl.DataFrame : the reduced, renamed metadata table.
    """
    master = pl.read_csv(master_file, infer_schema_length=10000)
    if filter_col is not None:
        col_dtype = master.schema[filter_col]
        try:
            if col_dtype in (pl.Float32, pl.Float64):
                typed_value = float(filter_value)
            elif col_dtype in (pl.Int8, pl.Int16, pl.Int32, pl.Int64,
                                pl.UInt8, pl.UInt16, pl.UInt32, pl.UInt64):
                typed_value = int(float(filter_value))
            else:
                typed_value = filter_value
            master = master.filter(pl.col(filter_col) == typed_value)
        except (ValueError, TypeError):
            # fall back to comparing as strings if the value truly can't
            # be coerced to the column's dtype
            master = master.filter(pl.col(filter_col).cast(pl.Utf8) == str(filter_value))
    master = master.select(keep_cols)
    if id_rename:
        master = master.rename(id_rename)
    return master


# ---------------------------------------------------------------------------
# 3. Merge .raw with master metadata
# ---------------------------------------------------------------------------

def merge_raw_with_master(
    raw_file: str,
    master_df: pl.DataFrame,
    id_col: str = "IID",
    drop_cols: tuple[str, ...] = ("FID",),
) -> tuple[pl.DataFrame, list[str]]:
    """
    Left-join a PLINK --recode A .raw file onto master metadata on id_col.

    Returns
    -------
    (merged_df, meta_cols) : the merged wide table, and the explicit list
        of non-variant ("metadata") column names in it -- master_df's
        columns plus PLINK's fixed columns (PAT/MAT/SEX/PHENOTYPE). This
        list is computed from actual column names, not a positional slice,
        so it stays correct regardless of how many metadata columns the
        master file has.
    """
    raw = pl.read_csv(raw_file, separator=" ", null_values="NA", infer_schema_length=10000)
    raw = raw.drop([c for c in drop_cols if c in raw.columns])
    merged = master_df.join(raw, on=id_col, how="left")

    meta_cols = list(master_df.columns)
    for c in PLINK_RAW_FIXED_COLS:
        if c not in meta_cols and c in merged.columns:
            meta_cols.append(c)
    return merged, meta_cols


# ---------------------------------------------------------------------------
# 4. Split by ancestry/label
# ---------------------------------------------------------------------------

def split_by_label(
    df: pl.DataFrame,
    label_col: str,
    out_dir: str,
    out_prefix: str,
    labels: list[str] | None = None,
) -> dict[str, str]:
    """
    Split df into one TSV per value of label_col.

    labels : explicit list of values to emit (e.g. the 10 GP2 ancestry
        codes). If None, every distinct non-null value present in
        label_col is used instead.

    Returns
    -------
    dict[label, path] for the files written.
    """
    if labels is None:
        labels = sorted(df[label_col].drop_nulls().unique().to_list())
    Path(out_dir).mkdir(parents=True, exist_ok=True)

    paths = {}
    for label in labels:
        tmp = df.filter(pl.col(label_col) == label)
        out_file = str(Path(out_dir) / f"{out_prefix}_{label}.raw")
        tmp.write_csv(out_file, separator="\t", include_header=True)
        paths[label] = out_file
    return paths


# ---------------------------------------------------------------------------
# 5. Long-format carrier table + annotation join
# ---------------------------------------------------------------------------

def carrier_long_format(
    label_raw_file: str,
    annotation_file: str,
    meta_cols: list[str],
    out_file: str,
    variant_suffix_pattern: str = r"_[^_]+$",
    join_key: str = "Variants",
) -> str:
    """
    Unpivot a wide per-label .raw (metadata cols + one col per variant) to
    long format, keep only carriers (Genotype > 0), strip the trailing
    allele suffix off the variant column name (e.g. "chr6:123:C:T_T" ->
    "chr6:123:C:T"), and left-join the Step 1 variant annotation table.

    meta_cols : explicit non-variant column names (from
        merge_raw_with_master's second return value) -- everything else in
        label_raw_file is treated as a variant genotype column.

    Returns
    -------
    str : out_file path.
    """
    raw = pl.read_csv(label_raw_file, separator="\t", null_values="NA", infer_schema_length=10000)
    value_vars = [c for c in raw.columns if c not in meta_cols]

    long = raw.unpivot(index=meta_cols, on=value_vars, variable_name="var", value_name="Genotype")
    long = long.with_columns(pl.col("Genotype").cast(pl.Float64, strict=False))
    long = long.filter(pl.col("Genotype") > 0)
    long = long.with_columns(pl.col("var").str.replace(variant_suffix_pattern, "").alias(join_key))

    annotation = pl.read_csv(annotation_file, separator="\t", null_values=["NA", "."], infer_schema_length=10000)
    results = long.join(annotation, on=join_key, how="left")
    results.write_csv(out_file, separator="\t", include_header=True)
    return out_file


# ---------------------------------------------------------------------------
# 6. Optional: per-ancestry case-only recode (for MAF-in-cases calculation)
# ---------------------------------------------------------------------------

def run_ancestry_case_maf(
    genes: list[str],
    labels: list[str],
    bfile_template: str,
    plink19_path: str = "plink1.9",
    filter_cases: bool = True,
    max_workers: int = 8,
) -> dict[tuple[str, str], str]:
    """
    Run plink19_recode_raw() for every (gene, label) combination, in
    parallel (thread pool -- these are I/O/subprocess-bound, mirroring the
    original bash "command & ... wait" fan-out per ancestry).

    bfile_template : a str.format() template for each ancestry's plink1
        fileset prefix, using {gene} and {label} placeholders, e.g.
        "data/wgs/{label}/{gene}/{label}_release12_vwb_{gene}_updated"
        (no .bed/.bim/.fam extension).

    Returns
    -------
    dict[(gene, label), raw_path] for every case-only .raw produced.
    """
    tasks = [(gene, label) for label in labels for gene in genes]

    def _one(gene: str, label: str) -> tuple[tuple[str, str], str]:
        bfile = bfile_template.format(gene=gene, label=label)
        out_prefix = f"{bfile}.case"
        raw_path = plink19_recode_raw(bfile, out_prefix, plink19_path, filter_cases=filter_cases)
        return (gene, label), raw_path

    results: dict[tuple[str, str], str] = {}
    with ThreadPoolExecutor(max_workers=max_workers) as ex:
        futures = [ex.submit(_one, g, l) for g, l in tasks]
        for f in as_completed(futures):
            key, path = f.result()
            results[key] = path
    return results


# ---------------------------------------------------------------------------
# 7. Optional: per-ancestry variant count / MAF summary
# ---------------------------------------------------------------------------

def summarize_variant_counts(
    case_raw_file: str,
    label: str,
    out_file: str,
    id_cols: tuple[str, ...] = ("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"),
    variant_suffix_pattern: str = r"_[^_]+$",
) -> str:
    """
    Summarize a PLINK1.9 --recode A .raw (typically case-only) into
    per-variant N_OBS / MISSING / HET / HOM_ALT / ALT_AF / MAF, columns
    suffixed with the ancestry label, keyed by Variants.

    Returns
    -------
    str : out_file path.
    """
    raw = pl.read_csv(case_raw_file, separator=" ", null_values="NA", infer_schema_length=10000)
    variant_cols = [c for c in raw.columns if c not in id_cols]

    long = raw.select(variant_cols).unpivot(variable_name="VAR", value_name="geno")
    counts = long.group_by("VAR").agg([
        pl.col("geno").is_not_null().sum().alias("N_OBS"),
        pl.col("geno").is_null().sum().alias("MISSING"),
        (pl.col("geno") == 1).sum().alias("HET"),
        (pl.col("geno") == 2).sum().alias("HOM_ALT"),
    ])
    counts = counts.with_columns(
        ((pl.col("HET") + 2 * pl.col("HOM_ALT")) / (2 * pl.col("N_OBS"))).alias("ALT_AF")
    )
    counts = counts.with_columns(
        pl.min_horizontal([pl.col("ALT_AF"), 1 - pl.col("ALT_AF")]).alias("MAF")
    )
    counts = counts.with_columns(
        pl.col("VAR").str.replace(variant_suffix_pattern, "").alias("Variants")
    )
    counts = counts.select(["Variants", "N_OBS", "MISSING", "HET", "HOM_ALT", "ALT_AF", "MAF"])
    counts = counts.rename({
        "N_OBS": f"N_OBS_{label}", "MISSING": f"MISSING_{label}",
        "HET": f"HET_{label}", "HOM_ALT": f"HOM_ALT_{label}",
        "ALT_AF": f"ALT_AF_{label}", "MAF": f"MAF_PD_{label}",
    })
    counts.write_csv(out_file, separator="\t", include_header=True)
    return out_file


def summarize_all_variant_counts(
    genes: list[str],
    labels: list[str],
    case_raw_template: str,
    out_template: str,
) -> dict[tuple[str, str], str]:
    """
    Run summarize_variant_counts() for every (gene, label) combination.

    case_raw_template, out_template : str.format() templates using
        {gene}/{label}, pointing at each case-only .raw (from
        run_ancestry_case_maf) and where to write each count file.

    Returns
    -------
    dict[(gene, label), out_path]
    """
    out = {}
    for gene in genes:
        for label in labels:
            case_raw = case_raw_template.format(gene=gene, label=label)
            out_file = out_template.format(gene=gene, label=label)
            Path(out_file).parent.mkdir(parents=True, exist_ok=True)
            summarize_variant_counts(case_raw, label, out_file)
            out[(gene, label)] = out_file
    return out


# ---------------------------------------------------------------------------
# 8. Combine Step 1 annotation + all ancestry counts -> ONE variant info file
# ---------------------------------------------------------------------------

def merge_annotation_with_counts(
    annotation_file: str,
    count_files: dict[str, str],
    out_file: str,
    join_key: str = "Variants",
) -> str:
    """
    Left-join the Step 1 variant annotation table with every ancestry's
    count file (from summarize_variant_counts) on join_key, producing one
    combined variant-level reference table (annotation + per-ancestry
    MAF_PD_<label> etc. all in one row per variant).

    count_files : {label: path}, e.g. output of summarize_all_variant_counts
        filtered/collapsed to one file per label (if summarize_all_variant_counts
        was run per-gene, pass only that gene's files here).

    Returns
    -------
    str : out_file path.
    """
    base = pl.read_csv(annotation_file, separator="\t", null_values=["NA", "."], infer_schema_length=10000)
    for _, path in count_files.items():
        counts = pl.read_csv(path, separator="\t", null_values=["NA", "."], infer_schema_length=10000)
        base = base.join(counts, on=join_key, how="left")
    base.write_csv(out_file, separator="\t", include_header=True)
    return out_file


# ---------------------------------------------------------------------------
# 9. Optional: double-het flagging
# ---------------------------------------------------------------------------

def find_double_het(
    carrier_file: str,
    annotation_counts_file: str,
    out_file: str,
    freq_col: str = "gnomad41_genome_fafmax_faf95_max",
    freq_threshold: float = 0.001,
    pathogenicity_col: str = "Franklin",
    exclude_values: tuple[str, ...] = ("Likely Benign", "Benign"),
    join_key: str = "Variants",
) -> str:
    """
    Flag individuals carrying more than one qualifying variant (rare --
    freq_col null or < freq_threshold -- AND not in exclude_values on
    pathogenicity_col), and attach the combined "Variant(Genotype), ..."
    string plus a doublet_het_status column.

    NOTE on join direction: this replicates the original script's join
    order exactly -- the double-het individuals table is joined as the
    LEFT side onto the full carrier table, so the OUTPUT ONLY CONTAINS
    ROWS FOR DOUBLE-HET INDIVIDUALS (each of their carried variants, via
    fan-out), not the full carrier table with a status column added. If
    you actually want every carrier row preserved with doublet_het_status
    defaulting to "No"/"." for everyone else, swap the join direction
    (results.join(grouped, ...) instead of grouped.join(results, ...)) --
    flag if you'd like that version instead.

    Returns
    -------
    str : out_file path.
    """
    results = pl.read_csv(carrier_file, separator="\t", null_values=["NA", "."], infer_schema_length=10000)
    annot = pl.read_csv(annotation_counts_file, separator="\t", null_values=["NA", "."], infer_schema_length=10000)
    annot = annot.with_columns(pl.col(freq_col).cast(pl.Float64, strict=False))

    qualifying = annot.filter(
        (pl.col(freq_col).is_null() | (pl.col(freq_col) < freq_threshold))
        & (~pl.col(pathogenicity_col).is_in(list(exclude_values)))
    ).select(join_key)

    double_het = results.filter(
        pl.col("IID").is_duplicated() & pl.col(join_key).is_in(qualifying[join_key])
    ).select(["IID", join_key, "Genotype"])
    double_het = double_het.with_columns(
        (pl.col(join_key) + "(" + pl.col("Genotype").cast(pl.Utf8) + ")").alias("Variant_Genotype")
    )

    grouped = (
        double_het.group_by("IID")
        .agg(pl.col("Variant_Genotype").str.join(", ").alias("doublet_het_vars"))
        .with_columns(pl.lit("Yes").alias("doublet_het_status"))
    )

    merged = grouped.join(results, on="IID", how="left").fill_null(".")
    merged.write_csv(out_file, separator="\t", include_header=True)
    return out_file


# ---------------------------------------------------------------------------
# 10. Optional: monogenic variant report merge
# ---------------------------------------------------------------------------

def add_monogenic_data(
    carrier_file: str,
    variant_report_file: str,
    out_file: str,
    report_keep_cols: list[str] | None = None,
    report_id_rename: dict[str, str] | None = None,
    final_column_order: list[str] | None = None,
) -> str:
    """
    Left-join a GP2-style monogenic variant report onto the carrier table
    on IID.

    report_keep_cols : which columns to keep from variant_report_file,
        using ITS own column names. Defaults to the GP2 variant-report
        columns used in the original workflow if not given:
        ["GP2ID", "Variant_Name", "Zygosity", "Data_type",
         "Pathogenicity", "Variant_interpretation", "comment"]
    report_id_rename : {original: renamed} mapping applied after
        selection. Defaults to GP2ID->IID plus a "Mono_" prefix on the
        other report columns, matching the original workflow, if not given.
    final_column_order : optional explicit column order/selection for the
        output. If None, all columns from the merge are kept, in whatever
        order the join produced (carrier columns first, then the Mono_*
        columns) -- pass an explicit list if you want the original
        project's fixed 80-ish-column layout reproduced exactly.

    Returns
    -------
    str : out_file path.
    """
    if report_keep_cols is None:
        report_keep_cols = ["GP2ID", "Variant_Name", "Zygosity", "Data_type",
                             "Pathogenicity", "Variant_interpretation", "comment"]
    if report_id_rename is None:
        report_id_rename = {
            "GP2ID": "IID", "Variant_Name": "Mono_Variant_Name",
            "Zygosity": "Mono_Zygosity", "Data_type": "Mono_Data_type",
            "Pathogenicity": "Mono_Pathogenicity",
            "Variant_interpretation": "Mono_Variant_interpretation",
            "comment": "Mono_comment",
        }

    report = pl.read_csv(variant_report_file, infer_schema_length=10000)
    report = report.select(report_keep_cols).rename(report_id_rename)

    carriers = pl.read_csv(carrier_file, separator="\t", null_values=["NA", "."], infer_schema_length=10000)
    merged = carriers.join(report, on="IID", how="left")

    if final_column_order is not None:
        merged = merged.select(final_column_order)

    merged.write_csv(out_file, separator="\t", include_header=True)
    return out_file


# ---------------------------------------------------------------------------
# 11. Final combine: all-ancestry carrier file
# ---------------------------------------------------------------------------

def combine_carrier_files(carrier_files: dict[str, str] | list[str], out_file: str) -> str:
    """
    Concatenate every per-ancestry carrier file into one final carrier-info
    file. Uses a relaxed diagonal concat so files with slightly different
    columns (e.g. some ancestries had no double-het/monogenic matches and
    thus differently-typed null columns) still combine cleanly.

    Returns
    -------
    str : out_file path.
    """
    paths = list(carrier_files.values()) if isinstance(carrier_files, dict) else carrier_files
    dfs = [pl.read_csv(p, separator="\t", null_values=["NA", "."], infer_schema_length=10000) for p in paths]
    combined = pl.concat(dfs, how="diagonal_relaxed")
    combined.write_csv(out_file, separator="\t", include_header=True)
    return out_file


def reorder_columns(in_file: str, column_order: list[str], out_file: str) -> str:
    """
    Select/reorder columns in a TSV to an explicit final layout. Any column
    in column_order that isn't present in the file (e.g. Mono_* columns
    when add_monogenic wasn't run, or MAF_PD_<label> columns when
    compute_ancestry_maf wasn't run) is added as an all-null column rather
    than raising -- so one fixed target column list works regardless of
    which optional pipeline stages were actually run. Columns present in
    the file but NOT in column_order are dropped.

    Returns
    -------
    str : out_file path.
    """
    df = pl.read_csv(in_file, separator="\t", null_values=["NA", "."], infer_schema_length=10000)
    missing = [c for c in column_order if c not in df.columns]
    if missing:
        df = df.with_columns([pl.lit(None).alias(c) for c in missing])
    df = df.select(column_order)
    df.write_csv(out_file, separator="\t", include_header=True)
    return out_file


def split_by_variant_function(
    carrier_file: str,
    out_dir: str,
    file_prefix: str,
    function_col: str = "Func.refGene",
) -> dict[str, str]:
    """
    Split a carrier file into one file per distinct value of function_col
    (e.g. Func.refGene: "exonic", "intronic", "UTR3", ...), writing each to
    f"{out_dir}/{file_prefix}_{safe_value}.txt". Null values in
    function_col are skipped (not written to any file). Non-alphanumeric
    characters in the value are replaced with "_" to make a filesystem-safe
    filename (e.g. "splicing;exonic" -> "splicing_exonic").

    Returns
    -------
    dict[function_value, path] for every file written.
    """
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    df = pl.read_csv(carrier_file, separator="\t", null_values=["NA", "."], infer_schema_length=10000)

    paths = {}
    for vtype in df[function_col].unique().to_list():
        if vtype is None:
            continue
        safe_vtype = "".join(c if c.isalnum() or c == "_" else "_" for c in vtype)
        subset = df.filter(pl.col(function_col) == vtype)
        out_file = str(Path(out_dir) / f"{file_prefix}_{safe_vtype}.txt")
        subset.write_csv(out_file, separator="\t", include_header=True)
        paths[vtype] = out_file
    return paths


# ---------------------------------------------------------------------------
# Orchestrator
# ---------------------------------------------------------------------------

def run_individual_level_pipeline(
    bfile_prefix: str,
    gene: str,
    work_dir: str,
    master_file: str,
    master_keep_cols: list[str],
    master_id_rename: dict[str, str],
    label_col: str,
    labels: list[str],
    annotation_file: str,
    plink19_path: str = "plink1.9",
    master_filter_col: str | None = None,
    master_filter_value=None,
    # -- optional: ancestry case-only MAF --
    compute_ancestry_maf: bool = False,
    ancestry_bfile_template: str | None = None,
    genes_for_maf: list[str] | None = None,
    # -- optional: double het --
    compute_double_het: bool = False,
    double_het_freq_col: str = "gnomad41_genome_fafmax_faf95_max",
    double_het_freq_threshold: float = 0.001,
    double_het_pathogenicity_col: str = "Franklin",
    double_het_exclude_values: tuple[str, ...] = ("Likely Benign", "Benign"),
    # -- optional: monogenic --
    add_monogenic: bool = False,
    monogenic_report_file: str | None = None,
    monogenic_final_column_order: list[str] | None = None,
    # -- optional: split by variant function --
    split_by_function: bool = False,
    function_col: str = "Func.refGene",
    function_split_out_dir: str | None = None,
    # -- final output naming/location --
    out_dir: str | None = None,
    out_name: str | None = None,
    variant_out_name: str | None = None,
    final_column_order: list[str] | None = DEFAULT_FINAL_COLUMN_ORDER,
) -> dict[str, str]:
    """
    Run the full individual-level pipeline for one gene: recode -> merge
    with master -> split by ancestry -> (optional) ancestry case-only MAF
    -> carrier long format (joined against the MAF-enriched annotation
    table, if computed) -> (optional) double het -> (optional) monogenic
    -> (optional) split by variant function -> combine -> final column
    reorder (applied by default -- see final_column_order below).

    out_dir : directory the FINAL combined carrier file (and, if
        compute_ancestry_maf, the final variant info file) are written to.
        Defaults to work_dir if not given. work_dir still holds all
        per-ancestry intermediate files regardless.
    out_name : base filename (no ext) for the final combined carrier file;
        written to f"{out_dir}/{out_name}.txt". Defaults to
        f"{gene}_carrier_info_combined" if not given.
    variant_out_name : base filename (no ext) for the final variant info
        file (only produced if compute_ancestry_maf=True); written to
        f"{out_dir}/{variant_out_name}.txt". Defaults to
        f"{gene}_variant_info" if not given.
    split_by_function : if True, split each label's fully-annotated
        carrier file into one file per distinct function_col value (e.g.
        "exonic", "intronic", ...) using split_by_variant_function(). This
        runs on the per-label files AFTER double het/monogenic (whichever
        ran), so the split files carry all annotation already joined in.
        It does NOT change what goes into the combined carrier file --
        it's an additional, separate breakdown.
    function_col : column to split on (default "Func.refGene").
    function_split_out_dir : where the split files are written. Defaults
        to f"{work_dir}/by_variant_type" if not given.
    final_column_order : explicit column list/order applied to
        the FINAL combined carrier file (after all other stages,
        including combine_carrier_files). Defaults to
        DEFAULT_FINAL_COLUMN_ORDER (the project's standard layout) --
        pass final_column_order=None explicitly to skip reordering and
        leave columns in whatever order the pipeline produced, or pass
        your own list to use a different fixed layout. Any column in the
        list that isn't present in a given run (e.g. Mono_* columns when
        add_monogenic=False) is added as an all-null column so the
        selection doesn't fail -- so one layout works regardless of which
        optional stages were run.

    Returns
    -------
    dict with keys "carrier_files" (dict[label, final per-label path]),
    "combined_carrier_file", and, if compute_ancestry_maf, "variant_info_file",
    and, if split_by_function, "split_by_function_files"
    (dict[label, dict[function_value, path]]).
    """
    work_dir_p = Path(work_dir)
    work_dir_p.mkdir(parents=True, exist_ok=True)
    out_dir_p = Path(out_dir) if out_dir else work_dir_p
    out_dir_p.mkdir(parents=True, exist_ok=True)
    out_name = out_name or f"{gene}_carrier_info_combined"
    variant_out_name = variant_out_name or f"{gene}_variant_info"

    # 1. recode
    raw_prefix = str(work_dir_p / f"{gene}_updated")
    raw_file = plink19_recode_raw(bfile_prefix, raw_prefix, plink19_path)

    # 2. master metadata
    master_df = build_master_metadata(
        master_file, master_keep_cols, master_id_rename,
        master_filter_col, master_filter_value,
    )

    # 3. merge
    merged_df, meta_cols = merge_raw_with_master(raw_file, master_df)

    # 4. split by ancestry
    label_paths = split_by_label(merged_df, label_col, str(work_dir_p), f"{gene}_updated", labels)

    result: dict[str, str] = {}

    # 5 (was 6-8, moved earlier): optional ancestry case-only MAF ->
    # combined variant info file. This MUST run before carrier_long_format
    # below, so the carrier table gets joined against the MAF-enriched
    # annotation table (not just the plain Step 1 annotation) -- otherwise
    # the MAF_PD_<label> etc. columns never make it into the carrier file.
    if compute_ancestry_maf:
        assert ancestry_bfile_template is not None, "ancestry_bfile_template required when compute_ancestry_maf=True"
        genes_for_maf = genes_for_maf or [gene]
        run_ancestry_case_maf(genes_for_maf, labels, ancestry_bfile_template, plink19_path)
        count_template = str(work_dir_p / "{label}_{gene}_var_counts.txt")
        case_raw_template = ancestry_bfile_template + ".case.raw"
        count_paths = summarize_all_variant_counts(genes_for_maf, labels, case_raw_template, count_template)
        this_gene_counts = {label: path for (g, label), path in count_paths.items() if g == gene}
        variant_info_file = str(out_dir_p / f"{variant_out_name}.txt")
        merge_annotation_with_counts(annotation_file, this_gene_counts, variant_info_file)
        result["variant_info_file"] = variant_info_file
        annotation_file = variant_info_file  # carrier_long_format below uses the enriched table

    # 6 (was 5): carrier long format, per label -- now joined against the
    # MAF-enriched annotation_file if compute_ancestry_maf ran above.
    carrier_files = {}
    for label, path in label_paths.items():
        out_file = str(work_dir_p / f"{gene}_{label}_carrier_annot.txt")
        carrier_long_format(path, annotation_file, meta_cols, out_file)
        carrier_files[label] = out_file

    result["carrier_files"] = dict(carrier_files)

    # 9. optional: double het (per label)
    if compute_double_het:
        for label, path in list(carrier_files.items()):
            out_file = str(work_dir_p / f"{gene}_{label}_carrier_annot.doublehet.txt")
            find_double_het(
                path, annotation_file, out_file,
                freq_col=double_het_freq_col, freq_threshold=double_het_freq_threshold,
                pathogenicity_col=double_het_pathogenicity_col, exclude_values=double_het_exclude_values,
            )
            carrier_files[label] = out_file

    # 10. optional: monogenic
    if add_monogenic:
        assert monogenic_report_file is not None, "monogenic_report_file required when add_monogenic=True"
        for label, path in list(carrier_files.items()):
            out_file = str(work_dir_p / f"{gene}_{label}_carrier_annot.final.txt")
            add_monogenic_data(path, monogenic_report_file, out_file,
                                final_column_order=monogenic_final_column_order)
            carrier_files[label] = out_file

    # 10b. optional: split by variant function, per label (does not affect
    # the combined carrier file below -- an additional, separate breakdown)
    if split_by_function:
        split_out_dir = function_split_out_dir or str(work_dir_p / "by_variant_type")
        split_files: dict[str, dict[str, str]] = {}
        for label, path in carrier_files.items():
            split_files[label] = split_by_variant_function(
                path, split_out_dir, f"{gene}_{label}", function_col=function_col,
            )
        result["split_by_function_files"] = split_files

    # 11. combine across ancestries
    combined_file = str(out_dir_p / f"{out_name}.txt")
    combine_carrier_files(carrier_files, combined_file)
    result["carrier_files"] = dict(carrier_files)

    # 12. optional: final column reorder, applied to the combined file
    if final_column_order is not None:
        reorder_columns(combined_file, final_column_order, combined_file)

    result["combined_carrier_file"] = combined_file

    return result


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _cli():
    p = argparse.ArgumentParser(description="Step 2: individual-level carrier/genotype pipeline.")
    p.add_argument("--bfile-prefix", required=True, help="plink1 fileset prefix (.bed/.bim/.fam, no extension)")
    p.add_argument("--gene", required=True)
    p.add_argument("--work-dir", required=True)
    p.add_argument("--master-file", required=True)
    p.add_argument("--master-keep-cols", required=True, help="Comma-separated list of columns to keep from the master file")
    p.add_argument("--master-id-rename", help="Comma-separated old:new pairs, e.g. 'GP2ID:IID'")
    p.add_argument("--master-filter-col")
    p.add_argument("--master-filter-value")
    p.add_argument("--label-col", required=True, help="Ancestry/label column name in the master file")
    p.add_argument("--labels", required=True, help="Comma-separated list of labels to split by")
    p.add_argument("--annotation-file", required=True, help="Step 1 variant annotation table")
    p.add_argument("--plink19-path", default="plink1.9")

    p.add_argument("--compute-ancestry-maf", action="store_true")
    p.add_argument("--ancestry-bfile-template", help="e.g. 'data/wgs/{label}/{gene}/{label}_release12_vwb_{gene}_updated'")
    p.add_argument("--genes-for-maf", help="Comma-separated gene list for the ancestry MAF step (default: just --gene)")

    p.add_argument("--compute-double-het", action="store_true")
    p.add_argument("--double-het-freq-col", default="gnomad41_genome_fafmax_faf95_max")
    p.add_argument("--double-het-freq-threshold", type=float, default=0.001)
    p.add_argument("--double-het-pathogenicity-col", default="Franklin")
    p.add_argument("--double-het-exclude-values", default="Likely Benign,Benign")

    p.add_argument("--add-monogenic", action="store_true")
    p.add_argument("--monogenic-report-file")

    p.add_argument("--split-by-function", action="store_true", help="Split each label's final carrier file into one file per Func.refGene value (or --function-col)")
    p.add_argument("--function-col", default="Func.refGene")
    p.add_argument("--function-split-out-dir", help="Where split files are written. Default: <work-dir>/by_variant_type")

    p.add_argument("--out-dir", help="Directory for the FINAL combined carrier file (and variant info file, if --compute-ancestry-maf). Default: work-dir")
    p.add_argument("--out-name", help="Base filename (no ext) for the FINAL combined carrier file. Default: <gene>_carrier_info_combined")
    p.add_argument("--variant-out-name", help="Base filename (no ext) for the FINAL variant info file. Default: <gene>_variant_info")
    p.add_argument("--final-column-order", help="Comma-separated explicit column list/order for the FINAL combined carrier file, overriding the built-in default layout (DEFAULT_FINAL_COLUMN_ORDER). Missing columns (e.g. Mono_* if --add-monogenic wasn't used) are added as all-null rather than erroring.")
    p.add_argument("--no-final-column-order", action="store_true", help="Skip column reordering entirely and leave columns in whatever order the pipeline produced (default layout is applied otherwise).")

    args = p.parse_args()

    keep_cols = [c.strip() for c in args.master_keep_cols.split(",")]
    id_rename = None
    if args.master_id_rename:
        id_rename = dict(pair.split(":") for pair in args.master_id_rename.split(","))
    labels = [l.strip() for l in args.labels.split(",")]
    genes_for_maf = [g.strip() for g in args.genes_for_maf.split(",")] if args.genes_for_maf else None
    exclude_values = tuple(v.strip() for v in args.double_het_exclude_values.split(","))
    if args.no_final_column_order:
        final_column_order = None
    elif args.final_column_order:
        final_column_order = [c.strip() for c in args.final_column_order.split(",")]
    else:
        final_column_order = DEFAULT_FINAL_COLUMN_ORDER

    result = run_individual_level_pipeline(
        bfile_prefix=args.bfile_prefix, gene=args.gene, work_dir=args.work_dir,
        master_file=args.master_file, master_keep_cols=keep_cols, master_id_rename=id_rename,
        label_col=args.label_col, labels=labels, annotation_file=args.annotation_file,
        plink19_path=args.plink19_path,
        master_filter_col=args.master_filter_col, master_filter_value=args.master_filter_value,
        compute_ancestry_maf=args.compute_ancestry_maf,
        ancestry_bfile_template=args.ancestry_bfile_template, genes_for_maf=genes_for_maf,
        compute_double_het=args.compute_double_het,
        double_het_freq_col=args.double_het_freq_col,
        double_het_freq_threshold=args.double_het_freq_threshold,
        double_het_pathogenicity_col=args.double_het_pathogenicity_col,
        double_het_exclude_values=exclude_values,
        add_monogenic=args.add_monogenic, monogenic_report_file=args.monogenic_report_file,
        split_by_function=args.split_by_function, function_col=args.function_col,
        function_split_out_dir=args.function_split_out_dir,
        out_dir=args.out_dir, out_name=args.out_name, variant_out_name=args.variant_out_name,
        final_column_order=final_column_order,
    )
    print(f"Combined carrier file: {result['combined_carrier_file']}")
    if "variant_info_file" in result:
        print(f"Variant info file: {result['variant_info_file']}")
    if "split_by_function_files" in result:
        n_files = sum(len(v) for v in result["split_by_function_files"].values())
        split_dir = args.function_split_out_dir or str(Path(args.work_dir) / "by_variant_type")
        print(f"Split-by-function: {n_files} files written to {split_dir}")
        for label, files in result["split_by_function_files"].items():
            print(f"  {label}: {list(files.keys())}")


if __name__ == "__main__":
    _cli()
