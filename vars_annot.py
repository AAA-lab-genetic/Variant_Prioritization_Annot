"""
variant_annotation_pipeline.py

Two-stage variant-level annotation pipeline (ANNOVAR + VEP, optional Franklin merge).

Stage 1 - prepare_variant_subset()
    Accepts a PLINK1 (.bed/.bim/.fam), PLINK2 (.pgen/.pvar/.psam), or VCF input.
    Order of operations matters here:
      1. Optional region extraction (chrom/start/end).
      2. Optional MAF filter -- computed on ALL samples in the input, since
         allele frequency is only meaningful across the full cohort. This
         happens BEFORE sample subsetting.
      3. Sample subsetting to N samples (default 10). This step exists
         purely to speed up annotation -- ANNOVAR/VEP only read variant
         records, not genotypes -- it is NOT meant to change which variants
         are present, which is why MAF filtering must happen first, on
         everyone.
    Produces a bgzipped, indexed VCF ready for annotation.

Stage 2 - run_variant_annotation()
    Runs ANNOVAR (table_annovar.pl) and VEP (with the CADD plugin) on that
    VCF, reduces each output to a fixed set of columns, merges them, and
    optionally merges in a Franklin export. Uses polars for all table work.

Both stages are plain functions -- call them from a notebook/script, or run
this file directly with the CLI at the bottom.
"""

from __future__ import annotations

import argparse
import gzip
import shutil
import subprocess
from pathlib import Path

import polars as pl

# ---------------------------------------------------------------------------
# Column sets kept from each annotation source (edit here, if needed)
# ---------------------------------------------------------------------------

ANNOVAR_KEEP_COLS = [
    "Chr", "Start", "End", "Ref", "Alt", "Otherinfo6", "avsnp147",
    "Func.refGene", "Gene.refGene", "GeneDetail.refGene",
    "ExonicFunc.refGene", "AAChange.refGene", "CLNSIG",
    "CADD_phred", "EVE_score", "AlphaMissense_score", "AlphaMissense_pred",
    "SIFT_score", "SIFT_pred", "SIFT4G_score", "SIFT4G_pred",
    "Polyphen2_HVAR_score", "Polyphen2_HVAR_pred", "REVEL_score",
    "gnomad41_genome_fafmax_faf95_max",
    "gnomad41_genome_AF_afr", "gnomad41_genome_AF_ami", "gnomad41_genome_AF_amr",
    "gnomad41_genome_AF_asj", "gnomad41_genome_AF_eas", "gnomad41_genome_AF_fin",
    "gnomad41_genome_AF_mid", "gnomad41_genome_AF_nfe",
    "gnomad41_genome_AF_remaining", "gnomad41_genome_AF_sas",
    "GTEx_V8_eQTL_gene", "GTEx_V8_eQTL_tissue",
    "GTEx_V8_sQTL_gene", "GTEx_V8_sQTL_tissue", "eQTLGen_snp_id",
]

# Otherinfo6 (the VCF ID column ANNOVAR carries through) is renamed to
# "Variants" throughout, since that's the join key against VEP and Franklin.
ANNOVAR_RENAME = {"Otherinfo6": "Variants"}

# CSQ subfields pulled out of the VEP output (index is looked up dynamically
# from the VCF header at runtime -- see _vep_csq_field_index()).
VEP_KEEP_FIELDS = ["Consequence", "IMPACT", "SIFT", "PolyPhen", "CADD_PHRED", "CADD_RAW"]


# ---------------------------------------------------------------------------
# Small subprocess helper
# ---------------------------------------------------------------------------

def _log_tool_resolution(bcftools_path: str, plink_path: str, bgzip_path: str, tabix_path: str) -> None:
    """
    Print the absolute path each tool name resolves to via PATH, before
    anything runs. If bcftools misbehaves (e.g. rejects modern flags like
    -O/-o), the resolved path here is the first thing to check -- it may
    not be the binary you get when typing `bcftools` interactively, since
    Python's subprocess does not see shell aliases/functions, only PATH.
    """
    print("[tool resolution]")
    for name, path in [("plink", plink_path), ("bcftools", bcftools_path),
                        ("bgzip", bgzip_path), ("tabix", tabix_path)]:
        resolved = shutil.which(path)
        print(f"  {name}: '{path}' -> {resolved or 'NOT FOUND ON PATH'}")

def _run(cmd: list[str], **kwargs) -> None:
    """
    Run a command, capturing stdout/stderr rather than letting it stream
    unbounded to the terminal. Normal-sized tool logs (plink2/bcftools
    progress messages, sample/variant counts, etc.) are printed in full, as
    before. If a command's output is unexpectedly huge (e.g. a `bcftools
    view` called without an output file, which would dump the entire VCF
    body to stdout) it is capped and truncated instead of flooding the
    terminal. On failure, the captured output is printed (truncated the
    same way) so the actual error is visible, then the error is re-raised.
    """
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
# Stage 1a: region + MAF filter on the FULL sample set
# ---------------------------------------------------------------------------

def _filter_region_and_maf(
    input_file: str,
    file_type: str,
    tmp_prefix: str,
    chrom: str | None,
    start: int | None,
    end: int | None,
    maf_min: float | None,
    maf_max: float | None,
    plink_path: str,
    bcftools_path: str,
    bgzip_path: str,
    tabix_path: str,
) -> str:
    """
    Apply optional region extraction and optional MAF filter(s), using every
    sample present in input_file. Returns a bgzipped, tabix-indexed,
    all-sample VCF path. This ALWAYS runs before any sample subsetting.

    maf_min keeps variants with MAF >= maf_min (drops rarer variants).
    maf_max keeps variants with MAF <= maf_max (drops more common variants).
    Set only maf_max (maf_min=None) to keep rare variants below a cutoff;
    set only maf_min to keep common variants above a cutoff; set both for
    a MAF range.
    """
    out_vcf = f"{tmp_prefix}.allsamples.vcf.gz"

    if file_type in ("plink1", "plink2"):
        cmd = [plink_path]
        cmd += ["--bfile", input_file] if file_type == "plink1" else ["--pfile", input_file]
        if chrom is not None and start is not None and end is not None:
            cmd += ["--chr", str(chrom), "--from-bp", str(start), "--to-bp", str(end)]
        if maf_min is not None:
            # computed across every sample in this plink call -- --keep is
            # deliberately NOT used here, so this happens on the full cohort.
            cmd += ["--maf", str(maf_min)]
        if maf_max is not None:
            cmd += ["--max-maf", str(maf_max)]
        # Write a plain (uncompressed) VCF -- 'vcf-iid' writes IID-only
        # sample names (no FID_IID) -- then bgzip/tabix it with the same
        # tools bcftools relies on internally. plink2's own internal `bgz`
        # export writer has been unreliable with `bcftools index` on this
        # setup; bgzip+tabix is the confirmed-working path.
        plink_out_prefix = f"{tmp_prefix}.allsamples"
        cmd += ["--export", "vcf-iid", "--out", plink_out_prefix]
        _run(cmd)
        plink_vcf = f"{plink_out_prefix}.vcf"
        _run([bgzip_path, "-f", plink_vcf])
        _run([tabix_path, "-f", "-p", "vcf", out_vcf])
        return out_vcf

    if file_type == "vcf":
        # single bcftools call, writing straight to out_vcf -- no
        # intermediate rename.
        cmd = [bcftools_path, "view"]
        if chrom is not None and start is not None and end is not None:
            cmd += ["-r", f"{chrom}:{start}-{end}"]
        if maf_min is not None:
            cmd += ["-q", f"{maf_min}:minor"]
        if maf_max is not None:
            cmd += ["-Q", f"{maf_max}:minor"]
        cmd += ["-Oz", "-o", out_vcf, input_file]
        _run(cmd)
        _run([tabix_path, "-f", "-p", "vcf", out_vcf])
        return out_vcf

    raise ValueError(f"Unknown file_type: {file_type}")


# ---------------------------------------------------------------------------
# Stage 1b: sample subsetting (speed-only -- runs on the MAF-filtered VCF)
# ---------------------------------------------------------------------------

def _vcf_sample_names(vcf_file: str) -> list[str]:
    """
    Read sample names directly from the #CHROM header line, without
    shelling out to bcftools. BGZF is valid (concatenated) gzip, so
    Python's stdlib gzip module can stream-decompress it directly -- this
    sidesteps any bcftools/htslib version mismatch on this system.
    """
    opener = gzip.open if vcf_file.endswith(".gz") else open
    with opener(vcf_file, "rt") as fh:
        for line in fh:
            if line.startswith("#CHROM"):
                fields = line.rstrip("\n").split("\t")
                return fields[9:]
            if not line.startswith("#"):
                break
    raise ValueError(f"No #CHROM header line found in {vcf_file}")


def _first_n_sample_ids(vcf_file: str, n: int) -> list[str]:
    return _vcf_sample_names(vcf_file)[:n]


def _subset_samples(
    vcf_file: str,
    output_prefix: str,
    n_samples: int | None,
    sample_list: str | None,
    bcftools_path: str,
    tabix_path: str,
) -> str:
    out_vcf = f"{output_prefix}.vcf.gz"

    if sample_list is None and n_samples is None:
        # keep everyone -- no subsetting requested. Write directly with
        # bcftools rather than os.rename(), since a rename isn't always
        # atomic/immediately visible on mounted filesystems (e.g. GCS FUSE).
        _run([bcftools_path, "view", "-Oz", "-o", out_vcf, vcf_file])
        _run([tabix_path, "-f", "-p", "vcf", out_vcf])
        return out_vcf

    if sample_list is not None:
        keep_file = sample_list
    else:
        ids = _first_n_sample_ids(vcf_file, n_samples)
        keep_file = f"{output_prefix}.keep_samples.txt"
        with open(keep_file, "w") as fh:
            for s in ids:
                fh.write(f"{s}\n")

    _run([bcftools_path, "view", "-S", keep_file, "-Oz", "-o", out_vcf, vcf_file])
    _run([tabix_path, "-f", "-p", "vcf", out_vcf])
    return out_vcf



def prepare_variant_subset(
    input_file: str,
    file_type: str,          # "plink1" | "plink2" | "vcf"
    out_dir: str,
    out_name: str,
    chrom: str | None = None,
    start: int | None = None,
    end: int | None = None,
    maf_min: float | None = None,
    maf_max: float | None = None,
    n_samples: int | None = 10,
    sample_list: str | None = None,
    plink_path: str = "plink2",
    bcftools_path: str = "bcftools",
    bgzip_path: str = "bgzip",
    tabix_path: str = "tabix",
) -> str:
    """
    Extract an optional region, apply an optional MAF filter (on every
    sample in the input), then subset to N samples purely to speed up the
    downstream annotation step. Returns the final bgzipped, indexed VCF.

    Parameters
    ----------
    input_file : path to .bed/.pgen/.vcf(.gz) (plink prefix or vcf path).
    file_type : "plink1", "plink2", or "vcf".
    out_dir : directory the output VCF (and intermediate files) are written
        to. Created if it doesn't exist.
    out_name : base filename (no extension) for the output VCF; final
        result is written to f"{out_dir}/{out_name}.vcf.gz".
    chrom, start, end : optional region to extract. If all three are None,
        the input is assumed to already be region-subset.
    maf_min : optional minimum minor allele frequency -- keeps variants
        with MAF >= maf_min (drops rarer variants). Computed across ALL
        samples in input_file, before any sample subsetting -- this is a
        variant filter, not a decision about which people to look at.
        None means no minimum applied.
    maf_max : optional maximum minor allele frequency -- keeps variants
        with MAF <= maf_max (drops more common variants). Use this alone
        (maf_min=None) to keep only rare variants below a cutoff, e.g.
        maf_max=0.01. Also computed on ALL samples before subsetting.
        Set both maf_min and maf_max for a MAF range.
    n_samples : how many samples to keep for the annotation step (default
        10). Set to None together with sample_list=None to keep all samples.
    sample_list : optional path to a file of sample IDs to keep instead of
        the first N (one VCF sample name per line). Overrides n_samples.

    Returns
    -------
    str : path to the resulting f"{out_dir}/{out_name}.vcf.gz"
    """
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    output_prefix = str(Path(out_dir) / out_name)

    _log_tool_resolution(bcftools_path, plink_path, bgzip_path, tabix_path)

    allsample_vcf = _filter_region_and_maf(
        input_file, file_type, output_prefix, chrom, start, end,
        maf_min, maf_max, plink_path, bcftools_path, bgzip_path, tabix_path,
    )
    final_vcf = _subset_samples(
        allsample_vcf, output_prefix, n_samples, sample_list, bcftools_path, tabix_path,
    )
    return final_vcf


# ---------------------------------------------------------------------------
# Stage 2: ANNOVAR + VEP annotation, column reduction, merge
# ---------------------------------------------------------------------------

def _run_annovar(vcf_file: str, gene: str, work_dir: str, annovar_pl: str,
                  annovar_db: str, genome_build: str) -> str:
    out_prefix = str(Path(work_dir) / f"{gene}.annovar")
    cmd = [
        "perl", annovar_pl, vcf_file, annovar_db,
        "-buildver", genome_build,
        "-out", out_prefix,
        "-remove",
        "-protocol", "refGene,cytoBand,exac03,avsnp147,dbnsfp47a,clinvar_20250721,gnomad41_genome",
        "-operation", "g,r,f,f,f,f,f",
        "-nopolish",
        "-nastring", ".",
        "-vcfinput",
    ]
    _run(cmd)
    return f"{out_prefix}.{genome_build}_multianno.txt"


def _reduce_annovar(multianno_file: str) -> pl.DataFrame:
    df = pl.read_csv(multianno_file, separator="\t", infer_schema_length=0)
    missing = [c for c in ANNOVAR_KEEP_COLS if c not in df.columns]
    if missing:
        raise KeyError(f"Columns missing from ANNOVAR multianno output: {missing}")
    return df.select(ANNOVAR_KEEP_COLS).rename(ANNOVAR_RENAME)


def _run_vep(vcf_file: str, gene: str, work_dir: str, vep_cache_dir: str,
             vep_plugin_dir: str, cadd_snv_file: str, cadd_indel_file: str,
             genome_build: str, vep_cache_version: str) -> str:
    out_vcf = str(Path(work_dir) / f"{gene}.VEP.vcf")
    assembly = "GRCh38" if genome_build == "hg38" else "GRCh37"
    cmd = [
        "vep",
        "--input_file", vcf_file,
        "--output_file", out_vcf,
        "--cache",
        "--cache_version", vep_cache_version,
        "--offline",
        "--assembly", assembly,
        "--dir_cache", vep_cache_dir,
        "--dir_plugins", vep_plugin_dir,
        "--plugin", f"CADD,{cadd_snv_file},{cadd_indel_file}",
        "--vcf",
        "--everything",
        "--force_overwrite",
    ]
    _run(cmd)
    return out_vcf


def _vep_csq_field_index(vep_vcf: str) -> dict[str, int]:
    """Parse the CSQ Format string from the VCF header into a {field: 0-based index} map."""
    with open(vep_vcf) as fh:
        for line in fh:
            if line.startswith("##INFO=<ID=CSQ"):
                fmt = line.split("Format: ")[1].rstrip('">\n')
                fields = fmt.split("|")
                return {name: i for i, name in enumerate(fields)}
            if not line.startswith("#"):
                break
    raise ValueError("No CSQ INFO header found in VEP output VCF.")


def _reduce_vep(vep_vcf: str, canonical_only: bool = True) -> pl.DataFrame:
    """
    Parse the VEP VCF and keep one row per variant (Chr, Start, Ref, Alt,
    Variants=ID) with the fields in VEP_KEEP_FIELDS.

    A variant can have multiple transcript annotation blocks in CSQ. By
    default this keeps the block flagged CANONICAL=YES; if no block is
    flagged canonical, it falls back to the first block for that variant.
    """
    idx = _vep_csq_field_index(vep_vcf)
    for f in VEP_KEEP_FIELDS:
        if f not in idx:
            raise KeyError(f"'{f}' not found in VEP CSQ format string.")
    canonical_idx = idx.get("CANONICAL")

    rows = []
    with open(vep_vcf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            chrom, pos, vid, ref, alt = fields[0], fields[1], fields[2], fields[3], fields[4]
            info = fields[7]
            csq_field = next((f for f in info.split(";") if f.startswith("CSQ=")), None)
            if csq_field is None:
                continue
            transcripts = csq_field[len("CSQ="):].split(",")

            chosen = None
            if canonical_only and canonical_idx is not None:
                for t in transcripts:
                    parts = t.split("|")
                    if len(parts) > canonical_idx and parts[canonical_idx] == "YES":
                        chosen = parts
                        break
            if chosen is None:
                chosen = transcripts[0].split("|")

            row = {"Chr": chrom.replace("chr", ""), "Start": pos, "Ref": ref, "Alt": alt, "Variants": vid}
            for f in VEP_KEEP_FIELDS:
                i = idx[f]
                row[f] = chosen[i] if i < len(chosen) else ""
            rows.append(row)

    return pl.DataFrame(rows)


def run_variant_annotation(
    vcf_file: str,
    gene: str,
    work_dir: str,
    annovar_pl: str,
    annovar_db: str,
    vep_cache_dir: str,
    vep_plugin_dir: str,
    cadd_snv_file: str,
    cadd_indel_file: str,
    franklin_file: str | None = None,
    genome_build: str = "hg38",
    vep_cache_version: str = "114",
    franklin_join: str = "left",
    out_dir: str | None = None,
    out_name: str | None = None,
) -> str:
    """
    Run ANNOVAR + VEP on vcf_file, reduce both to the fixed column sets,
    merge on the Variants (VCF ID) key, optionally merge in Franklin, and
    write one finalized annotation table.

    Parameters
    ----------
    work_dir : where ANNOVAR/VEP intermediate files are written
        (f"{gene}.annovar.*", f"{gene}.VEP.vcf").
    out_dir : directory the FINAL annotation table is written to. Defaults
        to work_dir if not given.
    out_name : base filename (no extension) for the final annotation table;
        written to f"{out_dir}/{out_name}.txt". Defaults to
        f"{gene}_annotation_final" if not given.
    franklin_join : "left" (default, keeps every ANNOVAR+VEP row and adds
        Franklin scores where available) or "inner" (keeps only variants
        present in the Franklin export).

    Returns
    -------
    str : path to the finalized annotation table.
    """
    Path(work_dir).mkdir(parents=True, exist_ok=True)
    out_dir = out_dir or work_dir
    out_name = out_name or f"{gene}_annotation_final"
    Path(out_dir).mkdir(parents=True, exist_ok=True)

    multianno_file = _run_annovar(vcf_file, gene, work_dir, annovar_pl, annovar_db, genome_build)
    annovar_df = _reduce_annovar(multianno_file)
    print(f"ANNOVAR reduced rows: {annovar_df.height}")

    vep_vcf = _run_vep(vcf_file, gene, work_dir, vep_cache_dir, vep_plugin_dir,
                        cadd_snv_file, cadd_indel_file, genome_build, vep_cache_version)
    vep_df = _reduce_vep(vep_vcf)
    print(f"VEP reduced rows: {vep_df.height}")

    merged = annovar_df.join(
        vep_df.drop(["Chr", "Start", "Ref", "Alt"]),
        on="Variants", how="left",
    )
    print(f"ANNOVAR+VEP merged rows: {merged.height}")

    if franklin_file is not None:
        franklin = pl.read_csv(franklin_file, separator="\t")
        franklin = franklin.rename({"Variants": "var"})
        franklin = franklin.with_columns(
            pl.col("var").str.replace("chr", "", literal=True).alias("tmp")
        )
        franklin = franklin.with_columns(
            pl.concat_str([pl.lit("chr"), pl.col("tmp")]).alias("Variants")
        ).drop(["var", "tmp"])
        merged = merged.join(franklin.select(["Variants", "Franklin"]), on="Variants", how=franklin_join)
        print(f"Franklin rows: {franklin.height}")
        print(f"Final merged rows (with Franklin): {merged.height}")

    out_file = str(Path(out_dir) / f"{out_name}.txt")
    merged.write_csv(out_file, separator="\t", include_header=True)
    return out_file


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _cli():
    p = argparse.ArgumentParser(description="Two-stage variant annotation pipeline.")
    p.add_argument("--input-file", required=True)
    p.add_argument("--file-type", required=True, choices=["plink1", "plink2", "vcf"])
    p.add_argument("--gene", required=True)
    p.add_argument("--work-dir", required=True, help="Where ANNOVAR/VEP intermediate files are written")
    p.add_argument("--chrom")
    p.add_argument("--start", type=int)
    p.add_argument("--end", type=int)
    p.add_argument("--maf-min", type=float, help="Minimum MAF -- keeps variants >= this (drops rare), computed on all samples before subsetting")
    p.add_argument("--maf-max", type=float, help="Maximum MAF -- keeps variants <= this (drops common), e.g. 0.01 to keep only rare variants; computed on all samples before subsetting")
    p.add_argument("--n-samples", type=int, default=10)
    p.add_argument("--sample-list")
    p.add_argument("--subset-out-dir", help="Directory for the Stage 1 subset VCF (default: work-dir)")
    p.add_argument("--subset-out-name", help="Base filename (no ext) for the Stage 1 subset VCF (default: <gene>.subset)")
    p.add_argument("--annovar-pl", required=True)
    p.add_argument("--annovar-db", required=True)
    p.add_argument("--vep-cache-dir", required=True)
    p.add_argument("--vep-plugin-dir", required=True)
    p.add_argument("--cadd-snv-file", required=True)
    p.add_argument("--cadd-indel-file", required=True)
    p.add_argument("--franklin-file")
    p.add_argument("--genome-build", default="hg38")
    p.add_argument("--vep-cache-version", default="114")
    p.add_argument("--out-dir", help="Directory for the FINAL annotation table (default: work-dir)")
    p.add_argument("--out-name", help="Base filename (no ext) for the FINAL annotation table (default: <gene>_annotation_final)")
    p.add_argument("--plink-path", default="plink2", help="Path to the plink2 binary (default: 'plink2', resolved via PATH)")
    p.add_argument("--bcftools-path", default="bcftools", help="Path to the bcftools binary. If bcftools misbehaves (e.g. 'invalid option' errors on modern flags like -O/-o), PATH is probably resolving to a stale/wrong binary -- pass the absolute path here instead, e.g. output of `which bcftools` run inside your activated env.")
    p.add_argument("--bgzip-path", default="bgzip", help="Path to the bgzip binary (default: 'bgzip', resolved via PATH)")
    p.add_argument("--tabix-path", default="tabix", help="Path to the tabix binary (default: 'tabix', resolved via PATH)")
    args = p.parse_args()

    Path(args.work_dir).mkdir(parents=True, exist_ok=True)
    subset_out_dir = args.subset_out_dir or args.work_dir
    subset_out_name = args.subset_out_name or f"{args.gene}.subset"
    vcf = prepare_variant_subset(
        input_file=args.input_file, file_type=args.file_type,
        out_dir=subset_out_dir, out_name=subset_out_name,
        chrom=args.chrom, start=args.start, end=args.end,
        maf_min=args.maf_min, maf_max=args.maf_max, n_samples=args.n_samples,
        sample_list=args.sample_list,
        plink_path=args.plink_path, bcftools_path=args.bcftools_path,
        bgzip_path=args.bgzip_path, tabix_path=args.tabix_path,
    )

    out_file = run_variant_annotation(
        vcf_file=vcf, gene=args.gene, work_dir=args.work_dir,
        annovar_pl=args.annovar_pl, annovar_db=args.annovar_db,
        vep_cache_dir=args.vep_cache_dir, vep_plugin_dir=args.vep_plugin_dir,
        cadd_snv_file=args.cadd_snv_file, cadd_indel_file=args.cadd_indel_file,
        franklin_file=args.franklin_file, genome_build=args.genome_build,
        vep_cache_version=args.vep_cache_version,
        out_dir=args.out_dir, out_name=args.out_name,
    )
    print(f"Final annotation table: {out_file}")


if __name__ == "__main__":
    _cli()
