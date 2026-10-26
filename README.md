# IRTools2 v2.1

IRTools2 ((I)ntron retention and (R)ead-through (Tools), version 2.1) analyzes intron retention and read-through in short-read RNA-seq and ONT direct RNA-seq libraries. This guide describes the v2.1 source code in this repository.

It reflects the code as it exists now, including:

- GTF and GFF3 support in `annotation`
- annotation-time mappability preprocessing
- `quant -g` file-or-genome-key behavior
- multifile BAM support in `quant`
- ONT direct RNA-seq intron retention and read-through counting from splice-aware BAMs
- `CIR_ratio.txt` summary export
- `RT` quantification for read-through index

## What's new in v2.1

- `ont -q IRI` and `ont -q RT` quantify ONT direct RNA-seq splice-aware BAM alignments using aligned-base coverage.
- `diff` accepts ONT IRI and RT results, with automatic NGS/ONT detection and an explicit `--input-format` option. A comparison must use one input format throughout.
- Regression tests cover NGS/ONT differential input detection, metric-column selection, mixed-format rejection, and ONT descriptive comparisons.

See [RELEASE_NOTES.md](RELEASE_NOTES.md) for the release notes.

## Overview

IRTools provides four subcommands:

```bash
IRTools annotation
IRTools quant
IRTools ont
IRTools diff
```

Main workflow:

1. Build an IRTools annotation file from a reference annotation.
2. Quantify intron retention or read-through from short-read BAM files with `IRI`, `IRC`, or `RT`, or quantify ONT direct RNA-seq intron retention or read-through with `ont`.
3. Optionally compare groups with `diff`.

## Installation

Recommended environment:

- Python `3.11+`
- install dependencies in a clean environment first

Example with `micromamba`:

```bash
micromamba create -n ir python=3.11 -y
micromamba activate ir
micromamba install -n ir -c conda-forge numpy scipy pandas networkx bx-python htseq pysam statsmodels -y
```

Install IRTools from source:

```bash
cd /path/to/IRTools-main
python -m pip install --no-build-isolation --no-deps .
```

If you are updating an existing install:

```bash
micromamba activate ir
python -m pip uninstall -y IRTools
cd /path/to/IRTools-main
python -m pip install --no-build-isolation --no-deps .
```

Verify the installed version after installation or upgrade:

```bash
IRTools --version
# Expected output: IRTools 2.1
```

## Command Summary

```bash
IRTools [-h] [-v] {annotation,quant,ont,diff} ...
```

## annotation

Build an IRTools annotation file from a source GTF or GFF3.

### Input

- `-g/--GTF-file`
  Input annotation file. Both GTF and GFF3 are supported.

- `-o/--annotation-file`
  Output IRTools annotation filename.

- `--outdir`
  Output directory.

### Optional mappability preprocessing

Mappability is now handled at annotation time, not quant time.

- `-u/--map-file`
  Mappability source. Supported forms:
  - genome shortcut: `hg19`, `hg38`, `mm9`, `mm10`, `mm11`
  - local `.bw` / `.bigWig`
  - local `.bedgraph` / `.bedgraph.gz`
  - local `.bed` / `.bed.gz`

- `--map-read-length`
  Supported values: `24`, `36`, `50`, `100`

- `--map-cutoff`
  Numeric mappability threshold. Default: `0.1`

- `--map-cache-dir`
  Optional writable directory for generated low-mappability BED cache files.
  By default, IRTools uses `IRTools/data/mappability_cache` inside the installed package.
  If that package cache is not writable when a new cache must be built, rerun with this option.

Rule used internally:

- score `< map-cutoff` -> low mappability -> excluded from effective intron length
- score `>= map-cutoff` -> retained

Species shortcuts use Hoffman/Bismap-hosted resources and build a cached blacklist BED under the installed IRTools package.

### Example

Without mappability:

```bash
IRTools annotation \
  -g annotation.gff3 \
  -o annotation.irtools.gff3 \
  --outdir .
```

With mappability:

```bash
IRTools annotation \
  -g annotation.gff3 \
  -o annotation.irtools.readlength50.cutoff0.1.gff3 \
  -u hg38 \
  --map-read-length 50 \
  --map-cutoff 0.1 \
  --outdir .
```

If the installed package cache is not writable, choose a project cache explicitly:

```bash
IRTools annotation \
  -g annotation.gff3 \
  -o annotation.irtools.readlength50.cutoff0.1.gff3 \
  -u hg38 \
  --map-read-length 50 \
  --map-cutoff 0.1 \
  --map-cache-dir /path/to/project/irtools_mappability_cache \
  --outdir .
```

### Notes

- Protein-level annotation is no longer required for IR annotation generation.
- The annotation builder now focuses on exon-derived RNA structure needed for IR analysis.
- If you use mappability, build the annotation once and reuse it for all matching BAMs.

## quant

Quantify intron retention or read-through from BAM files.

### Key changes in current version

- `quant` no longer uses `-u`
- `quant` no longer uses `--map-read-length`
- `quant` no longer uses `-e`
- `-g` is always required
- `-g` can be:
  - an actual IRTools annotation file path
  - a bundled genome key: `hg38` or `mm39`

### Bundled non-mappability annotations

This repository includes two prebuilt IRTools GFF3 annotations generated by `IRTools annotation` without `-u/--map-file`, so their constitutive intron effective lengths are not mappability-adjusted:

- `IRTools/data/gencode.v49.basic.annotation.level1.irtools.gff3.gz`
  Human GENCODE v49 basic level-1 annotation.
- `IRTools/data/gencode.vM38.basic.annotation.level1.irtools.gff3.gz`
  Mouse GENCODE vM38 basic level-1 annotation.

For `quant`, these can be used through bundled genome keys:

```bash
-g hg38
-g mm39
```

You can also pass the annotation file path directly if you prefer to be explicit.

### Required arguments

- `-q/--quant-type`
  `IRI`, `IRC`, or `RT`

- `-i/--alt-file`
  BAM file path, or a wildcard pattern such as `*.bam`

- `-p/--read-type`
  `single` or `paired`

- `-s/--library-type`
  `fr-unstranded`, `fr-firststrand`, or `fr-secondstrand`

- `-g/--annotation-file`
  IRTools annotation file path, or bundled genome key

- `-n/--name`
  Required by CLI parser
  - single-file mode: used as output prefix
  - multifile mode: ignored and replaced by each BAM stem

- `--outdir`
  Output directory

### IRC-specific argument

- `-m/--min_overlap`
  Default: `8`

### Single-file examples

IRI:

```bash
IRTools quant \
  -q IRI \
  -i sample.bam \
  -p single \
  -s fr-unstranded \
  -g annotation.irtools.gff3 \
  -n sample \
  --outdir sample_out
```

IRC:

```bash
IRTools quant \
  -q IRC \
  -i sample.bam \
  -p single \
  -s fr-unstranded \
  -g annotation.irtools.gff3 \
  -n sample \
  --outdir sample_out
```

RT:

```bash
IRTools quant \
  -q RT \
  -i sample.bam \
  -p single \
  -s fr-unstranded \
  -g annotation.irtools.gff3 \
  -n sample \
  --outdir sample_out
```

### Multifile example

Wildcard BAM input is supported directly by `quant`:

```bash
IRTools quant \
  -q IRI \
  -i "/path/to/*.bam" \
  -p single \
  -s fr-unstranded \
  -g annotation.irtools.gff3 \
  -n batch \
  --outdir /path/to/output
```

In multifile mode:

- each BAM is processed in one `IRTools quant` run
- sample name becomes the BAM filename stem
- output files are written into the same `--outdir`

Example outputs:

```text
sample1.quant.IRI.introns.txt
sample1.quant.IRI.genes.txt
sample2.quant.IRI.introns.txt
sample2.quant.IRI.genes.txt
CIR_ratio.txt
```

### quant outputs

IRI:

- `NAME.quant.IRI.introns.txt`
- `NAME.quant.IRI.genes.txt`

IRC:

- `NAME.quant.IRC.junctions.txt`
- `NAME.quant.IRC.introns.txt`
- `NAME.quant.IRC.genes.txt`

RT:

- `NAME.quant.RT.genes.txt`
- RT uses the gene's constitutive exonic regions as the gene-expression denominator, so its exon-side normalization is consistent with `IRI`.
- For each gene, RT is computed as:
  - `RT_RPKM = RT_read_count / (RT_length / 1000) / (eligible_read_count / 1000000)`
  - `gene_CER_RPKM = gene_CER_read_count / (gene_CER_length / 1000) / (eligible_read_count / 1000000)`
  - `gene_RT_index = RT_RPKM / gene_CER_RPKM`

### IRI genome-wide summary file

When running `IRI`, IRTools also writes:

```text
CIR_ratio.txt
```

This file is tab-separated and contains:

- `sample_name`
- `Total CIR`
- `Total CER`
- `CIR ratio`
- `CER ratio`
- `Eligible read count (total uniquely mapped reads/read-pairs)`

`CIR ratio` and `CER ratio` are the percentages printed on screen in the genome-wide summary.

## ont

Calculate intron retention or read-through from ONT direct RNA-seq BAM files generated by splice-aware aligners such as minimap2. Select the analysis with `-q IRI` or `-q RT`; the default is `IRI`.

ONT IRI uses accumulated aligned-base coverage so its primary `intron_IRI` and `gene_IRI` values follow the same CIR-versus-CER structure as short-read `quant -q IRI`:

- aligned bases are accumulated separately in each constitutive intronic region (CIR) and constitutive exonic region (CER)
- `CIR_RPKM = CIR_aligned_bases / (CIR_length / 1000) / (total_aligned_bases / 1000000)`
- `adjacent_CER_RPKM` is calculated from the CERs connected to the CIR by constitutive junctions
- `intron_IRI = CIR_RPKM / adjacent_CER_RPKM`
- `gene_IRI = gene_CIR_RPKM / gene_CER_RPKM`

ONT coverage uses the actual genomic CIR and CER lengths. Use the bundled non-mappability-adjusted annotation, or another annotation generated without `-u/--map-file`, for a consistent ONT-versus-NGS coverage comparison.

For easier comparison with short-read output, the reported `*_read_count` columns are mean-read-length equivalent counts. The actual accumulated evidence is also reported in `*_aligned_bases` columns.

The previous ONT event-based measurements remain as QC columns:

- a read is counted as `spliced` when a CIGAR `N` skipped interval matches the CIR boundaries within `--splice-tolerance`
- a read is counted as `retained` when its aligned blocks overlap at least `--min-overlap` intronic bases
- `intron_IR_event_ratio = retained_reads / (retained_reads + spliced_reads)`
- `gene_IR_event_ratio = gene_retained_reads / (gene_retained_reads + gene_spliced_reads)`

Example:

```bash
IRTools ont \
  -q IRI \
  -i sample.ont.minimap2.bam \
  -s fr-unstranded \
  -g hg38 \
  -n sample \
  --outdir sample_ont_out
```

Batch mode accepts wildcard BAM input:

```bash
IRTools ont \
  -q IRI \
  -i "/path/to/ont_bams/*.bam" \
  -s fr-unstranded \
  -g annotation.irtools.gff3 \
  -n batch \
  --outdir /path/to/ont_out
```

In multifile mode, each BAM is processed in one `IRTools ont` run, the sample name becomes the BAM filename stem, and all outputs are written into the same `--outdir`.

Main ONT-specific options:

- `-q/--quant-type`
  `IRI` or `RT`. Default: `IRI`
- `--min-overlap`
  For IRI event-count QC columns, minimum aligned bases inside a constitutive intronic region required to count a retained read. It does not affect coverage-based IRI. Default: `20`
- `--splice-tolerance`
  For IRI event-count QC columns, maximum bp distance between a CIGAR `N` skipped interval and a constitutive intron boundary to count a spliced read. It does not affect coverage-based IRI. Default: `10`
- `--min-mapq`
  Minimum mapping quality. Default: `0`
- `--include-secondary`
  Include secondary alignments. By default, secondary alignments are skipped.
- `--include-supplementary`
  Include supplementary alignments. By default, supplementary alignments are skipped.
- `--require-unique`
  Require `NH=1` when an `NH` tag exists.

ONT IRI outputs:

```text
NAME.ont.IRI.introns.txt
NAME.ont.IRI.genes.txt
ONT_IRI_ratio.txt
```

For ONT read-through quantification:

```bash
IRTools ont \
  -q RT \
  -i "/path/to/ont_bams/*.bam" \
  -s fr-unstranded \
  -g annotation.irtools.gff3 \
  -n batch \
  --outdir /path/to/ont_rt_out
```

ONT RT uses the same strand-aware read-through regions as short-read RT: the region begins at the furthest annotated transcript end and extends up to 10 kb, stopping before the next annotated gene. Aligned bases are accumulated in each gene's constitutive exonic regions and read-through region. Base totals are converted to mean-read-length equivalent counts for the familiar count columns, and RPKM values are normalized by total aligned bases. Therefore:

- `gene_CER_RPKM = gene_CER_aligned_bases / (gene_CER_length / 1000) / (total_aligned_bases / 1000000)`
- `RT_RPKM = RT_aligned_bases / (RT_length / 1000) / (total_aligned_bases / 1000000)`
- `gene_RT_index = RT_RPKM / gene_CER_RPKM`

ONT RT outputs:

```text
NAME.ont.RT.genes.txt
ONT_RT_ratio.txt
```

## diff

`diff` compares intron retention or read-through results between two sample groups.

Main inputs:

- `-q/--quant-type`
  `IRI`, `IRC`, or `RT`

- `--indir`
  directory containing quant result files

- `--input-format`
  `auto`, `NGS`, or `ONT`. The default `auto` detects `*.quant.*` or `*.ont.*` files. Specify the format when both versions exist for the same sample.

- `-s1/--s1-files`
  Required comma-separated sample or replicate names for sample group 1.

- `-s2/--s2-files`
  Required comma-separated sample or replicate names for sample group 2.

- `-t/--analysis-type`
  `P` for paired, `U` for unpaired. Paired analysis requires the same number of replicates in both groups.

- `-m/--method`
  `minus` or `log2fc`. Default `minus` reports `mean(S2) - mean(S1)`. `log2fc` reports `log2((mean(S2)+1e-6)/(mean(S1)+1e-6))`.

If either group has only one sample, `diff` writes a descriptive comparison with `PValue` and `FDR` set to `NA`.

SciPy t-test precision warnings are suppressed and the returned p-value is kept. If a p-value is truly not computable, `diff` keeps the feature and writes `PValue` and `FDR` as `NA`.

For paired analysis, if a feature row is missing from one replicate quant output, `diff` treats that missing value as `0` so the replicate pairs remain aligned.

- `-n/--name`
  output prefix

Supported result levels:

- `IRI`
  - intron level
  - gene level
  - accepts short-read `NAME.quant.IRI.*.txt` and ONT `NAME.ont.IRI.*.txt`
- `IRC`
  - junction level
  - intron level
  - gene level
- `RT`
  - gene level only
  - accepts short-read `NAME.quant.RT.genes.txt` and ONT `NAME.ont.RT.genes.txt`

RT example:

```bash
IRTools diff \
  -q RT \
  --indir /path/to/rt_out \
  -s1 sample1_rep1,sample1_rep2 \
  -s2 sample2_rep1,sample2_rep2 \
  -t U \
  -n sample1_vs_sample2 \
  --outdir /path/to/diff_out
```

This reads files such as:

- `sample1_rep1.quant.RT.genes.txt`
- `sample2_rep2.quant.RT.genes.txt`

and writes:

- `sample1_vs_sample2.diff.RT.genes.txt`

ONT IRI and RT use the same command structure. When `--indir` contains ONT results, `diff` automatically reads files such as `sample1_rep1.ont.IRI.introns.txt`, `sample1_rep1.ont.IRI.genes.txt`, or `sample1_rep1.ont.RT.genes.txt`. For example:

```bash
IRTools diff \
  -q IRI \
  --input-format ONT \
  --indir /path/to/ont_iri_out \
  -s1 control_rep1,control_rep2 \
  -s2 treated_rep1,treated_rep2 \
  -t U \
  -n control_vs_treated \
  --outdir /path/to/diff_out
```

All samples in one comparison must use the same format. ONT IRI differential analysis uses the coverage-based `intron_IRI` and `gene_IRI` columns; ONT RT uses `gene_RT_index`. The event-ratio QC columns are not used by `diff`.

## Practical notes

- BAM files should be coordinate-sorted.
- BAM index files should be present when your normal workflow expects them.
- Chromosome naming mismatches such as `chr1` vs `1` are handled more robustly than in the older version.
- `quant` now uses `pysam` for BAM iteration, while `HTSeq` is still used for annotation parsing and interval data structures.
- If using paired-end mode, validate one known sample after upgrading, because BAM-pair handling is one of the more sensitive code paths.
- `RT` uses `transcript_region` to find the furthest annotated TES per gene, then defines a downstream read-through region up to 10 kb or until the next annotated gene on the same chromosome.

## Recommended workflow

1. Prepare a clean environment.
2. Build one IRTools annotation for your reference annotation.
3. If needed, include mappability during `annotation`.
4. Reuse that IRTools annotation for all `quant` runs.
5. Use wildcard BAM input if multiple BAMs share the same output folder.

## Example end-to-end

Build annotation:

```bash
IRTools annotation \
  -g /path/to/annotation.gff3 \
  -o annotation.irtools.readlength50.cutoff0.1.gff3 \
  -u hg38 \
  --map-read-length 50 \
  --map-cutoff 0.1 \
  --outdir /path/to/annotation_out
```

Quantify multiple BAMs:

```bash
IRTools quant \
  -q IRI \
  -i "/path/to/*.bam" \
  -p single \
  -s fr-unstranded \
  -g /path/to/annotation_out/annotation.irtools.readlength50.cutoff0.1.gff3 \
  -n batch \
  --outdir /path/to/iri_out
```
