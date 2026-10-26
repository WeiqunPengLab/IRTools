# IRTools v2.1

## New features

- Add ONT direct RNA-seq intron retention and read-through quantification through `IRTools ont -q IRI` and `IRTools ont -q RT`.
- Use aligned-base coverage for the primary ONT IRI and RT metrics. ONT IRI retains retained/spliced event counts and event ratios as QC columns.
- Support ONT IRI and RT result files in `IRTools diff`, alongside existing short-read results. Automatically detect the input format or select it with `--input-format NGS` or `--input-format ONT`.
- Reject comparisons that mix NGS and ONT result formats, and require an explicit format when both formats exist for a sample.
- Add regression tests for differential input formats and ONT descriptive comparisons.

## Requirements and upgrade notes

- Python 3.11 or newer is required.
- Reinstall IRTools from the updated source folder, then run `IRTools --version` to confirm `IRTools 2.1`.
- For ONT coverage analysis, use the bundled non-mappability-adjusted `hg38` or `mm39` annotation, or an annotation generated without `-u/--map-file`.
- Existing `annotation`, `quant`, and `diff` commands remain available; `ont` provides the ONT workflow.
