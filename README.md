# Variant Interpretation Pipeline

End-to-end pipeline for annotating and classifying genomic variants (VCF) against
ACMG/AMP criteria, using InterVar, snpEff, ClinVar, gnomAD, dbNSFP, and RepeatMasker
annotations.

## Pipeline steps

Run in order (or via `scripts/run_pipeline_safe.sh`, which wraps all of them):

| Script | Purpose |
|---|---|
| `scripts/preflight_check.sh` | Environment / dependency checks before a run |
| `scripts/00_add_contigs.sh` | Add missing contig headers to input VCF |
| `scripts/01_validate.sh` | Validate VCF format |
| `scripts/02_normalize.sh` | Normalize variants (split multiallelics, left-align indels) |
| `scripts/03_annotate_snpeff.sh` | Functional annotation with snpEff |
| `scripts/04_clinvar_annotate.sh` | Annotate against ClinVar |
| `scripts/05_gnomad_annotate.py` | Add gnomAD population frequencies |
| `scripts/06_repeatmasker_annotate.py` | Flag repeat/low-complexity regions |
| `scripts/06c_dbnsfp_annotate.py` | Add dbNSFP in-silico predictions |
| `scripts/07_intervar_classify.sh` | ACMG/AMP classification via InterVar |
| `scripts/08_filter_prioritize.py` | Filter and prioritize classified variants |
| `scripts/09_report.py` | Generate final report |
| `diagnose_ba1.py` | Standalone diagnostic for BA1 (benign, stand-alone) calls |

## Usage

```bash
./scripts/preflight_check.sh
./scripts/run_pipeline_safe.sh
```

## Structure

```
.
├── scripts/                    # Pipeline steps (see table above)
├── diagnose_ba1.py             # BA1 diagnostic tool
├── test_data/
│   └── positive_controls/      # Known-pathogenic VCFs (HCM, LQTS, NF1, Coffin-Siris)
├── reports/                    # Per-gene research reports
└── docs/
    └── VariantInterpretationPipeline_Setup_Manual.md   # Full setup & installation guide
```

## Setup

See [`docs/VariantInterpretationPipeline_Setup_Manual.md`](docs/VariantInterpretationPipeline_Setup_Manual.md)
for full installation instructions (directory layout, required databases and tools,
and configuration).

## Test data

`test_data/positive_controls/` contains known-pathogenic variants for HCM, LQTS, NF1,
and Coffin-Siris syndrome, used as positive controls to validate pipeline output.

> Note: large reference databases (ClinVar, gnomAD, dbNSFP, InterVar db, etc.) are
> **not** checked into this repo — see the setup manual for how to fetch them.
