# Variant Interpretation Pipeline — Setup & Installation Manual

**Environment:** WSL2 (Ubuntu, Noble 24.04), user `eman`
**Target build:** GRCh38
**Purpose:** End-to-end pipeline for annotating and classifying genomic variants (VCF) against ACMG/AMP criteria, using InterVar, ANNOVAR/VEP, and curated clinical databases (ClinVar, gnomAD, dbNSFP, ClinGen, HPO, OMIM).

---

## 1. Directory Structure

The pipeline lives at `~/VariantInterpretationPipeline` and uses a fixed top-level layout so every script can reference predictable relative paths.

```bash
mkdir VariantInterpretationPipeline
cd VariantInterpretationPipeline
mkdir input output databases tools scripts annotation reports logs tmp
```

| Directory | Purpose |
|---|---|
| `input/` | Raw and intermediate VCFs at each pipeline stage |
| `output/` | Final annotated/classified outputs |
| `databases/` | All reference and annotation databases (ClinVar, gnomAD, dbNSFP, etc.) |
| `tools/` | Miscellaneous helper binaries |
| `scripts/` | Pipeline step scripts (00–09) |
| `annotation/` | Annotation-stage working files |
| `reports/` | Final human-readable variant reports |
| `logs/` | Per-step execution logs |
| `tmp/` | Scratch space, safe to clear between runs |

Later, `logs/qc`, `output/qc`, and `scripts/qc` subfolders were added for quality-control artifacts, and `software/` was introduced separately to hold externally cloned tools (InterVar, ensembl-vep) distinct from custom `scripts/`.

Final top-level structure after full setup:

```
.
├── annotation
├── databases
│   ├── cadd
│   ├── clingen
│   ├── clinvar
│   ├── dbnsfp
│   ├── gnomad
│   ├── hpo
│   ├── omim
│   ├── revel
│   └── spliceai
├── input
├── logs
│   └── qc
├── output
│   └── qc
├── reports
├── scripts
│   └── qc
├── software
│   ├── InterVar
│   └── ensembl-vep
├── tmp
├── tools
└── workflow
    └── pipeline_workflow.md

26 directories, 1 file
```

---

## 2. Prerequisite Tools

Before installing any pipeline-specific software, confirm the core toolchain is present and working. These are the minimum binaries the pipeline depends on at various stages.

### 2.1 Version Checks

```bash
java -version
python3 --version
bcftools --version
samtools --version
tabix --version
bgzip --version
wget --version
curl --version
git --version
```

**Confirmed working versions in this setup:**

| Tool | Version |
|---|---|
| OpenJDK | 21.0.11 |
| Python | 3.12.3 |
| bcftools | 1.19 (htslib 1.19) |
| samtools | 1.19.2 (htslib 1.19) |
| tabix | 1.19 (htslib) |
| bgzip | 1.19 (htslib) |
| wget | GNU Wget 1.21.4 |
| curl | 8.5.0 |
| git | 2.43.0 |

> **Why this matters:** `bcftools`/`samtools`/`tabix`/`bgzip` all share the same `htslib` backend (1.19 here). Mismatched htslib versions across tools is a common source of cryptic downstream errors (malformed index files, silent parsing failures), so confirming they all report the same htslib version up front is a useful sanity check.

### 2.2 System Package Dependencies

A broader set of system and Perl packages is required for downstream components (ANNOVAR/InterVar's Perl usage, and later VEP). Install these proactively rather than one at a time as errors surface:

```bash
sudo apt update
sudo apt install -y \
  git \
  curl \
  build-essential \
  cpanminus \
  libdbi-perl \
  libdbd-mysql-perl \
  libjson-perl \
  libarchive-zip-perl \
  libwww-perl \
  libmodule-build-perl \
  libbio-perl-perl
```

| Package | Why it's needed |
|---|---|
| `build-essential` | C/C++ compiler toolchain — required to build htslib and any C extensions (e.g. Bio::DB::HTS for VEP) |
| `cpanminus` | Perl module installer, simplifies pulling in any additional CPAN dependencies |
| `libdbi-perl`, `libdbd-mysql-perl` | Perl database interface — required by VEP's Ensembl API modules |
| `libjson-perl` | JSON parsing in Perl, used by various annotation tools |
| `libarchive-zip-perl` | Archive extraction, used by VEP's installer for cache files |
| `libwww-perl` | Provides `LWP::Simple`, required for VEP's self-test and any HTTP-based Perl calls |
| `libmodule-build-perl` | Build system used to compile `Bio::DB::HTS` against htslib |
| `libbio-perl-perl` | BioPerl — general bioinformatics Perl module set, dependency of several annotation tools |

> **Note:** If any of these are already satisfied by a previous `apt install`, `apt` reports `is already the newest version` rather than erroring — safe to run this block idempotently at any point in setup.

---

## 3. Database Setup

All reference/annotation databases live under `databases/`, one subdirectory per source:

```bash
mkdir -p databases
cd databases
mkdir clinvar gnomad dbnsfp cadd revel spliceai clingen hpo omim
```

### 3.1 ClinVar

ClinVar provides curated clinical significance calls (Pathogenic / Benign / VUS / etc.) per variant — this is the primary source of real-world clinical evidence the pipeline's ACMG classification step (InterVar) relies on.

**Files present:**
```
databases/clinvar/
├── clinvar.vcf.gz
└── clinvar.vcf.gz.tbi
```

**Inspecting the header** (confirms format version, build, and available INFO fields):
```bash
bcftools view -h clinvar.vcf.gz | less
```

Key header facts confirmed:
- `##fileformat=VCFv4.1`
- `##fileDate=2026-07-06`
- `##source=ClinVar`
- `##reference=GRCh38`

Relevant INFO fields available for downstream filtering/annotation:
- `CLNSIG` — aggregate clinical significance (Pathogenic, Likely_benign, Uncertain_significance, etc.)
- `CLNDN` — ClinVar's preferred disease name
- `GENEINFO` — gene symbol and ID
- `CLNREVSTAT` — review status / confidence tier of the classification
- `ORIGIN` — allele origin (germline, somatic, de novo, etc.)
- `RS` — dbSNP rsID, where available

**Querying specific fields directly** (useful for spot-checks without invoking the full pipeline):
```bash
bcftools query \
  -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/CLNSIG\t%INFO/CLNDN\t%INFO/GENEINFO\n' \
  clinvar.vcf.gz | head
```

Example output confirms the file is queryable and fields resolve correctly:
```
1   66926   AG  A   Uncertain_significance   Retinitis_pigmentosa   OR4F5:79501
1   69134   AG  G   Likely_benign             not_specified          OR4F5:79501
1   69241   C   T   Uncertain_significance    not_specified          OR4F5:79501
```

> **Update note:** the pipeline was later re-pointed from an older `clinvar_20210501` ANNOVAR-format snapshot to a fresher `clinvar_20240917` build (see §6.2) — the VCF-format ClinVar shown here (`databases/clinvar/clinvar.vcf.gz`) is a separate copy from the ANNOVAR-generic-format file used by `table_annovar.pl`; keep both in sync when refreshing.

### 3.2 HPO (Human Phenotype Ontology)

Used for phenotype-driven variant prioritization — mapping gene/variant hits to standardized clinical phenotype terms.

```bash
cd databases/hpo
wget https://purl.obolibrary.org/obo/hp.obo
```

This follows a redirect chain (`purl.obolibrary.org` → `github.com/obophenotype/human-phenotype-ontology` → versioned GitHub release asset) before landing on the actual file. Confirmed successful download:

```
hp.obo   100%[===================>]  10.70M  1.40MB/s   in 6.2s
-rw-r--r-- 1 eman eman 11M Jun 23 16:43 hp.obo
```

### 3.3 dbNSFP — VEP Plugins

dbNSFP scores (SIFT, PolyPhen, CADD, REVEL, etc.) are consumed via VEP's plugin architecture rather than a flat annotation file. The plugin *code* (not the dbNSFP data itself) is pulled from Ensembl's plugin repository into VEP's plugin directory:

```bash
ls ~/.vep          # confirm VEP's default working directory exists
cd ~/.vep/Plugins
git clone https://github.com/Ensembl/VEP_plugins.git
```

This pulls in ~8,300 objects covering the full VEP plugin ecosystem (`AlphaMissense.pm`, `CADD.pm`, `REVEL.pm`, `SpliceAI.pm`, `Blosum62.pm`, and dozens more) — dbNSFP-specific scoring is enabled by pairing `dbNSFP.pm` from this set with an actual downloaded dbNSFP data file (tracked separately under `databases/dbnsfp/`).

### 3.4 CADD, REVEL, SpliceAI, ClinGen, OMIM

Directories created and reserved (`databases/cadd`, `databases/revel`, `databases/spliceai`, `databases/clingen`, `databases/omim`) — population/functional-impact and gene-dosage data sources to be populated as needed per-analysis. Not all are required for every pipeline run; they extend evidence coverage for specific ACMG criteria (e.g. ClinGen dosage sensitivity data strengthens PVS1 calls for loss-of-function variants).

---

## 4. VEP (Ensembl Variant Effect Predictor)

VEP serves as an alternate/complementary consequence annotator alongside SnpEff, and is the delivery mechanism for dbNSFP/CADD/REVEL/SpliceAI plugin scoring.

### 4.1 Clone and Install

```bash
cd ~/VariantInterpretationPipeline/software
git clone https://github.com/Ensembl/ensembl-vep.git
cd ensembl-vep
perl INSTALL.pl
```

### 4.2 Installation Dependencies Encountered

The installer surfaces missing dependencies incrementally rather than all at once. In order encountered:

| Missing dependency | Fix |
|---|---|
| `DBI` Perl module | `sudo apt install -y libdbi-perl libdbd-mysql-perl libarchive-zip-perl libjson-perl build-essential` |
| `bzlib.h`, `lzma.h` headers (needed to build htslib/Bio::DB::HTS) | `sudo apt install -y libbz2-dev liblzma-dev` |
| `Module::Build` Perl module (needed to build Bio::DB::HTS) | `sudo apt install -y libmodule-build-perl` |
| `LWP::Simple` (needed for VEP's internal self-test) | `sudo apt install -y libwww-perl` |

> **Pattern:** each `Can't locate X.pm in @INC` error maps to `sudo apt install -y lib<module-name-lowercase>-perl`. Installing the full dependency set from §2.2 up front avoids most of these round-trips on a fresh setup.

### 4.3 Installer Flow

The interactive installer performs, in order:
1. Clones/updates the Ensembl API (`ensembl`, `ensembl-variation`, `ensembl-funcgen`, `ensembl-compara`, `ensembl-io`)
2. Builds `htslib` from source (via `gcc`, using the bz2/lzma headers above)
3. Builds `Bio::DB::HTS` Perl bindings against that htslib
4. Prompts to create `~/.vep` as the default cache directory — **answer `y`**
5. Prompts for cache file selection — see §4.4
6. Prompts for FASTA — point at the existing pipeline reference (`reference/GRCh38.fa`) rather than re-downloading
7. Prompts for plugins — optional at install time; can be added later
8. Runs a self-test (requires `LWP::Simple`, per above)

**Re-running the installer** after a partial/failed attempt: it will ask *"Destination directory already exists. Do you want to overwrite it?"* — answer `y`. This is safe; it does not affect any separate/existing Ensembl API install, and will skip re-doing steps (like htslib compilation) that already completed successfully.

### 4.4 Cache Selection

The installer lists all species/build caches available from Ensembl (521 options at the time of this setup, ~315GB combined). For a human GRCh38 pipeline, the relevant entries are:

| # | File | Size | Notes |
|---|---|---|---|
| 516 | `homo_sapiens_vep_116_GRCh37` | 24 GB | Wrong build — skip |
| 517 | `homo_sapiens_merged_vep_116_GRCh38` | 30 GB | Merged Ensembl+RefSeq — larger, not needed for basic use |
| 518 | `homo_sapiens_refseq_vep_116_GRCh37` | 23 GB | Wrong build — skip |
| 519 | `homo_sapiens_refseq_vep_116_GRCh38` | 26 GB | RefSeq-only transcripts |
| 520 | `homo_sapiens_vep_116_GRCh37` | 24 GB | Wrong build — skip |
| **521** | **`homo_sapiens_vep_116_GRCh38`** | **27 GB** | **Selected — standard Ensembl-transcript GRCh38 cache** |

**Selection made:** `521` (standard GRCh38 cache). This is the correct choice for matching the rest of the pipeline, which is built entirely around GRCh38 coordinates. Expect this single download to take anywhere from 20 minutes to a few hours depending on connection speed.

### 4.5 Verifying Installation

```bash
perl vep --help | head -20
```

Confirmed versions across all Ensembl API components at time of setup:

```
ensembl            : 116.0d85231
ensembl-compara     : 116.8057d0e
ensembl-funcgen     : 116.90049ea
ensembl-io          : 116.6afb5dc
ensembl-variation   : 116.2fb834b
ensembl-vep         : 116.0
```

---

## 5. InterVar

InterVar implements the ACMG/AMP 2015 variant classification guidelines programmatically, taking ANNOVAR-annotated output and producing a final Pathogenic/Likely Pathogenic/VUS/Likely Benign/Benign call per variant with supporting evidence codes (PVS1, PS1–4, PM1–6, PP1–5, BA1, BS1–4, BP1–7).

### 5.1 Clone

```bash
cd ~/VariantInterpretationPipeline/software
git clone https://github.com/WGLab/InterVar.git
```

### 5.2 Contents

```
InterVar/
├── Intervar.py          # main executable
├── README.md
├── config.ini            # paths to ANNOVAR, databases, build version
├── docs/
├── example/
│   └── ex1.avinput       # sample input for smoke-testing
└── intervardb/            # bundled evidence-criteria databases
    ├── BP1.genes.hg19 / .hg38
    ├── BS2_hom_het.hg19 / .hg38
    ├── PM1_domains_with_benigns.hg19 / .hg38
    ├── PP2.genes.hg19 / .hg38
    ├── PS1.AA.change.patho.hg19 / .hg38
    ├── PS4.variants.hg19 / .hg38
    ├── PVS1.LOF.genes.hg19 / .hg38
    ├── ext.variants.hg19 / .hg38
    ├── knownGeneCanonical.txt.hg19 / .hg38
    ├── mim_adultonset.txt
    ├── mim_domin.txt
    ├── mim_orpha.txt
    ├── mim_pheno.txt
    ├── mim_recessive.txt
    └── orpha.txt / orpha.txt.utf8
```

> **Note the hg19/hg38 pairing throughout `intervardb/`.** Every evidence-criteria file ships in both builds — this pipeline consistently uses the `.hg38` variants throughout (`config.ini` and `Intervar.py` are configured accordingly; confirm `-buildver hg38` is set wherever InterVar is invoked).

`config.ini` is the central place where InterVar is told where ANNOVAR lives, which databases to use (including dbNSFP, when enabled — see below), and which genome build to assume. When updating database versions (e.g. refreshing ClinVar), both `Intervar.py` and `config.ini` need corresponding edits to keep protocol names in sync — a mismatch here (e.g. old `clinvar_20210501` reference lingering after a database refresh) is a common source of silent annotation gaps.

---

## 6. Operational Notes & Lessons From Setup

A few non-obvious things worth documenting for future reference or re-setup:

1. **dbNSFP is not automatically wired into InterVar.** It must be explicitly re-enabled in `Intervar.py`'s protocol list, with matching operation flags — a count mismatch between protocols and operations breaks annotation silently. When re-enabling, verify with:
   ```bash
   grep -n "protocol refGene" software/InterVar/Intervar.py
   ```

2. **ClinVar versioning matters.** Clinically significant genes (MYH7, MYBPC3, TNNT2, etc.) are curated frequently — a several-year-old ClinVar snapshot can miss real pathogenic variants added since. Refresh periodically via:
   ```bash
   perl software/annovar/annotate_variation.pl -buildver hg38 -downdb -webfrom annovar clinvar_YYYYMMDD software/annovar/humandb/
   ```
   then update all references to the old filename across `scripts/`, `Intervar.py`, and `config.ini`.

3. **ANNOVAR's generic ClinVar format has no gene-symbol column** — only `Chr, Start, End, Ref, Alt, CLNALLELEID, CLNDN, CLNDISDB, CLNREVSTAT, CLNSIG, ...`. Filtering by gene requires cross-referencing against `hg38_refGene.txt` coordinates rather than grepping the ClinVar file directly.

4. **Chromosome naming must match the reference FASTA exactly.** This pipeline's `reference/GRCh38.fa` uses `chr`-prefixed contig names (`chr1`, `chr7`, `chr11`, ...). Any hand-built or externally-sourced VCF using bare numeric chromosome names (ANNOVAR's internal convention) will fail silently or with segfaults during `bcftools norm` unless corrected first.

5. **Hand-built VCFs need complete headers.** At minimum: `##contig` lines (with correct lengths) for every chromosome used, and `##FORMAT`/`##INFO` definitions for every field referenced in the data rows. Missing either causes `bcftools norm` to fail with `Invalid BCF, ... not present in the header` even when writing plain VCF output.

6. **dbNSFP/consequence-predictor tools only score missense SNVs.** Indels, frameshifts, and intronic/regulatory variants will never receive PP3/BP4-type evidence from these tools regardless of how many are installed — that evidence gap is expected, not a bug. Loss-of-function variants should instead rely on **PVS1**, which depends on correct consequence annotation (frameshift/stopgain calls) and ClinGen gene-dosage/haploinsufficiency data, not missense predictors.

---

## 7. Verified End-to-End Checklist

- [x] Core toolchain installed and version-confirmed (Java, Python, bcftools, samtools, tabix, bgzip, wget, curl, git)
- [x] System/Perl dependencies installed
- [x] Directory scaffold created
- [x] ClinVar (VCF format) downloaded, indexed, and spot-checked
- [x] HPO ontology (`hp.obo`) downloaded
- [x] VEP cloned, dependencies resolved, installed with GRCh38 cache (#521)
- [x] VEP_plugins repository cloned into `~/.vep/Plugins`
- [x] InterVar cloned, `.hg38` evidence databases confirmed present
- [ ] dbNSFP data file downloaded and wired into `Intervar.py` protocol list *(tracked separately per-analysis)*
- [ ] CADD / REVEL / SpliceAI / ClinGen / OMIM data populated *(reserved directories, not yet populated)*

---

*Document generated from setup session logs, WSL2/Ubuntu Noble environment, GRCh38 build.*
