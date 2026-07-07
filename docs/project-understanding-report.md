# Project understanding report: Cross-Species Contamination Detection Pipeline

_Last reviewed in this repository on 2026-07-07._

## How to read this report

This document separates three levels of confidence:

- **Repo-supported facts**: directly visible in code, tests, or documentation.
- **Reasonable inferences**: likely intent inferred from how files fit together.
- **Unknowns to confirm**: scientific or operational decisions that are not fully answerable from the repository alone.

The short version: this repository contains a Python package and Nextflow/HPC wrappers for finding non-human sequence signal in human whole-genome sequencing data. It extracts unmapped or optionally poorly mapped reads from human BAM/CRAM files, classifies those reads taxonomically with Kraken2, aggregates per-sample taxon counts into matrices, detects statistical outliers, and renders a static HTML cohort report.

---

## 1. Executive summary

### What problem is this project trying to solve?

Human sequencing projects, such as whole-genome sequencing (WGS), primarily aim to analyze DNA from a human participant. However, sequencing files can also contain reads from bacteria, viruses, fungi, laboratory reagents, sample swaps, environmental contamination, index hopping, library preparation artifacts, or other non-human sources. This repository tries to detect and summarize that non-human signal after a human alignment step.

The project is called **CSC**, short for **Cross-Species Contamination**. Its main workflow is:

```text
human BAM/CRAM
  → extract unmapped or low-quality-mapped reads
  → classify those reads against a Kraken2 taxonomic database
  → aggregate Kraken2 reports across samples
  → detect unusual sample/taxon outliers
  → generate HTML and TSV/JSON summaries
```

### What does the repo actually do?

The repo provides:

| Component | What it does | Main entry point |
|---|---|---|
| Extract | Uses `samtools` to pull unmapped reads, and optionally low-MAPQ reads, from BAM/CRAM into FASTQ. | `csc-extract` |
| Classify | Uses `Kraken2` to assign extracted reads to taxa. | `csc-classify` |
| Database helper | Downloads, validates, caches, and inspects Kraken2 databases. | `csc-db` |
| Aggregate | Parses Kraken2 reports into sample-by-taxon raw, CPM, per-rank, absolute-burden, and optional high-confidence matrices. | `csc-aggregate` |
| Detect | Applies statistical outlier detection to aggregated matrices. | `csc-detect` |
| Report | Builds a self-contained HTML report plus machine-readable sidecars. | `csc-report` |
| Nextflow | Runs extract → classify → aggregate → detect → summary at cohort scale. | `nextflow/main.nf` |
| 1000 Genomes scripts | SLURM/Apptainer scripts for a large remote-CRAM demo/workflow. | `tests/1000G/*.sh` |

### The most important practical point

This is not a raw FASTQ-to-variant-calls pipeline. It assumes you already have human-aligned BAM/CRAM files and uses reads that did not align well to the human reference as candidates for non-human classification.

---

## 2. Necessary background

### Human sequencing data

In WGS, a sequencer produces millions to billions of short DNA fragments called **reads**. For a human study, most reads should come from the human genome. The usual workflow aligns reads to a human **reference genome** such as GRCh38, producing an alignment file.

### Reads and pairs

A **read** is one sequenced DNA fragment. In **paired-end sequencing**, both ends of a DNA fragment are sequenced, producing read 1 (`R1`) and read 2 (`R2`). This repo preserves paired outputs where possible by writing separate R1/R2 FASTQ files plus singleton/other files.

### Reference genome and alignment

A **reference genome** is the coordinate system to which sequencing reads are aligned. An **alignment** describes where a read maps. Reads may be:

- confidently mapped to human;
- unmapped;
- poorly mapped, often indicated by low **MAPQ**;
- secondary/supplementary alignments; or
- duplicates.

This repo focuses on unmapped reads by default and can optionally include reads below a MAPQ threshold.

### FASTQ, SAM, BAM, and CRAM

| Format | Meaning | Role here |
|---|---|---|
| FASTQ | Text format containing reads and base-quality scores. Usually `.fastq` or `.fastq.gz`. | Output of extraction; input to Kraken2. |
| SAM | Text alignment format. | Conceptual format behind BAM/CRAM; not usually stored by this pipeline. |
| BAM | Binary compressed alignment file. | Input to `csc-extract`. |
| CRAM | Reference-compressed alignment file. Requires a reference FASTA for decoding. | Input to `csc-extract`; often used for large public datasets. |
| BAI/CRAI | Index files for BAM/CRAM. | Needed by `samtools idxstats`; important for efficient remote CRAM workflows. |
| VCF | Variant call format. | Not produced by this repo, but non-human reads can affect downstream variant calls if they remain in human analysis. |

### Non-human contamination

**Non-human contamination** means DNA signal in a human sequencing dataset that does not originate from the intended human genome. Possible sources include:

- microbes genuinely present in a biological sample;
- contamination during sample collection, DNA extraction, or library preparation;
- sequencing artifacts;
- reference/database artifacts;
- low-complexity reads that match many organisms poorly;
- human reads that fail to align because of reference differences or low quality.

The repo uses the term “kitome” for known reagent/environmental organisms that may appear broadly across many samples.

### Pre-alignment vs post-alignment filtering

- **Pre-alignment filtering** removes or classifies reads before human alignment, usually from raw FASTQ.
- **Post-alignment filtering** works after human alignment, often using unmapped reads from BAM/CRAM.

This repo is mainly a **post-alignment** workflow: it starts from BAM/CRAM and extracts reads that are unmapped or optionally poorly mapped.

### Kraken2, k-mers, and taxonomy

Kraken2 is a taxonomic classifier. It breaks reads into short substrings called **k-mers** and compares them to a database of k-mers associated with known organisms. Kraken2 assigns reads to taxonomic IDs using a lowest-common-ancestor style approach. This repo recommends a curated Kraken2 database called **PrackenDB**, described in the docs/code as one genome per species to reduce ambiguity from redundant genomes.

### Why non-human reads matter downstream

Non-human reads can cause problems if they are interpreted as human sequence:

- **SNV artifacts**: a microbial or contaminant read may align weakly to human and look like a single-nucleotide variant.
- **Structural variant artifacts**: split or discordant alignments can look like insertions, deletions, translocations, or other rearrangements.
- **Quality-control confusion**: broad low-level contamination may obscure whether a sample/library/batch is usable.
- **Biological misinterpretation**: a detected microbe might be real biology, lab contamination, or a database artifact.

This repo does not itself call variants; it helps quantify and report non-human signal that could affect downstream interpretation.

---

## 3. Repository map

| Path | Type | Why it matters |
|---|---|---|
| `README.md` | Main documentation | Provides quick start, module overview, CLI examples, database guidance, and Nextflow usage. |
| `pyproject.toml` | Python packaging config | Defines package name `csc`, version `0.2.0`, Python requirement, dependency on `pyyaml`, optional test dependency, and console scripts. |
| `Dockerfile` | Container build | Builds a Python 3.12 image with samtools, Kraken2, the CSC package, and tests. |
| `.github/workflows/ci.yml` | CI workflow | Runs tests across Python 3.9-3.12 and regenerates demo docs on main. |
| `.github/workflows/docker.yml` | CI/container workflow | Publishes or validates container image behavior. |
| `.github/workflows/remote-integration.yml` | Integration workflow | Exercises remote-data style workflows; includes cleanup steps. |
| `csc/default_config.yaml` | Default config | Central defaults for extract/classify/aggregate/detect/logging. |
| `csc/config.py` | Config loader | Loads defaults and recursively merges user YAML via `CSC_CONFIG` or CLI `--config`. |
| `csc/extract/` | Source code | Extraction module wrapping `samtools`. |
| `csc/classify/` | Source code | Kraken2 classification and database-management module. |
| `csc/aggregate/` | Source code | Kraken2 report parser, matrix builder, confidence-tier logic, taxonomy helpers. |
| `csc/detect/` | Source code | Outlier detection and detection report helpers. |
| `csc/report/` | Source code | Static HTML report generator and cohort-report visualization helpers. |
| `nextflow/` | Workflow code | DSL2 pipeline and modules for cohort-scale execution. |
| `docs/*.md` | User documentation | Module docs for extract, classify, aggregate, detect, report, configuration, pipeline, and testing. |
| `docs/index.html` | Generated demo report | Static example report generated from test/mock data. |
| `docs/per_sample_summary.tsv`, `docs/report_manifest.json` | Generated report sidecars | Example report outputs. |
| `tests/` | Tests and fixtures | Unit/integration tests, synthetic data generators, golden outputs. |
| `tests/1000G/` | SLURM demo scripts and manifest | Large-cohort 1000 Genomes remote CRAM extraction/classification/reporting workflow. |
| `research/current.md` | Notes | Present but empty or minimal in this checkout. |
| `LICENSE` | Legal | MIT license. |

### Documentation discrepancy to notice

`docs/README.md` still labels the `detect` module as a stub, while the main README and source code show `csc.detect` is implemented and tested. Treat the source code and current tests as more authoritative than that line.

---

## 4. Main workflows

### Workflow A: Python CLI, one sample at a time

**Purpose:** Analyze one BAM/CRAM or a small set of samples manually.

1. **Extract candidate non-human reads**

   ```bash
   csc-extract sample.bam -o extract/sample --threads 4
   ```

   This writes FASTQ files such as:

   - `sample.unmapped.R1.fastq.gz`
   - `sample.unmapped.R2.fastq.gz`
   - `sample.unmapped.singleton.fastq.gz`
   - `sample.unmapped.other.fastq.gz`
   - `sample.idxstats.tsv`
   - `sample.reads_summary.json`

2. **Classify extracted reads**

   ```bash
   csc-classify \
     extract/sample/sample.unmapped.R1.fastq.gz \
     extract/sample/sample.unmapped.R2.fastq.gz \
     --paired \
     --db /data/kraken2/prackendb \
     -o classify/sample \
     --sample-id sample \
     --threads 8
   ```

   This writes:

   - `sample.kraken2.report.txt`
   - `sample.kraken2.output.txt`

3. **Aggregate reports across samples**

   ```bash
   csc-aggregate classify/*/*.kraken2.report.txt \
     -o aggregate \
     --db-path /data/kraken2/prackendb \
     --kraken2-output classify/*/*.kraken2.output.txt
   ```

   This writes matrices such as `taxa_matrix_raw.tsv`, `taxa_matrix_cpm.tsv`, rank-filtered matrices, metadata, and optional confidence-tier matrices.

4. **Detect outliers**

   ```bash
   csc-detect aggregate/taxa_matrix_cpm.tsv -o detect --method all --rank-filter S G F
   ```

   This writes:

   - `flagged_samples.tsv`
   - `qc_summary.json`
   - `quarantine_list.txt`
   - optional per-rank and confidence-tier subdirectories

5. **Generate report**

   ```bash
   csc-report aggregate -o report/contamination_report.html --detect-dir detect
   ```

   This writes an HTML report plus sidecars such as `report_manifest.json` and per-sample/species TSV files.

### Workflow B: End-to-end Nextflow pipeline

**Purpose:** Automate extract → classify → aggregate → detect → summary for a cohort.

Input CSV example:

```csv
sample_id,file,reference
SAMPLE_001,/data/SAMPLE_001.bam,
SAMPLE_002,/data/SAMPLE_002.cram,/ref/GRCh38.fa
```

Command:

```bash
nextflow run nextflow/main.nf \
  --input_csv samples.csv \
  --kraken2_db /data/kraken2/prackendb \
  --db_path /data/kraken2/prackendb \
  --outdir results \
  -profile apptainer
```

Major automated steps:

1. Read sample sheet.
2. Run `csc-extract` per sample.
3. Run `csc-classify` per sample.
4. Collect all Kraken2 reports and outputs.
5. Run `csc-aggregate` once for the cohort.
6. Run `csc-detect` on the chosen matrix.
7. Produce summary outputs.

### Workflow C: 1000 Genomes HPC/SLURM workflow

**Purpose:** Run the project on 1000 Genomes high-coverage remote CRAMs using SLURM arrays and an Apptainer container.

The `tests/1000G` workflow is more operationally complex than the core package. It includes:

- a manifest of 3202 sample CRAM/CRAI URLs;
- SLURM array scripts for extraction and classification;
- submission wrappers;
- aggregation/detection/report jobs;
- remote CRAM index handling;
- Apptainer image pull/reuse.

A beginner should treat this as an HPC deployment example, not as the simplest starting point.

---

## 5. Input and output files

| Extension/file | Contains | Comes from | Used by | Produces/feeds |
|---|---|---|---|---|
| `.bam` | Binary human alignment records. | Human alignment workflow outside this repo. | `csc-extract`. | Unmapped/low-MAPQ FASTQ plus idxstats sidecars. |
| `.cram` | Reference-compressed human alignment records. | Human alignment workflow or public datasets. | `csc-extract`. | Same as BAM, but needs reference FASTA. |
| `.bai`, `.crai`, `.csi` | Alignment indexes. | Alignment/indexing tools. | `samtools idxstats`; remote workflows. | Fast read-count sidecars. |
| `.fa`, `.fasta` | Reference genome sequence. | Reference build such as GRCh38. | CRAM decoding. | Enables CRAM extraction. |
| `.fastq.gz`, `.fq.gz` | Sequencing reads and base qualities. | `csc-extract` output. | `csc-classify` / Kraken2. | Kraken2 reports and per-read outputs. |
| `.kraken2.report.txt` | Kraken2 taxonomic summary table. | `csc-classify`. | `csc-aggregate`. | Sample-by-taxon matrices. |
| `.kraken2.output.txt` | Per-read Kraken2 classifications. | `csc-classify`. | `csc-aggregate` confidence recomputation. | High-confidence tier matrices. |
| `taxonomy/nodes.dmp`, `taxonomy/names.dmp` | NCBI taxonomy files inside a Kraken2 DB. | Kraken2 database. | `csc-aggregate`, DB validation/helpers. | Domain annotation and confidence-tier recomputation. |
| `taxa_matrix_raw.tsv` | Raw direct-read counts by taxon and sample. | `csc-aggregate`. | `csc-detect`, `csc-report`. | Flags and report summaries. |
| `taxa_matrix_cpm.tsv` | Counts per million classified reads. | `csc-aggregate`. | `csc-detect`, `csc-report`. | Relative/compositional summaries. |
| `taxa_matrix_abs.tsv` | Counts per million total sequenced reads. | `csc-aggregate` when reads summaries are supplied. | `csc-detect`, `csc-report`. | Absolute burden and variant-impact sections. |
| `taxa_matrix_*_S.tsv`, `_G.tsv`, `_F.tsv` | Rank-filtered matrices. | `csc-aggregate`. | `csc-detect`, `csc-report`. | Species/genus/family-specific analyses. |
| `taxa_matrix_*_conf0p10.tsv` | High-confidence matrices. | `csc-aggregate` confidence tier. | `csc-detect`, `csc-report`. | Sensitive vs high-confidence comparison. |
| `aggregation_metadata.json` | Aggregation provenance and schema metadata. | `csc-aggregate`. | `csc-report`, auditing. | Reproducibility record. |
| `flagged_samples.tsv` | Outlier sample/taxon flags. | `csc-detect`. | User review, report. | Quarantine/review decisions. |
| `qc_summary.json` | Detection summary. | `csc-detect`. | Report and automation. | QC summary. |
| `quarantine_list.txt` | Unique flagged sample IDs. | `csc-detect`. | Human/automation review. | Candidate samples to investigate. |
| `.html` | Static interactive report. | `csc-report` or docs generation. | Humans. | Final review artifact. |

---

## 6. How the code works

### Package entry points

The package exposes console scripts in `pyproject.toml`:

- `csc-extract = csc.extract.cli:main`
- `csc-classify = csc.classify.cli:main`
- `csc-db = csc.classify.db_cli:main`
- `csc-aggregate = csc.aggregate.cli:main`
- `csc-detect = csc.detect.cli:main`
- `csc-report = csc.report.cli:main`

Each CLI parses arguments, loads default/user config where supported, and calls core Python functions.

### Extract module

Simplified logic:

1. Validate input path and extension (`.bam` or `.cram`).
2. Resolve the CRAM reference from `--reference` or `REF_PATH` if needed.
3. Build output FASTQ filenames.
4. Run `samtools idxstats` unless skipped, writing read-count sidecars.
5. Build a `samtools fastq` command.
6. If `--mapq` is supplied, pipe `samtools view -e 'flag.unmap || mapq < threshold'` into `samtools fastq`.
7. Exclude secondary, duplicate, and supplementary records to avoid double-counting.
8. Return paths and read counts.

Important technical detail: the code excludes SAM flags `0x100`, `0x400`, and `0x800`, meaning secondary alignments, PCR/optical duplicates, and supplementary alignments are not sent to Kraken2.

### Classify module

Simplified logic:

1. Check that Kraken2 is on `PATH`.
2. Validate the Kraken2 database directory contains `hash.k2d`, `opts.k2d`, and `taxo.k2d`.
3. Validate one FASTQ for single-end mode or exactly two FASTQs for paired mode.
4. Derive or accept a sample ID.
5. Run `kraken2 --db ... --output ... --report ... --confidence ... --threads ...`.
6. Return paths to the report and per-read output.

### Database helper module

The `csc-db` support code can:

- compute MD5/SHA-256 hashes;
- fetch a database from local path, URL, or S3 URI;
- unpack tar archives;
- validate required Kraken2 files;
- list and clean cached databases;
- estimate memory requirements;
- check for PrackenDB-compatible taxonomy files.

Important dependency assumption: S3 fetching requires the AWS CLI to be available.

### Aggregate module

Simplified logic:

1. Parse each Kraken2 report line into taxon records.
2. Use direct-read counts to build a sample-by-taxon matrix.
3. Filter low-count taxa using `min_reads`.
4. Write raw count matrices.
5. Normalize to CPM, meaning counts per million classified reads in that sample.
6. If reads-summary/idxstats sidecars are supplied, write absolute-burden matrices per million total sequenced reads.
7. Write per-rank matrices for species (`S`), genus (`G`), and family (`F`) by default.
8. If per-read Kraken2 outputs and database taxonomy are supplied, recompute confidence tiers such as `conf0p10` without rerunning Kraken2.
9. Write metadata JSON.

Important practical note: comments in the source say aggregation holds sample data in memory; `chunk_size` controls progress logging frequency, not memory use. For very large cohorts, memory behavior should be validated.

### Detect module

Simplified logic:

1. Load an aggregated matrix.
2. Optionally remove known kitome/environmental taxa.
3. Optionally subtract background per taxon.
4. For each taxon, compare each sample's value to the cohort distribution.
5. Flag outliers using one or more methods:
   - **MAD**: median absolute deviation, robust to outliers.
   - **IQR**: interquartile range fence.
   - **GMM**: two-component Gaussian mixture model for separated background/contaminated groups.
6. Write flagged sample/taxon rows, a QC summary JSON, and a quarantine list.
7. Discover and process sibling rank/tier/absolute matrices when present.

### Report module

Simplified logic:

1. Load aggregate matrices and metadata.
2. Load detect outputs if provided.
3. Compute species summaries, prevalence, burden, diversity, clustering, and variant-impact summaries.
4. Render SVG/HTML/JS/CSS without heavy plotting dependencies.
5. Write a self-contained HTML report plus sidecars such as `report_manifest.json`, `per_sample_summary.tsv`, and species/variant-impact TSVs.

The report code intentionally distinguishes **CPM** from **absolute burden** because those answer different biological questions.

---

## 7. How to run the project

### Required software

Minimum for core Python package:

```bash
python >=3.9
pip
samtools >=1.12
Kraken2
Kraken2 database
```

Optional/for scaling:

```bash
Nextflow >=22.10
Docker, Singularity, or Apptainer
SLURM for HPC scripts
AWS CLI for S3 database fetching
pysam and pytest for tests/synthetic data
MultiQC if you want to combine pipeline summaries with other QC outputs
```

### Install from source

```bash
# Create and activate your own environment first if desired.
pip install .

# For tests:
pip install pysam pytest-cov
pip install ".[test]"
```

### Verify command availability

```bash
csc-extract --version
csc-classify --version
samtools --version | head -1
kraken2 --version | head -1
```

### Recommended beginner test path

If you are new to this, do not start with all 1000 Genomes samples. Instead:

1. Run tests to confirm the package works locally.
2. Generate synthetic data if you want small BAM/CRAM fixtures.
3. Try one BAM file with `csc-extract`.
4. Classify a tiny extracted FASTQ against a small test Kraken2 database.
5. Aggregate two or three reports.
6. Generate a report.

### Nextflow run

```bash
nextflow run nextflow/main.nf \
  --input_csv samples.csv \
  --kraken2_db /data/kraken2/prackendb \
  --db_path /data/kraken2/prackendb \
  --outdir results \
  -profile docker
```

Notes:

- `--kraken2_db` is required for classification.
- `--db_path` is used by aggregation for taxonomy/domain and confidence-tier logic. It often should be the same directory as `--kraken2_db`.
- For CRAM input, include the `reference` column in the CSV.
- Use `-profile slurm`, `-profile apptainer`, or both depending on local Nextflow configuration and cluster rules.

### Common mistakes

- Running CRAM extraction without a reference FASTA.
- Supplying a Kraken2 database path that lacks `hash.k2d`, `opts.k2d`, or `taxo.k2d`.
- Confusing `--kraken2_db` and `--db_path` in Nextflow.
- Assuming CPM equals absolute sample burden. CPM is normalized among classified reads; absolute burden needs total sequenced-read denominators.
- Treating every microbial taxon as real biology without checking reagent/background patterns.
- Running high-confidence aggregation without per-read Kraken2 output files and taxonomy files.

---

## 8. External tools and dependencies

| Tool/package | What it does | Where it appears | What to know |
|---|---|---|---|
| Python | Main implementation language. | `csc/`, tests, packaging. | Requires Python >=3.9. |
| PyYAML | YAML parsing. | `pyproject.toml`, `csc/config.py`. | Only required runtime Python package in packaging metadata. |
| samtools | Reads BAM/CRAM, runs `idxstats`, extracts FASTQ. | `csc/extract/extract.py`, tests, Dockerfile. | CRAM often needs a reference FASTA; low-MAPQ mode uses samtools expression filtering. |
| Kraken2 | Taxonomic classifier. | `csc/classify/*`, Dockerfile, docs. | Requires large database; RAM can be high unless memory mapping is used. |
| PrackenDB | Recommended Kraken2 DB. | README, default config, DB helper. | Repo recommends one genome per species for species-level robustness. |
| Nextflow | Workflow engine. | `nextflow/*.nf`. | Useful for cohort-scale runs and resuming failed steps. |
| Docker | Container runtime/build system. | Dockerfile, Nextflow profile. | Container includes samtools, Kraken2, CSC. DB should be mounted. |
| Apptainer/Singularity | HPC container runtime. | `tests/1000G`, Nextflow profiles. | Common on clusters where Docker is not allowed. |
| SLURM | HPC scheduler. | `tests/1000G/*.sh`, Nextflow `slurm` profile. | Uses arrays and resource flags. |
| AWS CLI | S3 download helper. | `csc/classify/db.py`. | Needed only for `s3://` database fetches. |
| pytest | Test runner. | `tests/`, CI. | Optional test dependency. |
| pysam | Synthetic BAM/CRAM test generation. | Tests and Dockerfile. | Not listed as runtime dependency, installed separately in CI/Docker. |
| MultiQC | Aggregated QC reporting. | README/Nextflow docs. | Optional; consumes summary outputs. |

---

## 9. Current project status

### Appears complete or mature

- Python package structure and CLI entry points are implemented.
- Extract/classify/aggregate/detect/report modules all have tests.
- CI runs the test suite across multiple Python versions.
- Dockerfile builds an environment with samtools and Kraken2.
- Nextflow modules exist for the main workflow.
- Report generation is substantial, including cohort-focused HTML and sidecar files.
- Golden regression tests exist for aggregation/detection outputs.

### Appears experimental or operationally specific

- The 1000 Genomes scripts are substantial but tightly coupled to SLURM, Apptainer, EBI endpoints, and a specific output convention.
- High-confidence tier defaults are implemented, but scientific thresholds should still be validated for the actual dataset.
- GMM detection exists in code/tests, but statistical behavior should be reviewed on real cohorts.

### Documentation issues or discrepancies

- `docs/README.md` says detect is a stub, but code/tests show it is implemented.
- Some docs emphasize four modules, while the current workflow includes report and database helpers as first-class user tools.
- Nextflow `summary.nf` is described as MultiQC-compatible summary, while the richer HTML report appears to be produced by `csc-report` and/or the 1000G `generate_report.sh`; users should confirm exactly what the Nextflow summary process emits in their intended run.

### Search for TODO/FIXME/etc.

A repository-wide search found no major `TODO` or `FIXME` implementation markers in core code. It did find placeholders and comments such as temporary-file cleanup, a Nextflow placeholder-channel comment, and a stale `stub` mention in `csc/__init__.py` and `docs/README.md`.

---

## 10. Scientific and technical assumptions

| Assumption | Evidence/where | Why it matters |
|---|---|---|
| Input data are already aligned human BAM/CRAM files. | Extract module and README start from BAM/CRAM. | The workflow will miss non-human reads that aligned confidently to human unless low-MAPQ extraction captures them. |
| Unmapped reads are enriched for non-human sequence. | Extract-first design. | Some unmapped reads are human but hard to align; some contaminant reads may map to human-like regions. |
| CRAM decoding needs the correct reference FASTA. | Extract code requires reference for CRAM. | Wrong or missing reference can break extraction or produce incorrect reads. |
| Secondary/supplementary/duplicate alignments should be excluded. | Extract flag filter. | Prevents double-counting but may remove records someone else expected to inspect. |
| Kraken2 database quality strongly affects results. | PrackenDB recommendation. | Taxonomic calls are only as good as the database. |
| Species-level detection benefits from one-genome-per-species DB design. | PrackenDB docs/code comments. | Other databases may inflate or blur species-level signal. |
| Default Kraken2 confidence is sensitive (`0.0`). | Config and docs. | More sensitivity can mean more false positives. |
| High-confidence tier `0.1` is useful by default. | Config and docs. | Threshold should be validated for the actual read lengths, database, and study design. |
| Paired-end mode expects exactly two FASTQ files. | Classify validation. | Singleton files may require separate handling if scientifically important. |
| Sample IDs can be inferred from filenames unless supplied. | Extract/classify/aggregate helpers. | Naming collisions or inconsistent suffixes can silently confuse cohort matrices. |
| CPM denominator is classified reads, not total sequencing. | Aggregate/report docs/code. | CPM can exaggerate apparent burden in samples with very few classified reads. |
| Absolute burden requires reads-summary/idxstats sidecars. | Extract/aggregate design. | If sidecars are missing, variant-impact reporting is less complete. |
| Cohort outlier detection assumes enough comparable samples. | Detect design. | Small cohorts or strongly batch-structured data may produce unstable flags. |
| Kitome taxa must be supplied by user/context. | Config default empty list. | Known reagent taxa may otherwise be repeatedly flagged or misread as biology. |
| HPC scripts assume SLURM, Apptainer, network access, and specific endpoints. | `tests/1000G`. | May need local cluster adaptation. |

---

## 11. Risks, limitations, and things to verify

### Scientific risks

- **False positives:** Kraken2 may classify low-complexity, human-derived, or database-contaminant reads as microbial.
- **False negatives:** Real non-human reads that align confidently to human or are absent from the database may be missed.
- **Database dependence:** Changing databases can change taxon calls substantially.
- **Taxonomic ambiguity:** Closely related species can share k-mers; species-level assignments may be uncertain.
- **Kitome/background effects:** Reagent contaminants may appear in many samples and should not automatically be interpreted biologically.
- **Batch effects:** A sequencing lane, extraction plate, or reagent lot may create cohort patterns that statistical outlier methods do or do not capture.

### Technical risks

- **Memory scaling:** Aggregation stores data in memory; verify on very large cohorts.
- **Remote CRAM fragility:** EBI/FTP/network behavior can fail independently of code correctness.
- **Filename/sample-ID matching:** High-confidence per-read outputs are matched by filename conventions.
- **Reference mismatch:** CRAM extraction depends on reference compatibility.
- **Unpinned database version:** If PrackenDB or Kraken2 DB URLs change, results may not be reproducible unless database hashes/versions are recorded.
- **Stale docs:** Some documentation is behind source code.

### Validation questions before relying on results

- Does the pipeline remove or ignore real human reads that should remain available for downstream analysis?
- Are non-human reads being counted once per original fragment/read, not multiple times?
- Do high-confidence tier results agree with the sensitive tier for important taxa?
- Are known negative controls clean?
- Are positive spike-in or known-contamination controls detected?
- Are flagged samples correlated with batch, extraction date, plate, library prep, sequencing center, or ancestry/population labels?
- Are outputs reproducible if rerun with the same container and database?

---

## 12. Glossary

| Term | Plain-language definition |
|---|---|
| Absolute burden | Non-human read count normalized by total sequenced reads, often per million total reads. |
| Alignment | Placement of a sequencing read onto a reference genome. |
| BAM | Binary format for aligned sequencing reads. |
| BAI/CRAI/CSI | Index files that allow tools to quickly access parts of BAM/CRAM files. |
| Clade | A taxonomic group including an organism and its descendants. |
| Contamination | DNA signal not intended to be part of the primary sample/genome under study. |
| CPM | Counts per million; here usually per million classified reads in a sample. |
| CRAM | Compressed alignment format that usually needs a reference genome to decode. |
| FASTQ | Text format containing sequencing reads and quality scores. |
| GMM | Gaussian mixture model; a statistical model that tries to separate data into multiple normal-like groups. |
| HPC | High-performance computing cluster, often using a scheduler such as SLURM. |
| IQR | Interquartile range; the range between the 25th and 75th percentiles. |
| Kitome | Informal term for organisms associated with lab kits/reagents/environmental background. |
| Kraken2 | A k-mer-based taxonomic classifier for sequencing reads. |
| k-mer | A DNA substring of length k, used as a matching unit by classifiers. |
| MAPQ | Mapping quality; a score describing alignment confidence. |
| MAD | Median absolute deviation; robust spread statistic used for outlier detection. |
| Nextflow | Workflow engine for running pipelines reproducibly across local/HPC/cloud systems. |
| Paired-end | Sequencing approach where both ends of a DNA fragment are read. |
| PrackenDB | Recommended Kraken2 database in this repo; intended to contain one reference genome per species. |
| Rank | Taxonomic level such as species, genus, family. |
| Read | A sequenced DNA fragment. |
| Reference genome | Genome assembly used as the coordinate system for alignment. |
| SAM | Text format for alignments; BAM is its compressed binary version. |
| SLURM | Common HPC job scheduler. |
| Taxon/taxa | A named organism or taxonomic group. |
| VCF | Variant Call Format; not produced here, but relevant downstream. |

---

## 13. Suggested learning path

### Stage 1: Understand the biological goal

Learn:

- what aligned human BAM/CRAM files contain;
- why reads can be unmapped;
- what contamination means in WGS;
- why “microbial read” does not automatically mean infection or biology.

Inspect first:

- `README.md`
- `docs/pipeline.md`
- `docs/extract.md`
- `docs/classify.md`

### Stage 2: Learn the file formats

Practice commands:

```bash
samtools view -H sample.bam | head
samtools idxstats sample.bam | head
samtools flagstat sample.bam
zcat sample.unmapped.R1.fastq.gz | head
```

Learn what SAM flags, MAPQ, paired-end reads, and CRAM references mean.

### Stage 3: Run the smallest possible pipeline

Try:

```bash
csc-extract sample.bam -o work/extract --threads 2
csc-classify work/extract/sample.unmapped.R1.fastq.gz work/extract/sample.unmapped.R2.fastq.gz \
  --paired --db /path/to/db -o work/classify --sample-id sample --threads 2
csc-aggregate work/classify/*.kraken2.report.txt -o work/aggregate
csc-detect work/aggregate/taxa_matrix_cpm.tsv -o work/detect
csc-report work/aggregate -o work/report/contamination_report.html --detect-dir work/detect
```

### Stage 4: Read the code in pipeline order

Read these files in order:

1. `csc/extract/extract.py`
2. `csc/classify/classify.py`
3. `csc/classify/db.py`
4. `csc/aggregate/aggregate.py`
5. `csc/aggregate/confidence.py`
6. `csc/detect/detect.py`
7. `csc/report/report.py`
8. `nextflow/main.nf`
9. `nextflow/modules/*.nf`

### Stage 5: Learn workflow/HPC basics

Learn:

- what Nextflow channels/processes are;
- what SLURM `sbatch` and array jobs are;
- what Apptainer/Singularity containers are;
- how cluster resource requests work (`cpus`, `memory`, `time`);
- how to resume failed Nextflow runs.

### Stage 6: Ask scientific validation questions

Before interpreting a real cohort, discuss controls, database versions, thresholds, kitome taxa, and expected contamination sources with your colleague.

---

## 14. Questions for my colleague

### Project purpose and intended use

1. Is this repo meant for production QC, exploratory research, or both?
2. What specific biological question motivated the pipeline?
3. Is the primary goal to remove contaminated samples, annotate them, or discover real microbial signal?
4. Are the outputs intended to influence downstream variant calling or only QC reporting?

### Input data and references

5. Which human reference build should be assumed: GRCh37, GRCh38, T2T, or something else?
6. Are the BAM/CRAM files coordinate-sorted and indexed before this pipeline starts?
7. Should low-MAPQ mapped reads be included, or only unmapped reads?
8. For CRAMs, where should the exact reference FASTA come from?

### Kraken2/database choices

9. Which Kraken2 database version should be used for real analyses?
10. Is PrackenDB mandatory or just recommended?
11. Are database hashes/version dates recorded for reproducibility?
12. Are viral, fungal, bacterial, archaeal, protist, vector, and UniVec sequences all intended to be included?
13. How should human reads classified by Kraken2 be interpreted?

### Thresholds and detection

14. Why is the default high-confidence threshold `0.1` appropriate for this project’s datasets?
15. Should MAD, IQR, GMM, or `all` be the default detection method in production?
16. What `kitome_taxa` list should be used?
17. What constitutes a sample that should be quarantined?
18. Are there known positive and negative controls?

### Outputs and interpretation

19. Should users focus on raw counts, CPM, or absolute burden?
20. What report sections are considered decision-making sections versus exploratory visuals?
21. How should broad low-level taxa across many samples be handled?
22. How should rare high-burden taxa be validated?

### Operations and reproducibility

23. Which workflow is preferred: Python CLI, Nextflow, or 1000G SLURM scripts?
24. What cluster profile/resource settings are known to work?
25. Should runs be done in Docker, Apptainer, or a conda/pip environment?
26. Are output directories and filenames considered stable APIs for downstream tools?
27. What parts of the code are still under active development?
28. Are there any known bugs or datasets where the pipeline behaves poorly?

---

## Appendix: files inspected for this report

Important files inspected include:

- `README.md`
- `pyproject.toml`
- `Dockerfile`
- `.github/workflows/ci.yml`
- `.github/workflows/docker.yml`
- `.github/workflows/remote-integration.yml`
- `csc/default_config.yaml`
- `csc/config.py`
- `csc/extract/extract.py`
- `csc/extract/cli.py`
- `csc/classify/classify.py`
- `csc/classify/cli.py`
- `csc/classify/db.py`
- `csc/classify/db_cli.py`
- `csc/aggregate/aggregate.py`
- `csc/aggregate/confidence.py`
- `csc/aggregate/taxonomy.py`
- `csc/aggregate/cli.py`
- `csc/detect/detect.py`
- `csc/detect/cli.py`
- `csc/detect/report.py`
- `csc/report/report.py`
- `csc/report/cohort.py`
- `csc/report/cohort_report.py`
- `csc/report/interactive.py`
- `csc/report/svg.py`
- `csc/report/cli.py`
- `nextflow/main.nf`
- `nextflow/nextflow.config`
- `nextflow/modules/*.nf`
- `docs/*.md`
- `docs/index.html`
- `tests/*.py`
- `tests/golden/*`
- `tests/1000G/*.sh`
- `tests/1000G/README.md`
