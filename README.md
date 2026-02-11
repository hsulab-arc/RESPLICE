# RESPLICE

[![DOI](https://zenodo.org/badge/746389362.svg)](https://doi.org/10.5281/zenodo.17488465)

Computational tools for analyzing RESPLICE RNA-seq data. This repository accompanies the paper:

> Chandrasekaran, S.\*, **Tau, C.\*** *et al.* Rewriting endogenous human transcripts with dual CRISPR-guided 3' trans-splicing. *Cell Systems* (2025). [doi:10.1016/j.cels.2025.101487](https://doi.org/10.1016/j.cels.2025.101487)

## Overview

RESPLICE (RNA-guided trans-splicing with CRISPR) is a system for rewriting endogenous human transcripts using dual CRISPR-guided 3' trans-splicing. This repository provides the computational pipeline for processing, aligning, and quantifying trans-splicing events from RNA-seq data, including measurement of on-target efficiency and off-target specificity.

## Repository Structure

```
├── build_indicies.sh        # Build STAR alignment indices for each target gene
├── staralign.sh             # End-to-end RNA-seq alignment and chimeric read detection
├── rnaseq_analysis.py       # Trans-splicing efficiency and specificity analysis
├── ratio_errorprop.py       # Fluorescence ratio analysis with error propagation
├── cleaned_filter_plot.ipynb # Data filtering and visualization notebook
├── NewReporterFull.gtf      # Reporter construct gene annotation
├── ReporterFullCis.fa       # Reporter construct sequence (cis configuration)
├── ReporterFullTrans.fa     # Reporter construct sequence (trans configuration)
├── ITGB1.gtf / ITGB1_trans.fa       # Integrin Beta 1 annotations and sequences
├── SMARCA4.gtf / SMARCA4_trans.fa   # SMARCA4 annotations and sequences
└── TFRC.gtf / TFRC_trans.fa         # Transferrin Receptor annotations and sequences
```

## Pipeline

### 1. Build alignment indices

`build_indicies.sh` downloads the hg38 human genome reference and constructs gene-specific STAR indices by concatenating the reference with custom GTF annotations and FASTA sequences for each target (Reporter, ITGB1, SMARCA4, TFRC, CD40L).

```bash
bash build_indicies.sh
```

### 2. Align and detect chimeric reads

`staralign.sh` runs the full alignment pipeline for a given sample:

1. Downloads FASTQ files from Google Cloud Storage
2. Trims adapters and filters by quality (Q10, min length 40 bp) with `cutadapt`
3. Aligns with STAR in chimeric read detection mode (two-pass mapping)
4. Extracts and indexes chimeric alignments with `samtools`
5. Runs downstream analysis via `rnaseq_analysis.py`

```bash
bash staralign.sh <sample_prefix> <gene_mode> <testing_flag>
```

**Arguments:**
- `sample_prefix` — sample identifier for file naming
- `gene_mode` — target gene (`ITGB1`, `SMARCA4`, `TFRC`, or default `Reporter`)
- `testing_flag` — set to `real` to export results to cloud storage

### 3. Quantify trans-splicing efficiency and specificity

`rnaseq_analysis.py` processes the aligned data to compute:

- **Efficiency** — ratio of on-target trans-spliced reads to total reads
- **Specificity** — ratio of on-target to off-target chimeric junctions
- Off-target mapping via BLAST against the hg38 genome

```bash
python rnaseq_analysis.py <gene_mode> <use_existing_flag>
```

### 4. Fluorescence ratio analysis

`ratio_errorprop.py` analyzes promoter screening flow cytometry data, computing green-to-blue fluorescence ratios with propagated uncertainties.

## Dependencies

- [STAR](https://github.com/alexdobin/STAR) — splice-aware RNA-seq aligner
- [samtools](http://www.htslib.org/) — BAM/SAM manipulation
- [cutadapt](https://cutadapt.readthedocs.io/) — adapter trimming
- [BLAST+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html) — sequence alignment (with hg38 database)
- Python 3 with `pandas`, `numpy`, `biopython`, `openpyxl`

## Supported Targets

| Target | GTF | FASTA | Description |
|--------|-----|-------|-------------|
| Reporter | `NewReporterFull.gtf` | `ReporterFullCis.fa`, `ReporterFullTrans.fa` | Fluorescent reporter construct |
| ITGB1 | `ITGB1.gtf` | `ITGB1_trans.fa` | Integrin Beta 1 |
| SMARCA4 | `SMARCA4.gtf` | `SMARCA4_trans.fa` | SWI/SNF chromatin remodeler |
| TFRC | `TFRC.gtf` | `TFRC_trans.fa` | Transferrin Receptor |

## Citation

If you use this code, please cite:

```bibtex
@article{chandrasekaran2025resplice,
  title={Rewriting endogenous human transcripts with dual CRISPR-guided 3' trans-splicing},
  author={Chandrasekaran, Sita and Tau, Cyrus and others},
  journal={Cell Systems},
  year={2025},
  doi={10.1016/j.cels.2025.101487}
}
```

## License

See the [LICENSE](LICENSE) file for details.
