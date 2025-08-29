# SAS: an alignment-based metric for assembly errors

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.6+](https://img.shields.io/badge/python-3.6+-blue.svg)](https://www.python.org/downloads/)

## Overview

SAS (Sufficient Alignment Support) is a comprehensive toolkit for systematically evaluating genome assembly quality using multi-technology data integration. Unlike traditional methods that identify inconsistencies between assembly and aligned reads, SAS detects assembly errors by flagging erroneous windows without sufficient support from mapped reads. It combines short-read (Illumina), HiFi (PacBio), and ONT (Oxford Nanopore) sequencing data to identify and classify both Base-level Errors (BE) and Structural Errors (SE) with high sensitivity and accuracy.

The SAS algorithm employs a dual-approach strategy with technology-specific optimizations:
- **BE (Base-level Errors)**: Detected in 50bp windows using high-precision reads from all platforms
- **SE (Structural Errors)**: Identified in 2kb windows using ultra-long ONT reads 



### Detailed Detection Methodology

#### 1. **BE (Base-level Errors) Detection**

**Window Strategy**: 
- **Window Size**: 50bp sliding windows across diploid assembly
- **Support Threshold**: >10% of overall mean sequencing depth from any platform
- **Data Sources**: High-precision reads from Illumina, HiFi, and ONT platforms

**Adaptive Thresholds for Special Regions**:
The support threshold is relaxed to >10% of local depth in challenging regions where high-accuracy platforms show bias:

1. **Telomeric Regions**: Within 10kb of chromosomal ends
2. **Degenerative Repeats**: GnA or CnT repeats where GA/CT% ≥90% and GC% ≥60%
3. **High GC Content**: Regions with GC% ≥90%

*Rationale*: In these regions, Element and HiFi platforms show strong bias, requiring reliance on ONT reads for support validation.

#### 2. **SE (Structural Errors) Detection**  

**Window Strategy**:
- **Window Size**: 2kb sliding windows across the genome
- **Data Source**: Ultra-long ONT reads exclusively  
- **Analysis Approach**: Multi-signal integration with position-aware criteria

**Multi-Signal Detection Framework**:

**Core Detection Modules**:
- **INDEL Events**: Large insertions/deletions (≥50bp) identified through CIGAR analysis
- **Clipping Analysis**: Soft/hard clipping events indicating potential misjoins (excluding chromosome termini)
- **Depth Anomalies**: Coverage deviations (<20% or >200% of mean) indicating collapses/duplications
- **Integration Logic**: Rule-based combination of signals with different criteria for chromosome ends vs. internal regions

**Position-Aware Analysis**:
- **Internal Regions**: Apply both depth-based and event-based criteria ( >50% of window reads)
- **Terminal Regions** (±10kb from chromosome ends): Relaxed criteria accounting for telomeric complexity, focus on event-based detection

**Key Features**:
- Multi-window allocation for spanning events
- False positive mitigation at chromosome termini  
- Context-specific thresholds for different genomic regions

##### Non-spanning Reads Detection

- For each 2kb step window, check whether any ONT read fully spans a 30kb window starting at that position (the last window always covers the chromosome end, with length 28-30kb). If no spanning read is found, the region is considered non-spanning.

- **Terminal correction:** After merging, if a region starts at the chromosome beginning, the last 28kb is removed; if a region ends at the chromosome end, the first 28kb is removed. This ensures only valid internal non-spanning regions are reported and avoids false positives at chromosome boundaries.

- **BED file:** All merged and terminal-corrected non-spanning regions are output as `se.non_spanning_regions.merged.bed` (0-based, BED format).


### Quality Value Calculation

The final assembly quality is quantified using the Phred-scaled Quality Value:

```
QV_SAS = -10 × log₁₀((SE_windows + BE_windows) / Assembly_Length)
```

**Where**:
- **SE_windows**: Number of 2kb windows flagged as structural errors
- **BE_windows**: Number of 50bp windows flagged as base-level errors  
- **Assembly_Length**: Total length of the reference assembly (in bases)




## Installation

**Required External Tools:**
- `samtools` (v1.10+)
- `bedtools` (v2.25+)
- `mosdepth` (v0.3.0+, for SE analysis)

### System Dependencies

**Ubuntu/Debian:**
```bash
# Update package lists
sudo apt-get update

# Install basic tools
sudo apt-get install -y samtools bedtools

# Install mosdepth (for SE analysis)
wget https://github.com/brentp/mosdepth/releases/download/v0.3.3/mosdepth
chmod +x mosdepth
sudo mv mosdepth /usr/local/bin/

# Verify installations
samtools --version
bedtools --version
mosdepth --version
```

**CentOS/RHEL/Rocky Linux:**
```bash
# Install EPEL repository
sudo yum install -y epel-release

# Install basic tools
sudo yum install -y samtools bedtools

# Install mosdepth manually
wget https://github.com/brentp/mosdepth/releases/download/v0.3.3/mosdepth
chmod +x mosdepth
sudo mv mosdepth /usr/local/bin/
```

**macOS with Homebrew:**
```bash
# Install basic tools
brew install samtools bedtools

# Install mosdepth
brew install mosdepth

# Verify installations
samtools --version
bedtools --version
mosdepth --version
```

**Alternative: Using Conda/Mamba (Recommended):**
```bash
# Create a new environment with all dependencies
conda create -n sas-pipeline python=3.8 samtools bedtools mosdepth -c bioconda
conda activate sas-pipeline

# Or if using mamba (faster)
mamba create -n sas-pipeline python=3.8 samtools bedtools mosdepth -c bioconda
mamba activate sas-pipeline
```

### Python Dependencies

**Basic installation:**
```bash
# Install core dependencies
pip install -r requirements.txt

# Or install specific packages
pip install "numpy>=1.16.0" "pandas>=1.0.0" "pysam>=0.15.0"
```


### SAS Pipeline Installation

```bash
# Clone the repository
git clone https://github.com/KANGYUlab/sas-pipeline.git
cd sas-pipeline

# Make scripts executable
chmod +x scripts/*.sh

# Test installation
python src/sas_pipeline.py --help


```

### Verification

Test your installation with a simple command:
```bash
# Check Python dependencies
python -c "import pysam, pandas, numpy; print('✓ Python dependencies OK')"

# Check external tools
samtools --version && echo "✓ samtools OK"
bedtools --version && echo "✓ bedtools OK"  
mosdepth --version && echo "✓ mosdepth OK"

# Check SAS Pipeline
python src/sas_pipeline.py --version
```

## Quick Start

### Basic Usage

```bash
# Run complete pipeline (BE + SE analysis)
python src/sas_pipeline.py \
  --short-read-bam illumina.bam \
  --hifi-bam hifi.bam \
  --ont-bam ont.bam \
  --reference-fa reference.fa \
  --output-dir results \
  --threads 64
```

### Advanced Usage

```bash
# Run only BE (base-level errors) analysis
python src/sas_pipeline.py \
  --short-read-bam illumina.bam \
  --hifi-bam hifi.bam \
  --ont-bam ont.bam \
  --reference-fa reference.fa \
  --output-dir results \
  --no-se-analysis \
  --threads 32

# Run only SE (structural errors) analysis with pre-filtered BAM
python src/sas_pipeline.py \
  --short-read-bam illumina.bam \
  --hifi-bam hifi.bam \
  --ont-bam ont_filtered.bam \
  --reference-fa reference.fa \
  --output-dir results \
  --no-be-analysis \
  --skip-filter \
  --threads 32

# Custom analysis with specific parameters
python src/sas_pipeline.py \
  --short-read-bam illumina.bam \
  --hifi-bam hifi.bam \
  --ont-bam ont.bam \
  --reference-fa reference.fa \
  --output-dir results \
  --threads 64 \
  --se-window-size 2000 \
  --se-min-indel-length 50
```

## Parameters Reference

### Required Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `--short-read-bam` | string | **Required.** Path to short-read BAM file (Illumina, Element, etc.). |
| `--hifi-bam` | string | **Required.** Path to HiFi BAM file (PacBio HiFi).  |
| `--ont-bam` | string | **Required.** Path to ONT BAM file (Oxford Nanopore).  |
| `--reference-fa` | string | **Required.** Path to reference genome FASTA file. Will be indexed automatically if `.fai` not provided. |
| `--output-dir` | string | **Required.** Output directory for all results. Will be created if it doesn't exist. |

### Optional Parameters

#### Basic Configuration

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--reference-fai` | string | auto-generated | Path to reference genome FASTA index file (.fai). If not provided, will be generated automatically from `--reference-fa`. |
| `--threads` | integer | 8 | Number of threads to use for parallel processing. Recommended: 32-128. |
| `--scripts-dir` | string | `{pipeline_base}/scripts` | Path to scripts directory containing analysis tools. |
| `--src-dir` | string | `{pipeline_base}/src` | Path to source directory containing pipeline modules. |

#### Pipeline Control

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--skip-filter` | flag | False | Skip BAM filtering step if input BAM files are already filtered. Use when BAMs are pre-processed. |
| `--no-be-analysis` | flag | False | Skip base error (BE) analysis completely. Only runs SE analysis. |
| `--no-se-analysis` | flag | False | Skip structural error (SE) analysis completely. Only runs BE analysis. |

#### SE Analysis Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--se-window-size` | integer | 2000 | Window size for SE analysis in base pairs. Standard is 2kb windows. |
| `--se-min-indel-length` | integer | 50 | Minimum INDEL length for SE analysis in base pairs. Events smaller than this are ignored. |

#### Bed Masking Options

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--bed-files` | string(s) | None | Paths to bed files for masking (e.g., `chr1.bed chr2.bed`). Multiple files are allowed. Regions in these files will be excluded from analysis. |

#### Information

| Parameter | Type | Description |
|-----------|------|-------------|
| `--version` | flag | Display version information and exit. |
| `--help` | flag | Display help information and exit. |




## Output Structure

```
results/
├── be_analysis/                      # BE (Base-level Errors) analysis results
│   ├── final_be_results.txt         # Final BE candidates
│   ├── smallerror.bed               # Initial error regions from short reads
│   ├── hifi_filtered_regions.txt    # HiFi cross-validation results
│   ├── ont_filtered_regions.txt     # ONT cross-validation results
│   └── merged_results.txt           # Multi-technology consensus results
├── se_analysis/                      # SE (Structural Errors) analysis results
│   ├── se.anomalous_regions.bed           # All anomalous regions (BED)
│   ├── se.anomalous_regions.tsv           # All anomalous regions (table)
│   ├── se.anomaly_summary.tsv             # Per-chromosome anomaly summary
│   ├── se.mosdepth.global.dist.txt        # Global depth distribution
│   ├── se.mosdepth.region.dist.txt        # Per-region depth distribution
│   ├── se.mosdepth.summary.txt            # Depth summary statistics
│   ├── se.non_spanning.bed                # All 2kb non-spanning windows
│   ├── se.non_spanning_regions.merged.bed # Merged non-spanning regions (BED)
│   ├── se.regions.bed.gz                  # All 2kb windows (gzipped BED)
│   ├── se.regions.bed.gz.csi              # Index for regions.bed.gz
│   ├── se.result.bed                      # Merged anomalous regions (BED)
│   └── se.window_stats.tsv                # Window-level statistics
├── sas_pipeline_summary.txt         # Overall QV assessment and statistics
├── pipeline_config.json             # Configuration backup
└── sas_pipeline.log                 # Detailed execution log
```

## Citation

If you use the SAS method in your research, please cite:

Yanan Chu, Zhuo Huang, Changjun Shao, Shuming Guo, Xinyao Yu, Jian Wang, Yabin Tian, Jing Chen, Ran Li, Yukun He, Jun Yu, Jie Huang, Zhancheng Gao, Yu Kang. Approaching an Error-Free Diploid Human Genome. bioRxiv. doi: https://doi.org/10.1101/2025.08.01.667781

