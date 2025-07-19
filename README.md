# SAS: an alignment-based metric for assembly errors

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.6+](https://img.shields.io/badge/python-3.6+-blue.svg)](https://www.python.org/downloads/)

## Overview

SAS (Sufficient Alignment Support) is a comprehensive toolkit for systematically evaluating genome assembly quality using multi-technology data integration. Unlike traditional methods that identify inconsistencies between assembly and aligned reads, SAS detects assembly errors by flagging erroneous windows without sufficient support from mapped reads. It combines short-read (Illumina/Element), HiFi (PacBio), and ONT (Oxford Nanopore) sequencing data to identify and classify both Base-level Errors (BE) and Structural Errors (SE) with high sensitivity and accuracy.

The SAS algorithm employs a dual-approach strategy with technology-specific optimizations:
- **BE (Base-level Errors)**: Detected in 50bp windows using high-precision reads from all platforms
- **SE (Structural Errors)**: Identified in 2kb windows using ultra-long ONT reads 



### Detailed Detection Methodology

#### 1. **BE (Base-level Errors) Detection**

**Window Strategy**: 
- **Window Size**: 50bp sliding windows across diploid assembly
- **Support Threshold**: >10% of overall mean sequencing depth from any platform
- **Data Sources**: Reads from Illumina/Element, HiFi, and ONT platforms

**Adaptive Thresholds for Special Regions**:
The support threshold is relaxed to >10% of local depth in challenging regions devoid of reads from high-accuracy platforms:

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
- **Terminal Regions** (within 10kb from chromosome ends): Relaxed criteria accounting for telomeric complexity, focus on event-based detection

**Key Features**:
- Multi-window allocation for spanning events
- False positive mitigation at chromosome termini  
- Context-specific thresholds for different genomic regions

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

## Input Requirements

### Required Files

1. **Short-read BAM**: High-quality Illumina or similar short-read data
2. **HiFi BAM**: PacBio HiFi long-read data
3. **ONT BAM**: Oxford Nanopore long-read data
4. **Reference FASTA**: Reference genome sequence

### File Format Requirements

- All BAM files must be sorted and indexed
- Reference FASTA should be indexed (`.fai` file will be generated automatically if not provided)
- Chromosome names must be consistent across all files

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
│   ├── structural_errors.anomalous_regions.tsv    # Detected structural anomalies
│   ├── structural_errors.window_stats.tsv         # Window-based statistics
│   └── structural_errors.anomaly_summary.tsv      # Summary of anomaly types
├── sas_pipeline_summary.txt         # Overall QV assessment and statistics
├── pipeline_config.json             # Configuration backup
└── sas_pipeline.log                 # Detailed execution log
```


## Configuration Options

### Basic Parameters

- `--threads`: Number of processing threads (default: 8)
- `--skip-filter`: Skip BAM filtering if files are pre-filtered
- `--no-be-analysis`: Skip BE (base-level errors) analysis
- `--no-se-analysis`: Skip SE (structural errors) analysis

### BE Analysis Parameters

- `--be-window-size`: Window size for BE analysis (default: 50bp)
- `--be-support-threshold`: Support threshold percentage (default: 10%)
- `--be-local-threshold`: Local depth threshold for special regions (default: 10%)

### SE Analysis Parameters

- `--se-window-size`: Window size for SE analysis (default: 2000bp)
- `--se-min-indel-length`: Minimum INDEL length for detection (default: 50bp)
- `--se-depth-range`: Acceptable depth range as percentage of mean (default: 20-200%)
- `--se-event-threshold`: Event support threshold for structural detection (default: 25%)
- `--se-terminal-region`: Distance from chromosome ends for special handling (default: 10kb)

### Special Region Handling

- `--telomere-region`: Distance from chromosome ends (default: 10kb)
- `--gc-threshold`: High GC content threshold (default: 90%)
- `--repeat-gc-threshold`: GC threshold for repeat regions (default: 60%)
- `--repeat-composition-threshold`: GA/CT composition threshold (default: 90%)


---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
