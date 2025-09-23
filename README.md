# OxyMetaG

Oxygen metabolism profiling from metagenomic data using Pfam domains. OxyMetaG predicts the percent relative abundance of aerobic bacteria in metagenomic reads based on the ratio of abundances of a set of 20 Pfams.

## Installation

### Using Conda (Recommended)

```bash
# Create environment with dependencies
conda env create -f environment.yml
conda activate oxymetag

# Install OxyMetaG
pip install oxymetag
```

### Using Pip

First install external dependencies:
- Kraken2
- DIAMOND
- KrakenTools
- R with mgcv and dplyr packages

Then install OxyMetaG:
```bash
pip install oxymetag
```

## Quick Start

### 1. Setup Kraken2 database
```bash
oxymetag setup
```

### 2. Extract bacterial reads
```bash
oxymetag extract -i sample1_R1.fastq.gz sample1_R2.fastq.gz -o BactReads -t 48
```

### 3. Profile samples
```bash
oxymetag profile -i BactReads -o diamond_output -t 8
```

### 4. Predict aerobe levels
```bash
# For modern DNA
oxymetag predict -i diamond_output -o results.tsv -m modern

# For ancient DNA
oxymetag predict -i diamond_output -o results.tsv -m ancient

# Custom cutoffs
oxymetag predict -i diamond_output -o results.tsv -m custom --idcut 50 --bitcut 30 --ecut 0.01
```

## Commands

### oxymetag setup
**Function:** Sets up the standard Kraken2 database for taxonomic classification.

**What it does:** Downloads and builds a comprehensive Kraken2 database containing bacterial, archaeal, and viral genomes (fungi excluded). This database is used by the `extract` command to identify and isolate bacterial sequences from metagenomic samples.

**Time:** 2-4 hours depending on internet speed and system performance.

**Output:** Creates a `kraken2_db/` directory with the complete database.

---

### oxymetag extract
**Function:** Extracts bacterial reads from metagenomic samples using taxonomic classification.

**What it does:** 
1. Runs Kraken2 to classify all reads in your metagenomic samples
2. Uses KrakenTools to extract only the reads classified as bacterial
3. Outputs cleaned bacterial-only FASTQ files for downstream analysis

**Input:** Raw metagenomic FASTQ files (paired-end or merged)
**Output:** Bacterial-only FASTQ files in `BactReads/` directory

**Arguments:**
- `-i, --input`: Input fastq.gz files (paired-end or merged)
- `-o, --output`: Output directory (default: BactReads)
- `-t, --threads`: Number of threads (default: 48)
- `--kraken-db`: Kraken2 database path (default: kraken2_db)

---

### oxymetag profile
**Function:** Profiles bacterial reads against oxygen metabolism protein domains.

**What it does:**
1. Takes bacterial-only reads from the `extract` step
2. Uses DIAMOND blastx to search against a curated database of 20 Pfam domains related to oxygen metabolism
3. Identifies protein-coding sequences and their functional annotations
4. Creates detailed hit tables for each sample

**Input:** Bacterial FASTQ files (uses R1 or merged reads only)
**Output:** DIAMOND alignment files (TSV format) in `diamond_output/` directory

**Arguments:**
- `-i, --input`: Input directory with bacterial reads (default: BactReads)
- `-o, --output`: Output directory (default: diamond_output)
- `-t, --threads`: Number of threads (default: 4)
- `--diamond-db`: Custom DIAMOND database path (optional)

---

### oxymetag predict
**Function:** Predicts aerobe abundance from protein domain profiles using machine learning.

**What it does:**
1. Processes DIAMOND output files with appropriate quality filters
2. Normalizes protein domain counts by gene length (reads per kilobase)
3. Calculates aerobic/anaerobic domain ratios for each sample  
4. Applies a trained GAM (Generalized Additive Model) to predict percentage of aerobes
5. Outputs final predictions with confidence metrics

**Input:** DIAMOND output directory from `profile` step
**Output:** Tab-separated file with aerobe predictions for each sample

**Arguments:**
- `-i, --input`: Directory with DIAMOND output (default: diamond_output)
- `-o, --output`: Output file (default: oxygen_predictions.tsv)
- `-t, --threads`: Number of threads (default: 4)
- `-m, --mode`: Filtering mode - 'modern', 'ancient', or 'custom' (default: modern)
- `--idcut`: Custom identity cutoff (for custom mode)
- `--bitcut`: Custom bitscore cutoff (for custom mode)  
- `--ecut`: Custom e-value cutoff (for custom mode)

## Filtering Modes

OxyMetaG includes three pre-configured filtering modes optimized for different types of DNA:

### Modern DNA (default)
**Best for:** Contemporary environmental samples, cultured isolates, high-quality metagenomes
- Identity ≥ 60%
- Bitscore ≥ 50
- E-value ≤ 0.001

### Ancient DNA  
**Best for:** Archaeological samples, paleogenomic data, degraded environmental DNA
- Identity ≥ 45% (accounts for DNA damage)
- Bitscore ≥ 25 (accommodates shorter fragments)
- E-value ≤ 0.1 (more permissive for low-quality data)

### Custom
**Best for:** Specialized applications or when you want to optimize parameters
- Specify your own `--idcut`, `--bitcut`, and `--ecut` values
- Useful for method development or unusual sample types

## Output

The final output (`oxygen_predictions.tsv`) contains:
- `SampleID`: Sample identifier extracted from filenames
- `ratio`: Aerobic/anaerobic domain ratio
- `aerobe_pfams`: Number of aerobic Pfam domains detected
- `anaerobe_pfams`: Number of anaerobic Pfam domains detected  
- `Per_aerobe`: **Predicted percentage of aerobic bacteria (0-100%)**

## Biological Interpretation

The `Per_aerobe` value represents the predicted percentage of aerobic bacteria in your sample based on functional gene content:

- **0-20%**: Predominantly anaerobic community (e.g., sediments, anoxic environments)
- **20-40%**: Mixed anaerobic community with some aerobic components
- **40-60%**: Balanced aerobic/anaerobic community
- **60-80%**: Predominantly aerobic community
- **80-100%**: Highly aerobic community (e.g., surface soils, oxic water)

## Citation

If you use OxyMetaG in your research, please cite:

```
Bueno de Mesquita, C.P. et al. (2025). OxyMetaG: Oxygen metabolism profiling 
from metagenomic data using protein functional domains. 
```

## License

GPL-3.0 License

## Support

For questions, bug reports, or feature requests, please open an issue on GitHub:
https://github.com/cliffbueno/oxymetag/issues
