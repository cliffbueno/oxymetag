# OxyMetaG

Oxygen metabolism profiling from metagenomic data using Pfam domains. OxyMetaG predicts the percent relative abundance of aerobic bacteria in metagenomic reads based on the ratio of abundances of a set of 20 Pfams. It is recommended to use a HPC cluster or server rather than laptop to run OxyMetaG due to memory requirements, particularly for the step of extracting bacterial reads. If you already have bacterial reads, the "profile" and "predict" functions will run quickly on a laptop.

If you are working with modern metagenomes, we recommend first quality filtering the raw reads with your method of choice and standard practices, and then extracting bacterial reads with Kraken2 and KrakenTools, which is performed with the OxyMetaG extract function. For profiling modern metagenomes, use DIAMOND blastx with the default `-m diamond` mode for the profile step. You can also use `-m custom` in the predict step to test different hit cutoffs.

If you are working with ancient metagenomes, we recommend first quality filtering the raw reads with your method of choice and standard practices, and then extracting bacterial reads with a workflow optimized for ancient DNA, such as the read mapping approach employed by De Sanctis et al. (2025). For profiling ancient metagenomes, use MMseqs2 with `-m mmseqs2` for the profile step and `-m ancient` for the predict step. The ancient mode uses parameters optimized for ancient DNA along with 97 decoy Pfams to reduce instances of false positives. We are still working on optimizing the methods for ancient DNA, which will be released as v2.0.0. 

## Installation

First clone the repository.

```bash
git clone https://github.com/cliffbueno/oxymetag.git
cd oxymetag
```

### Using Conda (Recommended)

```bash
# Create environment with dependencies
conda env create -f environment.yml
conda activate oxymetag

# Install OxyMetaG
pip install -e .

# Index the MMseqs2 database (one-time setup, ~5-10 minutes)
mmseqs createindex oxymetag/data/oxymetag_pfams_n117_db tmp
```

**Note:** The MMseqs2 database indexing is optional but highly recommended for faster searches.

### Using Pip

First install external dependencies:
- Kraken2
- DIAMOND
- MMseqs2
- KrakenTools
- R with mgcv, dplyr, tidyr, and rlang packages

Then install OxyMetaG:
```bash
pip install oxymetag
```

## Quick Start

### Modern DNA workflow

```bash
# 1. Setup Kraken2 database (one-time)
oxymetag setup

# 2. Extract bacterial reads
oxymetag extract -i sample1_R1.fastq.gz -o BactReads -t 48

# 3. Profile samples with DIAMOND
oxymetag profile -i BactReads -o diamond_output -m diamond -t 8

# 4. Predict aerobe levels
oxymetag predict -i diamond_output -o per_aerobe_predictions.tsv -m modern
```

### Ancient DNA workflow

```bash
# 1. Extract bacterial reads with an ancient DNA-optimized workflow (not currently provided by OxyMetaG)

# 2. Profile samples with MMseqs2
oxymetag profile -i BactReads -o mmseqs_output -m mmseqs2 -t 8

# 3. Predict aerobe levels with ancient mode
oxymetag predict -i mmseqs_output -o per_aerobe_predictions.tsv -m ancient
```

### Custom parameters

```bash
oxymetag predict -i diamond_output -o per_aerobe_predictions.tsv -m custom --idcut 50 --bitcut 30 --ecut 0.01
```

## Commands

### oxymetag setup
**Function:** Sets up the standard Kraken2 database for taxonomic classification.

**What it does:** Downloads and builds the standard Kraken2 database containing bacterial, archaeal, and viral genomes. This database is used by the `extract` command to identify bacterial sequences from metagenomic samples.

**Time:** Depends on internet speed and system performance, but will likely take several hours (4-20 hours typical). 

**Output:** Creates a `kraken2_db/` directory with the standard database (~90 GB).

Make sure you run oxymetag setup from the directory where you want the database to live, or plan to always specify the --kraken-db path when running extract. The database is quite large, so choose a location with sufficient storage.

---

### oxymetag extract
**Function:** Extracts bacterial reads from metagenomic samples using taxonomic classification.

**What it does:** 
1. Runs Kraken2 to classify all reads in your metagenomic samples
2. Uses KrakenTools to extract only the reads classified as bacterial
3. Outputs cleaned bacterial-only FASTQ files for downstream analysis

**Input:** Quality filtered metagenomic read FASTQ files (paired-end or merged)  
**Output:** Bacterial-only FASTQ files in `BactReads/` directory

**Arguments:**
- `-i, --input`: Input fastq.gz files (paired-end or merged)
- `-o, --output`: Output directory (default: BactReads)
- `-t, --threads`: Number of threads (default: 48)
- `--kraken-db`: Kraken2 database path (default: kraken2_db)

---

### oxymetag profile
**Function:** Profiles bacterial reads against oxygen metabolism protein domains using DIAMOND or MMseqs2.

**What it does:**
1. Takes bacterial-only reads from the `extract` step
2. Uses DIAMOND blastx (for modern DNA) or MMseqs2 (for ancient DNA) to search against a curated database of Pfam domains related to oxygen metabolism
   - DIAMOND mode: 20 target Pfams
   - MMseqs2 mode: 20 target Pfams + 97 decoy Pfams (117 total) to reduce false positives
3. Identifies protein-coding sequences and their functional annotations
4. Creates detailed hit tables for each sample

**Input:** Bacterial FASTQ files (uses R1 or merged reads only)  
**Output:** Alignment files (TSV format) in `diamond_output/` or `mmseqs_output/` directory

**Arguments:**
- `-i, --input`: Input directory with bacterial reads (default: BactReads)
- `-o, --output`: Output directory (default: diamond_output or mmseqs_output depending on method)
- `-t, --threads`: Number of threads (default: 4)
- `-m, --method`: Profiling method - 'diamond' or 'mmseqs2' (default: diamond)
- `--diamond-db`: Custom DIAMOND database path (optional, for diamond method)
- `--mmseqs-db`: Custom MMseqs2 database path (optional, for mmseqs2 method)

---

### oxymetag predict
**Function:** Predicts aerobe abundance from protein domain profiles using machine learning.

**What it does:**
1. Processes DIAMOND or MMseqs2 output files with appropriate quality filters
2. Selects the top hit per read based on bitscore
3. For MMseqs2 (ancient mode): filters out decoy Pfams after selecting top hits
4. Normalizes protein domain counts by gene length (reads per kilobase)
5. Calculates aerobic/anaerobic domain ratios for each sample  
6. Applies a trained GAM (Generalized Additive Model) to predict percentage of aerobes
7. Outputs a table with the sampleID, # Pfams detected, and predicted % aerobic bacteria

**Input:** DIAMOND or MMseqs2 output directory from `profile` step  
**Output:** Tab-separated file with aerobe predictions for each sample

**Mode determines input type:**
- `-m modern`: Uses DIAMOND output (default input: diamond_output/)
- `-m ancient`: Uses MMseqs2 output (default input: mmseqs_output/)
- `-m custom`: Auto-detects DIAMOND or MMseqs2 files in input directory

**Arguments:**
- `-i, --input`: Directory with profiling output (default: diamond_output for modern, mmseqs_output for ancient)
- `-o, --output`: Output file (default: per_aerobe_predictions.tsv)
- `-t, --threads`: Number of threads (default: 4)
- `-m, --mode`: Filtering mode - 'modern', 'ancient', or 'custom' (default: modern)
- `--idcut`: Custom identity cutoff (for custom mode)
- `--bitcut`: Custom bitscore cutoff (for custom mode)  
- `--ecut`: Custom e-value cutoff (for custom mode)

## Filtering Modes

OxyMetaG includes three pre-configured filtering modes optimized for different types of DNA. In any case, it is always recommended to try several different parameters (using -m custom) to check how sensitive the results are to the cutoffs.

### Modern DNA (default)
**Best for:** Modern environmental metagenomes  
**Method:** DIAMOND blastx  
**Filters:**
- Identity ≥ 60%
- Bitscore ≥ 50
- E-value ≤ 0.001

**Usage:**
```bash
oxymetag profile -m diamond
oxymetag predict -m modern
```

### Ancient DNA
**Best for:** Archaeological samples, paleogenomic data, degraded environmental DNA  
**Method:** MMseqs2 with decoy Pfams  
**Filters:**
- Identity ≥ 86%
- Bitscore ≥ 50
- E-value ≤ 0.001

**Note:** The ancient mode uses stricter identity cutoffs but employs 97 decoy Pfams to reduce false positives from damaged DNA. Reads matching decoys better than target Pfams are filtered out.

**Usage:**
```bash
oxymetag profile -m mmseqs2
oxymetag predict -m ancient
```

### Custom
**Best for:** Specialized applications or when you want to optimize parameters
- Specify your own `--idcut`, `--bitcut`, and `--ecut` values
- Auto-detects whether input is from DIAMOND or MMseqs2
- Useful for method development or unusual sample types

**Usage:**
```bash
oxymetag predict -m custom --idcut 50 --bitcut 30 --ecut 0.01
```

## Output

The final output (`per_aerobe_predictions.tsv`) contains:
- `SampleID`: Sample identifier extracted from filenames
- `ratio`: Aerobic/anaerobic domain ratio
- `aerobe_pfams`: Number of aerobic Pfam domains detected (from 20 target Pfams)
- `anaerobe_pfams`: Number of anaerobic Pfam domains detected (from 20 target Pfams)
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
Bueno de Mesquita, C.P., Stallard-Olivera, E., Fierer, N. (2025). 
Quantifying the oxygen preferences of bacterial communities using a 
metagenome-based approach.
```

### Additional citations

If you use the **extract** function, also cite Kraken2 and KrakenTools:
```
Lu, J., Rincon, N., Wood, D.E. et al. Metagenome analysis using the Kraken 
software suite. Nat Protoc 17, 2815–2839 (2022). 
https://doi.org/10.1038/s41596-022-00738-y
```

If you use the **profile** function with DIAMOND (`-m diamond`), also cite:
```
Buchfink, B., Xie, C. & Huson, D. Fast and sensitive protein alignment using 
DIAMOND. Nat Methods 12, 59–60 (2015). 
https://doi.org/10.1038/nmeth.3176
```

If you use the **profile** function with MMseqs2 (`-m mmseqs2`), also cite:
```
Steinegger, M., Söding, J. MMseqs2 enables sensitive protein sequence 
searching for the analysis of massive data sets. Nat Biotechnol 35, 
1026–1028 (2017). 
https://doi.org/10.1038/nbt.3988
```

## License

GPL-3.0 License

## Support

For questions, bug reports, or feature requests, please open an issue on GitHub:
https://github.com/cliffbueno/oxymetag/issues
