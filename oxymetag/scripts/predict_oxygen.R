#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(rlang)
  library(mgcv)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 7) {
  stop("Usage: Rscript predict_oxygen.R <input_dir> <output_file> <package_data_dir> <mode> <idcut> <ecut> <bitcut>")
}

input_dir <- args[1]
output_file <- args[2] 
package_data_dir <- args[3]
mode <- args[4]
idcut <- as.numeric(args[5])
ecut <- as.numeric(args[6])
bitcut <- as.numeric(args[7])

predict_oxygen <- function(input_dir, output_file, package_data_dir, mode, idcut, ecut, bitcut) {
  
  # Load package data files
  map_file <- file.path(package_data_dir, "pfam_headers_table.txt")
  lengths_file <- file.path(package_data_dir, "pfam_lengths.tsv")
  model_file <- file.path(package_data_dir, "oxygen_model.rds")
  pfams_file <- file.path(package_data_dir, "Oxygen_pfams.csv")
  
  # Check if package data files exist
  if (!file.exists(map_file)) stop(paste("Package data file not found:", map_file))
  if (!file.exists(lengths_file)) stop(paste("Package data file not found:", lengths_file))
  if (!file.exists(model_file)) stop(paste("Package data file not found:", model_file))
  if (!file.exists(pfams_file)) stop(paste("Package data file not found:", pfams_file))
  
  # Load package data
  map <- read.table(map_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "") %>%
    separate(Header, into = c("Header", "Junk"), sep = " ") %>%
    select(-Junk) %>%
    filter(!duplicated(Header))
  
  pfam_gene_length <- read.delim(lengths_file)
  
  # Load the trained model and oxygen classifications
  oxygen_model <- readRDS(model_file)
  oxygen_pfams <- read.csv(pfams_file, stringsAsFactors = FALSE)
  aerobic_pfams <- oxygen_pfams %>% filter(Oxygen == "aerobic")
  anaerobic_pfams <- oxygen_pfams %>% filter(Oxygen == "anaerobic")
  
  # Determine method and file pattern based on mode
  if (mode == "modern") {
    method <- "diamond"
    file_pattern <- "*_diamond.tsv"
  } else if (mode == "ancient") {
    method <- "mmseqs2"
    file_pattern <- "*_mmseqs.tsv"
  } else if (mode == "custom") {
    # For custom mode, try to detect which files exist
    diamond_files <- list.files(input_dir, pattern = "*_diamond.tsv", full.names = TRUE)
    mmseqs_files <- list.files(input_dir, pattern = "*_mmseqs.tsv", full.names = TRUE)
    
    if (length(diamond_files) > 0 && length(mmseqs_files) == 0) {
      method <- "diamond"
      file_pattern <- "*_diamond.tsv"
    } else if (length(mmseqs_files) > 0 && length(diamond_files) == 0) {
      method <- "mmseqs2"
      file_pattern <- "*_mmseqs.tsv"
    } else if (length(diamond_files) > 0 && length(mmseqs_files) > 0) {
      stop("Both DIAMOND and MMseqs2 files found. Please specify input directory with only one type.")
    } else {
      stop("No profiling output files found in input directory.")
    }
  } else {
    stop(paste("Unknown mode:", mode))
  }
  
  files <- list.files(input_dir, pattern = file_pattern, full.names = TRUE)
  
  if (length(files) == 0) {
    stop(paste("No", file_pattern, "files found in", input_dir))
  }
  
  message("Processing ", length(files), " files using ", method, " output format")
  
  # Initialize results dataframe
  results <- data.frame(
    SampleID = character(length(files)),
    ratio = numeric(length(files)),
    aerobe_pfams = integer(length(files)),
    anaerobe_pfams = integer(length(files)),
    Per_aerobe = numeric(length(files)),
    stringsAsFactors = FALSE
  )
  
  # Process each file
  for (i in 1:length(files)) {
    
    # Extract sample ID from filename
    sample_id <- basename(files[i])
    sample_id <- gsub("_diamond.tsv$|_mmseqs.tsv$", "", sample_id)
    results$SampleID[i] <- sample_id
    
    # Read and filter profiling output
    if (file.size(files[i]) == 0) {
      message("Warning: Empty file ", files[i])
      next
    }
    
    # Read data based on method
    if (method == "diamond") {
      d <- read.table(files[i], stringsAsFactors = FALSE) %>%
        set_names(c("qseqid", "sseqid", "pident", "length", "qstart", "qend", 
                    "sstart", "send", "evalue", "bitscore")) %>%
        left_join(map, by = c("sseqid" = "Header")) %>%
        group_by(qseqid) %>%
        slice_max(bitscore, n = 1, with_ties = FALSE) %>%
        ungroup()
      
      # Apply filtering based on mode
      if (mode == "modern") {
        d <- d %>% filter(pident >= 60, evalue < 0.001, bitscore >= 50)
      } else if (mode == "custom") {
        d <- d %>% filter(pident >= idcut, evalue < ecut, bitscore >= bitcut)
      }
      
    } else if (method == "mmseqs2") {
      d <- read.table(files[i], stringsAsFactors = FALSE) %>%
        set_names(c("query", "target", "pident", "length", "mismatch", "gapopen", 
                    "qstart", "qend", "tstart", "tend", "evalue", "bitscore", 
                    "qlen", "tlen", "cigar", "qaln", "taln")) %>%
        mutate(pident = pident * 100) %>%
        left_join(map, by = c("target" = "Header")) %>%
        group_by(query) %>%
        slice_max(bitscore, n = 1, with_ties = FALSE) %>%
        ungroup() %>%
        filter(Pfam %in% oxygen_pfams$Pfam)
      
      # Apply filtering based on mode
      if (mode == "ancient") {
        d <- d %>% filter(pident >= 86, evalue < 0.001, bitscore >= 50)
      } else if (mode == "custom") {
        d <- d %>% filter(pident >= idcut, evalue < ecut, bitscore >= bitcut)
      }
    }
    
    if (nrow(d) == 0) {
      message("Warning: No significant hits for ", sample_id)
      results$ratio[i] <- 0
      results$aerobe_pfams[i] <- 0
      results$anaerobe_pfams[i] <- 0
      next
    }
    
    # Count pfams
    pf_count <- as.data.frame(table(d$Pfam))
    results$aerobe_pfams[i] <- sum(as.character(pf_count$Var1) %in% aerobic_pfams$Pfam)
    results$anaerobe_pfams[i] <- sum(as.character(pf_count$Var1) %in% anaerobic_pfams$Pfam)
    
    # Calculate gene hits and length correction
    gene.hits <- d %>% 
      group_by(Pfam) %>% 
      summarise(total_count = n(), .groups = 'drop')
    
    gene.hit.length.correction <- gene.hits %>%
      left_join(pfam_gene_length, by = "Pfam") %>%
      mutate(RPK = total_count / (Gene.length/1000)) %>%
      left_join(oxygen_pfams, by = "Pfam")
    
    # Sum by oxygen type
    oxygen_rpk <- gene.hit.length.correction %>%
      group_by(Oxygen) %>%
      summarize(RPKsum = sum(RPK, na.rm = TRUE), .groups = 'drop')
    
    # Calculate ratio (aerobe/anaerobe)
    aerobe_rpk <- oxygen_rpk$RPKsum[oxygen_rpk$Oxygen == "aerobic"]
    anaerobe_rpk <- oxygen_rpk$RPKsum[oxygen_rpk$Oxygen == "anaerobic"]
    
    if (length(anaerobe_rpk) == 0 || anaerobe_rpk == 0) {
      results$ratio[i] <- ifelse(length(aerobe_rpk) > 0 && aerobe_rpk > 0, Inf, 0)
    } else {
      results$ratio[i] <- aerobe_rpk / anaerobe_rpk
    }
    
    message("Processed sample ", i, "/", length(files), ": ", sample_id)
  }
  
  # Make predictions using the GAM model
  new_data <- data.frame(ratio = results$ratio)
  results$Per_aerobe <- predict(oxygen_model, newdata = new_data, type = "response")
  
  # Constrain predictions to 0-100% and set to 100% if ratio > 35
  results <- results %>%
    mutate(Per_aerobe = pmax(0, pmin(100, Per_aerobe))) %>%
    mutate(Per_aerobe = ifelse(ratio > 35, 100, Per_aerobe))
  
  # Save results
  write.table(results, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
  message("Results saved to: ", output_file)
  
  return(results)
}

# Run the function
predict_oxygen(input_dir, output_file, package_data_dir, mode, idcut, ecut, bitcut)
