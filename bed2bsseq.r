#!/usr/bin/env Rscript

# DNA Methylation Analysis Script in R
# 
# This script processes multiple BED files containing methylation data and
# outputs a combined BSseq object in R format.
#
# Usage:
#   Rscript methylation_analysis.R --input file1.bed,file2.bed --total-col 5 --meth-col 6 --output bsseq_data.rds

library(bsseq)
library(data.table)
library(argparse)

# Parse command line arguments
parser <- ArgumentParser(description='Process methylation data from BED files into a BSseq object')
parser$add_argument('--input', required=TRUE, help='Comma-separated list of input BED files')
parser$add_argument('--total-col', type='integer', required=TRUE, help='Column number for total reads (1-based)')
parser$add_argument('--meth-col', type='integer', required=TRUE, help='Column number for methylated reads (1-based)')
parser$add_argument('--output', required=TRUE, help='Output file for BSseq object (RDS format)')

args <- parser$parse_args()

# Split the input string into individual file paths
bed_files <- unlist(strsplit(args$input, ","))
sample_names <- basename(tools::file_path_sans_ext(bed_files))

# Function to read and process a BED file
read_bed_file <- function(file_path, total_col, meth_col) {
  # Set column names for the first 3 columns (chr, start, end)
  col_names <- c("chr", "start", "end")
  
  # Read the file using fread for speed
  bed_data <- fread(file_path, header = FALSE, select = c(1:3, total_col, meth_col))
  
  # Set column names
  colnames(bed_data)[1:3] <- c("chr", "start", "end")
  colnames(bed_data)[4:5] <- c("total", "meth")
  
  # Return data.table
  return(bed_data)
}

# Process all BED files
cat("Processing BED files...\n")
all_sites <- list()
all_M <- list()
all_Cov <- list()

for (i in 1:length(bed_files)) {
  file_path <- bed_files[i]
  sample_name <- sample_names[i]
  
  cat(sprintf("Processing %s...\n", file_path))
  
  # Read and process the BED file
  bed_data <- read_bed_file(file_path, args$total_col, args$meth_col)
  
  # Store the data
  all_sites[[i]] <- bed_data[, .(chr, start, end)]
  all_M[[i]] <- bed_data$meth
  all_Cov[[i]] <- bed_data$total
  
  cat(sprintf("Found %d methylation positions in %s\n", nrow(bed_data), sample_name))
}

# Combine all sites and create a unique set
cat("Combining sites from all samples...\n")
combined_sites <- rbindlist(all_sites)
unique_sites <- unique(combined_sites, by = c("chr", "start"))
setkey(unique_sites, chr, start)

# Create matrices for BSseq object
n_sites <- nrow(unique_sites)
n_samples <- length(bed_files)

M <- matrix(0, nrow = n_sites, ncol = n_samples)
Cov <- matrix(0, nrow = n_sites, ncol = n_samples)

# Fill the matrices with data from each sample
cat("Creating methylation matrices...\n")
for (i in 1:length(bed_files)) {
  sample_sites <- all_sites[[i]]
  setkey(sample_sites, chr, start)
  
  # Find indices of sample sites in the unique sites
  site_indices <- match(
    paste(sample_sites$chr, sample_sites$start), 
    paste(unique_sites$chr, unique_sites$start)
  )
  
  # Fill in the M and Cov matrices
  M[site_indices, i] <- all_M[[i]]
  Cov[site_indices, i] <- all_Cov[[i]]
}

# Create the BSseq object
cat("Creating BSseq object...\n")
bs <- BSseq(
  M = M, 
  Cov = Cov,
  pos = unique_sites$start,
  chr = unique_sites$chr,
  sampleNames = sample_names
)

# Save the BSseq object
cat(sprintf("Saving BSseq object to %s...\n", args$output))
saveRDS(bs, file = args$output)

cat("Done!\n")
