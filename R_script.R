# Load required libraries
library(DEXSeq)
library(DESeq2)
library(dplyr)
library(ggplot2)

# Set up the working directory
data_dir <- "/beegfs/scratch/ric.broccoli/ric.broccoli/PW_RNA_seq_deep"
working_dir <- "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Snords"
setwd(working_dir)

# Local Functions

compare_file_structures <- function(file_list) {
  results <- lapply(file_list, function(file_path) {
    n_lines <- length(readLines(file_path))
    file_size <- file.info(file_path)$size
    list(
      file = basename(file_path),
      n_lines = n_lines,
      size_mb = file_size / 1024 / 1024,
      avg_line_size = if (n_lines > 0) file_size / n_lines else 0
    )
  })
  do.call(rbind, results) %>% as.data.frame()
}

parse_dexseq_id <- function(feature_id) {
  if (startsWith(feature_id, "*")) {
    return(list(gene_id = feature_id, exon_num = NA))
  }
  tryCatch({
    parts <- strsplit(feature_id, ":")[[1]]
    gene_id <- parts[1]
    exon_num <- gsub('"', '', parts[2])
    list(gene_id = gene_id, exon_num = exon_num)
  }, error = function(e) {
    list(gene_id = feature_id, exon_num = NA)
  })
}

examine_dexseq_file <- function(file_path, n_head = 5, n_random = 5, n_tail = 5) {
  cat("\nExamining file:", basename(file_path), "\n")
  
  # Peek at raw file contents
  cat("\nFirst few lines (raw):\n")
  system(paste("head -n 5", file_path))
  
  # Read the file
  df <- tryCatch({
    read.table(file_path, sep = "\t", header = FALSE, 
               col.names = c("feature_id", "count"), stringsAsFactors = FALSE)
  }, error = function(e) {
    cat("Error reading file:", conditionMessage(e), "\n")
    return(NULL)
  })
  
  if (is.null(df)) return(NULL)
  
  # Process feature_id
  parsed_ids <- lapply(df$feature_id, parse_dexseq_id)
  df$gene_id <- sapply(parsed_ids, `[[`, "gene_id")
  df$exon_number <- sapply(parsed_ids, `[[`, "exon_num")
  
  # File information
  cat("\nFile size:", file.info(file_path)$size / 1024, "KB\n")
  cat("Number of lines:", nrow(df), "\n")
  
  # Basic DataFrame information
  cat("\nDataFrame Info:\n")
  str(df)
  
  # Print column names
  cat("\nColumn names:\n")
  print(colnames(df))
  
  # Show first few lines
  cat("\nFirst", n_head, "lines:\n")
  print(head(df, n_head))
  
  # Basic statistics for counts
  cat("\nCount statistics:\n")
  print(summary(df$count))
  
  # Additional information
  cat("\nNumber of unique genes:", length(unique(df$gene_id)), "\n")
  cat("\nTop 5 genes by total counts:\n")
  gene_counts <- aggregate(count ~ gene_id, data = df, sum)
  print(head(gene_counts[order(-gene_counts$count), ], 5))
  
  # Check for potential issues
  cat("\nChecking for potential issues:\n")
  issues <- character()
  if (any(is.na(df$count))) issues <- c(issues, "- Contains missing values in counts")
  if (any(df$count < 0)) issues <- c(issues, "- Contains negative counts")
  if (!is.numeric(df$count)) issues <- c(issues, "- Counts are not numeric")
  if (!is.character(df$feature_id)) issues <- c(issues, paste("- Feature ID column is not string type (current type:", class(df$feature_id), ")"))
  
  # Check for special entries
  special_entries <- df[startsWith(df$feature_id, "*"), ]
  if (nrow(special_entries) > 0) {
    cat("\nSpecial entries found:\n")
    print(special_entries)
  }
  
  if (length(issues) > 0) {
    cat("\nIssues found:\n")
    cat(paste(issues, collapse = "\n"), "\n")
  } else {
    cat("No major issues found\n")
  }
  
  return(df)
}

examine_genes_of_interest <- function(df, genes) {
  cat("\nExamining genes of interest:\n")
  for (gene in genes) {
    gene_data <- df[grepl(gene, df$gene_id, fixed = TRUE), ]
    if (nrow(gene_data) > 0) {
      cat("\n", gene, ":\n")
      cat("Number of exons:", nrow(gene_data), "\n")
      cat("Total counts:", sum(gene_data$count), "\n")
      cat("Mean counts per exon:", mean(gene_data$count), "\n")
      cat("\nExon-level data:\n")
      print(gene_data[order(gene_data$exon_number), ])
      
      # Create plot for this gene
      p <- ggplot(gene_data, aes(x = seq_along(count), y = count)) +
        geom_bar(stat = "identity") +
        ggtitle(paste('Exon counts for', gene)) +
        xlab('Exon number') +
        ylab('Counts') +
        theme_minimal()
      
      ggsave(paste0('exon_counts_', gene, '.pdf'), p, width = 10, height = 4)
    }
  }
}

prepare_sample_table <- function(sample_info) {
  samples_subset <- sample_info[sample_info$condition %in% c('EDO', 'ND1'), ]
  samples_subset$condition <- factor(samples_subset$condition)
  samples_subset$replicate <- factor(samples_subset$replicate)
  rownames(samples_subset) <- NULL
  samples_subset
}

create_dexseq_object <- function(sample_table) {
  formula <- ~ sample + exon + condition:exon
  reduced_formula <- ~ sample + exon
  
  dxd <- DEXSeqDataSetFromHTSeq(
    countfiles = sample_table$count_file,
    sampleData = sample_table,
    design = formula,
    flattenedfile = NULL
  )
  
  dxd
}

run_dexseq_analysis <- function(dxd) {
  dxd <- estimateSizeFactors(dxd)
  dxd <- estimateDispersions(dxd)
  dxd <- testForDEU(dxd)
  dxd <- estimateExonFoldChanges(dxd)
  dxd
}

extract_results <- function(dxd) {
  res <- DEXSeqResults(dxd)
  results_df <- as.data.frame(res)
  results_df$significant <- results_df$padj < 0.05
  results_df
}

filter_significant_results <- function(results_df) {
  significant_results <- results_df[
    results_df$padj < 0.05 & 
    results_df$log2fold_EDO_ND1 != 0,
  ]
  significant_results[order(significant_results$padj), ]
}

# DEXSeq Files

# Load Files
dexseq_dir <- file.path(working_dir, "DexSeq_counts")
count_files <- list.files(dexseq_dir, pattern = "\\.formatted\\.counts$", full.names = TRUE)

# Create sample information
sample_info <- data.frame(
  sample = sub("\\.formatted\\.counts$", "", basename(count_files)),
  condition = ifelse(startsWith(basename(count_files), "EDO"), "EDO",
                     ifelse(startsWith(basename(count_files), "ND1"), "ND1",
                            ifelse(startsWith(basename(count_files), "PW1"), "PW1", "Unknown"))),
  replicate = sub(".*_(\\d+)\\.formatted\\.counts$", "\\1", basename(count_files)),
  count_file = count_files,
  stringsAsFactors = FALSE
)

# Filter for EDO and ND1 samples (excluding PW1)
edo_nd1_samples <- sample_info[sample_info$condition %in% c('EDO', 'ND1'), ]
edo_nd1_samples$condition <- factor(edo_nd1_samples$condition)
edo_nd1_samples$replicate <- factor(edo_nd1_samples$replicate)
rownames(edo_nd1_samples) <- NULL

# Examine Files
cat("Comparing file structures:\n")
structure_comparison <- compare_file_structures(edo_nd1_samples$count_file)
print(structure_comparison)

# Examine first file in detail
first_file <- edo_nd1_samples$count_file[1]
df <- examine_dexseq_file(first_file)

# Plot overall count distribution
pdf("count_distribution.pdf", width = 15, height = 5)
par(mfrow = c(1, 3))
counts <- df$count[df$count > 0]
hist(counts, breaks = 50, col = 'skyblue', border = 'black', 
     xlim = c(0, quantile(counts, 0.99)), main = 'Count Distribution (> 0)',
     xlab = 'Counts', ylab = 'Frequency')
dev.off()

# Print outliers
outliers <- df[df$count > quantile(counts, 0.99), ]
cat("\nOutliers (counts above 99th percentile):\n")
print(head(outliers[, c('feature_id', 'count', 'gene_id')]))

# Plot outliers distribution
pdf("outlier_distribution.pdf", width = 15, height = 5)
par(mfrow = c(1, 3))
outlier_counts <- outliers$count
hist(outlier_counts, breaks = 50, col = 'salmon', border = 'black',
     main = 'Outlier Count Distribution', xlab = 'Counts', ylab = 'Frequency')
hist(log10(outlier_counts), breaks = 50, col = 'lightgreen', border = 'black',
     main = 'Log10 Outlier Count Distribution', xlab = 'Log10(Counts)', ylab = 'Frequency')
boxplot(outlier_counts, main = 'Outlier Count Boxplot', ylab = 'Counts')
dev.off()

# DEXSeq Analysis
create_dexseq_dataset <- function(sample_info, formatted_files) {
  tryCatch({
    sample_data <- sample_info[, c('sample', 'condition')]
    sample_data$sample <- as.character(sample_data$sample)
    sample_data$condition <- as.character(sample_data$condition)
    
    ordered_files <- sample_info$count_file
    
    cat("\nSample data for DEXSeq:\n")
    print(sample_data)
    
    cat("\nOrdered count files:\n")
    for (i in seq_along(sample_data$sample)) {
      cat(sample_data$sample[i], ": ", ordered_files[i], "\n")
    }
    
    gff_file <- "gencode.v31.basic.annotation.DEXSeq.gff"
    if (!file.exists(gff_file)) {
      cat("GFF file not found:", gff_file, "\n")
    } else {
      cat("GFF file exists:", gff_file, "\n")
      cat(readLines(gff_file, n = 5), sep = "\n")
    }
    
    dxd <- DEXSeqDataSetFromHTSeq(
      countfiles = ordered_files,
      sampleData = sample_data,
      design = ~ condition,
      flattenedfile = gff_file
    )
    
    return(dxd)
    
  }, error = function(e) {
    cat("Error creating DEXSeqDataSet:", conditionMessage(e), "\n")
    stop(e)
  })
}

dxd <- create_dexseq_dataset(sample_info, count_files)
