# Install required packages if needed
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(c("GenomicFeatures", "AnnotationHub", "rtracklayer", "GenomicRanges"))

library(GenomicFeatures)
library(AnnotationHub)
library(rtracklayer)
library(GenomicRanges)
library(dplyr)

# Function to create GRanges object from your peak data
create_peak_granges <- function(data) {
  # Ensure all required columns are numeric and not NA
  data$start <- as.numeric(data$start)
  data$end <- as.numeric(data$end)
  
  # Remove any rows with NA in essential columns
  data <- data[!is.na(data$start) & !is.na(data$end) & !is.na(data$seqnames), ]
  
  # Create GRanges object
  gr <- GRanges(
    seqnames = data$seqnames,
    ranges = IRanges(
      start = data$start,
      end = data$end
    ),
    transcript_id = data$transcriptId,
    gene_id = data$ENSEMBL,
    symbol = data$SYMBOL
  )
  
  return(gr)
}

# Read the Peaks.csv file
peaks_data <- read.csv("Peaks.csv", stringsAsFactors = FALSE)

# Check data structure
print("Data structure:")
str(peaks_data)

# Print first few rows to verify
print("First few rows:")
head(peaks_data)

# Create GRanges object with error handling
peaks_granges <- tryCatch({
  create_peak_granges(peaks_data)
}, error = function(e) {
  print(paste("Error creating GRanges:", e$message))
  return(NULL)
})

# Function to check specific GENCODE versions
check_specific_versions <- function(peaks_granges, versions = c("v19", "v24", "v29", "v31", "v32", "v38", "v43")) {
  if (is.null(peaks_granges)) {
    return(NULL)
  }
  
  # Initialize AnnotationHub
  ah <- AnnotationHub()
  
  # Create results storage
  results_list <- list()
  
  for(version in versions) {
    message("Checking GENCODE ", version)
    
    # Query for specific GENCODE version
    query_string <- paste0("GENCODE", version)
    gencode_query <- query(ah, c(query_string, "GTF", "Human", "GRCh38"))
    
    # Get GTF file
    tryCatch({
      gtf_id <- gencode_query$ah_id[1]
      gtf_data <- ah[[gtf_id]]
      
      # Extract transcript information
      transcripts <- gtf_data[gtf_data$type == "transcript"]
      
      # Check overlaps and matches
      matches <- findOverlaps(peaks_granges, transcripts)
      
      # Compile results
      results <- data.frame(
        peak_id = seq_along(peaks_granges),
        transcript_id = peaks_granges$transcript_id,
        gene_id = peaks_granges$gene_id,
        symbol = peaks_granges$symbol,
        gencode_transcript_id = transcripts$transcript_id[subjectHits(matches)],
        gencode_gene_id = transcripts$gene_id[subjectHits(matches)],
        gencode_gene_name = transcripts$gene_name[subjectHits(matches)]
      )
      
      # Add to results list
      results_list[[version]] <- results
      
    }, error = function(e) {
      message("Error processing GENCODE ", version, ": ", e$message)
    })
  }
  
  return(results_list)
}

# Run the analysis if peaks_granges was created successfully
if (!is.null(peaks_granges)) {
  gencode_results <- check_specific_versions(peaks_granges)
  
  # Function to generate detailed report
  generate_report <- function(results_list) {
    if (is.null(results_list)) {
      cat("No results to report - check data loading issues\n")
      return()
    }
    
    cat("Summary of GENCODE Version Comparison\n")
    cat("=====================================\n\n")
    
    for(version in names(results_list)) {
      cat(sprintf("\nMatches in GENCODE %s:\n", version))
      if (nrow(results_list[[version]]) > 0) {
        print(table(results_list[[version]]$transcript_id %in% 
                      results_list[[version]]$gencode_transcript_id))
        
        # Add transcript ID details
        cat("\nTranscript matches:\n")
        matched_transcripts <- results_list[[version]][
          results_list[[version]]$transcript_id %in% 
            results_list[[version]]$gencode_transcript_id,
        ]
        print(head(matched_transcripts[, c("transcript_id", "gene_id", "symbol")]))
      } else {
        cat("No matches found\n")
      }
    }
  }
  
  # Generate report
  generate_report(gencode_results)
  
  # Save results if they exist
  if (!is.null(gencode_results)) {
    # Save as CSV for easier viewing
    for(version in names(gencode_results)) {
      write.csv(gencode_results[[version]], 
                file = paste0("gencode_", version, "_comparison.csv"),
                row.names = FALSE)
    }
    # Save full results object
    saveRDS(gencode_results, "gencode_version_comparison.rds")
  }
} else {
  message("Analysis couldn't be performed due to data loading issues")
}