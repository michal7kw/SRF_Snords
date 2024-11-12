# Install required packages if needed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("GenomicFeatures", "AnnotationHub"))

library(GenomicFeatures)
library(AnnotationHub)
library(dplyr)

# Function to check transcript existence in different GENCODE versions
check_transcript_versions <- function(transcript_id) {
  # Connect to AnnotationHub
  ah <- AnnotationHub()
  
  # Query for GENCODE datasets
  gencode_queries <- query(ah, c("GENCODE", "Human"))
  
  # Get different versions of GENCODE
  gencode_versions <- gencode_queries[grep("TxDb", gencode_queries$title)]
  
  # Check each version
  results <- data.frame(version=character(),
                        found=logical(),
                        coordinates_match=logical(),
                        stringsAsFactors=FALSE)
  
  for(i in seq_along(gencode_versions)) {
    tryCatch({
      txdb <- gencode_versions[[i]]
      
      # Check if transcript exists
      tx <- transcripts(txdb, filter=list(tx_name=transcript_id))
      
      version_name <- names(gencode_versions)[i]
      found <- length(tx) > 0
      
      # Add to results
      results <- rbind(results,
                       data.frame(version=version_name,
                                  found=found,
                                  stringsAsFactors=FALSE))
    }, error=function(e) {
      message("Error processing version ", names(gencode_versions)[i])
    })
  }
  
  return(results)
}

# Example usage with your transcript
results <- check_transcript_versions("ENST00000384335.1")
print(results)

# To check multiple transcripts from your data
transcripts_to_check <- c("ENST00000384335.1", 
                          "ENST00000656407.1",
                          "ENST00000581463.1")

all_results <- lapply(transcripts_to_check, check_transcript_versions)
names(all_results) <- transcripts_to_check