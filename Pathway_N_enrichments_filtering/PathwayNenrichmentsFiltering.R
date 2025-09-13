# Setting working directory
setwd("C:/Users/aaisu/Box/Carter Lab/Projects/RNAseq_EthiopiaNUCI/Pathway_N_enrichments_filtering")

# Importing necessary libraries
library(tidyverse)

# R Script to Read Results from a Nested Directory Structure
# This version is designed to run from a script located in a sibling directory.

# --- 1. DEFINE TOP-LEVEL DIRECTORIES ---

# A character vector containing the relative paths to the main analysis folders.
# The '../' tells R to go up one level from the current working directory.
main_directories <- c("../GOenrichment_results", "../GSEA_results", "../KEGGS_results")


# --- 2. INITIALIZE THE FINAL LIST ---

# We'll create an empty list to store all the results.
all_results <- list()


# --- 3. LOOP THROUGH DIRECTORIES AND READ FILES ---

cat("Starting to read analysis results...\n")

# Outer loop: Iterates through each main directory path
for (main_dir_path in main_directories) {
  
  # Check if the main directory actually exists before proceeding
  if (!dir.exists(main_dir_path)) {
    cat("Warning: Main directory not found:", main_dir_path, "- Skipping.\n")
    next # Skips to the next main_dir_path
  }
  
  # Use basename() to get a clean name for messages and list keys
  main_dir_name <- basename(main_dir_path)
  cat("Processing directory:", main_dir_name, "\n")
  
  # Get a list of all subdirectories within the current main directory.
  sub_directories <- list.dirs(path = main_dir_path, full.names = TRUE, recursive = FALSE)
  
  # Create a temporary list to hold the data for the current main directory
  sub_dir_data_list <- list()
  
  # Middle loop: Iterates through the subdirectories found above
  for (sub_dir in sub_directories) {
    
    cat("  -> Reading from subdirectory:", basename(sub_dir), "\n")
    
    # Get a list of all files ending in .csv within the current subdirectory.
    csv_files <- list.files(path = sub_dir, pattern = "\\.csv$", full.names = TRUE)
    
    # Create another temporary list to hold the data frames for this subdirectory
    csv_data_list <- list()
    
    # Inner loop: Iterates through the CSV files found in the subdirectory
    for (csv_file in csv_files) {
      
      # Read the CSV file into a data frame
      current_df <- read.csv(csv_file)
      
      # Create a clean name for the data frame by removing the .csv extension
      df_name <- sub(pattern = "\\.csv$", replacement = "", x = basename(csv_file))
      
      # Add the data frame to our list for this subdirectory, using the clean name
      csv_data_list[[df_name]] <- current_df
    }
    
    # After reading all CSVs in a subdirectory, add that list of data frames
    # to the list for the main directory. Use the subdirectory's base name as the key.
    sub_dir_name <- basename(sub_dir)
    sub_dir_data_list[[sub_dir_name]] <- csv_data_list
  }
  
  # After processing all subdirectories, assign the collected data
  # to the final results list, using the main directory's clean base name as the key.
  all_results[[main_dir_name]] <- sub_dir_data_list
}

cat("...Finished reading all files.\n")


# --- 4. INSPECT THE RESULTS ---

# You can now inspect the structure of your new list
cat("\nStructure of the final 'all_results' list:\n")
str(all_results, max.level = 2)

# The way you access a specific data frame remains the same and is very clean:
# all_results[["KEGGS_results"]][["After_filtering"]][["YourFileName"]]


# --- 5. DEFINE INPUTS ---

# The path to your CSV file that contains the comparison pairs.
# This file MUST have three columns with headers "Set1", "Set2" and "Key".
combinations_file_path <- "comparisons.csv"


# --- 6. READ THE COMPARISON PLAN ---

# Check if the combinations file exists before trying to read it.
if (!file.exists(combinations_file_path)) {
  stop("Error: The combinations file was not found at the specified path: ", combinations_file_path)
}

combinations_df <- read.csv(combinations_file_path)


# --- 7. INITIALIZE RESULTS LIST ---

# Create an empty list to store the results of each comparison.
comparison_results_from_csv <- list()


# --- 8. LOOP, PARSE, AND COMPARE ---

cat("Starting automated comparisons based on", combinations_file_path, "\n")

# Loop through each row of the combinations data frame.
for (i in 1:nrow(combinations_df)) {
  
  # 4a. Get the paths for the current pair from the CSV
  set1_path <- combinations_df$Set1[i]
  set2_path <- combinations_df$Set2[i]
  
  # The name of the column in your data frames that contains the IDs to compare.
  # This is likely gene IDs, KEGG IDs, or GO term IDs.
  key_column <- combinations_df$Key[i] # <-- IMPORTANT: Change this if your ID column is named differently (e.g., "geneID")
  
  cat("\nProcessing pair", i, ":\n  Set1:", set1_path, "\n  Set2:", set2_path, "\n")
  
  # 4b. Parse the paths to access the data frames from the 'all_results' list
  # strsplit breaks the path into parts based on the "/" delimiter.
  parts1 <- strsplit(set1_path, "/")[[1]]
  parts2 <- strsplit(set2_path, "/")[[1]]
  
  # --- Retrieve DataFrames and Handle Potential Errors ---
  df1 <- NULL
  df2 <- NULL
  
  # A 'try' block will prevent the whole script from crashing if one path is invalid.
  try({
    df1 <- all_results[[parts1[1]]][[parts1[2]]][[parts1[3]]]
    df2 <- all_results[[parts2[1]]][[parts2[2]]][[parts2[3]]]
  }, silent = TRUE)
  
  # 4c. Validate that both data frames were successfully retrieved
  if (is.null(df1) || is.null(df2)) {
    cat("  -> ERROR: Could not retrieve one or both data frames. Please check paths in CSV. Skipping.\n")
    next # Skip to the next iteration
  }
  
  # 4d. Validate that the key column exists in both data frames
  if (!key_column %in% names(df1) || !key_column %in% names(df2)) {
    cat("  -> ERROR: Key column '", key_column, "' not found in one or both data frames. Skipping.\n", sep="")
    next # Skip to the next iteration
  }
  
  # 4e. Extract the vectors of IDs to compare
  ids1 <- df1[[key_column]]
  ids2 <- df2[[key_column]]
  
  # 4f. Perform the set operations
  intersection <- intersect(ids1, ids2)
  in_set1_only <- setdiff(ids1, ids2)
  in_set2_only <- setdiff(ids2, ids1)
  
  # 4g. Store the results
  # Create a clean name for this comparison's results
  comparison_name <- paste0(gsub("/", "_", set1_path), "___VS___", gsub("/", "_", set2_path))
  
  current_result <- list(
    intersection = intersection,
    in_set1_only = in_set1_only,
    in_set2_only = in_set2_only
  )
  
  comparison_results_from_csv[[comparison_name]] <- current_result
  cat("  -> Comparison successful. Results stored.\n")
}

cat("\n--- All comparisons complete. ---\n")


# --- 5. INSPECT THE FINAL RESULTS ---

# You can now inspect the structure of your new list
cat("\nStructure of the final 'comparison_results_from_csv' list:\n")
str(comparison_results_from_csv, max.level = 1)

