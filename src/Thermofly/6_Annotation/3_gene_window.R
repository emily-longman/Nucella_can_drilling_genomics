library(tidyverse)

# Read in arguments and set working directory
args <- commandArgs(trailingOnly = TRUE)
k <- args[1]
setwd(k)

# Load gene positions
gene_positions <- read_tsv("/gpfs1/home/d/s/dsadler1/scratch/SeaUrchinEnvPop/GDS/gene_positions.txt", col_names = c("chr", "start", "end", "gene_id"))

# Adjust gene positions to include +/- 1kb window
gene_positions <- gene_positions %>%
  mutate(window_start = pmax(0, start - 1000),  # Ensure positions don't go below 0
         window_end = end + 1000)

# List of RDS files to process
file_list <- list.files(path = k, pattern = "\\.rds$", full.names = TRUE)

# Initialize a variable to store the aggregated ratios
all_ratios <- NULL

# Function to update summary statistics incrementally
update_summary_stats <- function(existing_stats, new_ratios) {
  if (is.null(existing_stats)) {
    return(new_ratios)
  } else {
    return(rbind(existing_stats, new_ratios))
  }
}

# Process each chunk file
for (file in file_list) {
  # Extract chromosome name from file (adjust regex if filenames differ)
  current_chr <- str_extract(basename(file), "NW_\\d+\\.\\d+")

  # Filter gene list for the current chromosome
  current_genes <- gene_positions %>% filter(chr == current_chr)
  
  if (nrow(current_genes) == 0) {
    next  # Skip if no genes for the current chromosome
  }
  
  # Load the data
  chunktest <- readRDS(file)
  
  # Merge data with gene windows
  data_with_windows <- chunktest %>%
    mutate(pos = (start + end) / 2) %>%  # Adjust to actual position column names in your data
    inner_join(current_genes, by = c("chr")) %>%
    filter(pos >= window_start & pos <= window_end) %>%  # Filter positions within the window
    mutate(gene_window = gene_id)  # Assign the gene ID for grouping
  
  # Process data to calculate AIC differences and ratios
  model_aic <- data_with_windows %>%
    mutate(AICdiff = if_else(AIC_wea < AIC_null, 1, 0)) %>%
    arrange(desc(AICdiff))
  
  model_aic_1 <- model_aic %>%
    group_by(gene_window, variable, test_code) %>%
    summarise(rr = sum(AICdiff == 1, na.rm = TRUE), .groups = "drop")
  
  real_data <- model_aic_1 %>% filter(test_code == 0)
  perm_data <- model_aic_1 %>% filter(test_code > 0)
  
  ratios <- perm_data %>%
    left_join(real_data, by = c("gene_window", "variable"), suffix = c("_perm", "_real")) %>%
    mutate(rr_ratio = rr_real / rr_perm) %>%
    select(gene_window, variable, test_code_perm, rr_ratio)
  
  # Incrementally update the overall ratios
  all_ratios <- update_summary_stats(all_ratios, ratios)
  
  # Clear memory and proceed to the next chunk
  rm(chunktest, data_with_windows, model_aic, model_aic_1, real_data, perm_data, ratios)
  gc()  # Force garbage collection to free memory
}

# Save aggregated ratios
save(all_ratios, file = "all_ratios_genewindows2.RData")

# Calculate final summary statistics across all chunks
summary <- all_ratios %>%
  group_by(gene_window, variable) %>%
  drop_na() %>%
  summarise(mean = mean(rr_ratio),
            sd = sd(rr_ratio),
            n = n(),  # Number of observations
            se = sd / sqrt(n),  # Calculate standard error
            ci_low = mean - 1.96 * se,
            ci_high = mean + 1.96 * se,
            quantile_0.05 = quantile(rr_ratio, 0.05),
            quantile_0.95 = quantile(rr_ratio, 0.95),
            quantile_0.5 = quantile(rr_ratio, 0.5)) %>%
  ungroup()

# View summary statistics
print(summary)

# Save summary statistics
save(summary, file = "summary_genewindows_2.RData")
