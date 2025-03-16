# data_prep.R - Run this script once to prepare data for the app

# Load required libraries
library(dplyr)

# Extract unique elements from final_pairs
driver_genes <- unique(final_pairs$driver_gene)
partner_genes <- unique(final_pairs$partner_gene)
cancer_types <- unique(final_pairs$cancer_types)

# Create directory for app data
dir.create("app/data", recursive = TRUE, showWarnings = FALSE)

# Save final_pairs directly
saveRDS(final_pairs, "app/data/final_pairs.rds")

# Subset mutation matrix - only rows for drivers in final_pairs
mutation_matrix_subset <- mutation_matrix[driver_genes, ]
saveRDS(mutation_matrix_subset, "app/data/mutation_matrix.rds")

# Subset dependency profiles
dependencyProfile_matrix_subset <- dependencyProfile_matrix[intersect(rownames(dependencyProfile_matrix), driver_genes), ]
saveRDS(dependencyProfile_matrix_subset, "app/data/dependencyProfile_matrix.rds")

# Subset LFC and p-value matrices
lfc_matrix_subset <- lfc_matrix[intersect(rownames(lfc_matrix), driver_genes), ]
pvalue_matrix_subset <- pvalue_matrix[intersect(rownames(pvalue_matrix), driver_genes), ]
saveRDS(lfc_matrix_subset, "app/data/lfc_matrix.rds")
saveRDS(pvalue_matrix_subset, "app/data/pvalue_matrix.rds")

# Cancer-specific matrices
mut_spec_lfcs_subset <- list()
mut_spec_pvalues_subset <- list()
mut_spec_dependencyProfiles_subset <- list()

for (cancer_type in cancer_types) {
  if (cancer_type == "pan") next
  
  # Get drivers for this cancer type - extract as character vector
  cancer_drivers <- as.character(unique(final_pairs[final_pairs$cancer_types == cancer_type, "driver_gene"]))
  
  if (length(cancer_drivers) > 0) {
    # Check if this cancer type exists in original data
    if (!cancer_type %in% names(mut_spec_lfcs)) next
    
    # Get drivers that exist in the matrix
    existing_drivers <- intersect(rownames(mut_spec_lfcs[[cancer_type]]), cancer_drivers)
    
    if (length(existing_drivers) > 0) {
      # Create subsets for this cancer type
      mut_spec_lfcs_subset[[cancer_type]] <- mut_spec_lfcs[[cancer_type]][existing_drivers, ]
      mut_spec_pvalues_subset[[cancer_type]] <- mut_spec_pvalues[[cancer_type]][existing_drivers, ]
      
      if (cancer_type %in% names(mut_spec_dependencyProfiles)) {
        mut_spec_dependencyProfiles_subset[[cancer_type]] <- mut_spec_dependencyProfiles[[cancer_type]][existing_drivers, ]
      }
    }
  }
}

saveRDS(mut_spec_lfcs_subset, "app/data/mut_spec_lfcs.rds")
saveRDS(mut_spec_pvalues_subset, "app/data/mut_spec_pvalues.rds")
saveRDS(mut_spec_dependencyProfiles_subset, "app/data/mut_spec_dependencyProfiles.rds")

# Save cancer_types_filt and common_genes
cancer_types_filt <- cancer_types
common_genes <- driver_genes
saveRDS(cancer_types_filt, "app/data/cancer_types_filt.rds")
saveRDS(common_genes, "app/data/common_genes.rds")

# If using druggable_genes, save them too
saveRDS(druggable_genes, "app/data/druggable_genes.rds")


# First, identify all cell lines (depmap_ids) needed for your analysis
needed_cell_lines <- unique(c(
  # Cell lines with mutations in your driver genes
  colnames(mutation_matrix)[colSums(mutation_matrix[unique(as.character(final_pairs$driver_gene)), ]) > 0],
  # Additional cell lines that might be referenced
  as.character(CRISPR$depmap_id[CRISPR$gene_name %in% c(
    unique(as.character(final_pairs$driver_gene)), 
    unique(as.character(final_pairs$partner_gene))
  )])
))

# Subset METADATA to only needed cell lines
METADATA_subset <- METADATA[METADATA$depmap_id %in% needed_cell_lines, ]

# Keep only essential columns from METADATA
essential_columns <- c("depmap_id", "cell_line", "stripped_cell_line_name", 
                       "primary_disease", "lineage_subtype")

METADATA_subset <- METADATA_subset[, intersect(names(METADATA_subset), essential_columns)]

saveRDS(METADATA_subset, "app/data/METADATA.rds")

# Get all genes needed (both drivers and partners)
needed_genes <- unique(c(
  as.character(final_pairs$driver_gene),
  as.character(final_pairs$partner_gene)
))

# Get needed cell lines (same as for METADATA)
# Already defined above: needed_cell_lines

# Subset CRISPR to only include relevant genes and cell lines
CRISPR_subset <- CRISPR[CRISPR$gene_name %in% needed_genes & 
                          CRISPR$depmap_id %in% needed_cell_lines, ]

# Keep only essential columns
CRISPR_subset <- CRISPR_subset[, c("depmap_id", "gene_name", "dependency")]

saveRDS(CRISPR_subset, "app/data/CRISPR.rds")

# For TPM data, we need:
# 1. Expression data for all genes that appear in our analysis
# 2. Only for relevant cell lines

# For pathway analysis, we might need more genes than just the final pairs
# Include any genes in the gene sets referenced by correlation_dfs (if available)
if (exists("correlation_dfs") && length(correlation_dfs) > 0) {
  # Extract gene sets mentioned in correlations
  referenced_gene_sets <- unique(unlist(lapply(correlation_dfs, function(df) {
    if ("gene_set" %in% colnames(df)) {
      return(unique(df$gene_set))
    }
    return(NULL)
  })))
  
  # Get genes in these gene sets
  genes_in_sets <- unique(unlist(lapply(gene_sets[referenced_gene_sets], function(genes) {
    unlist(genes)
  })))
  
  # Add these genes to our needed_genes list
  needed_genes <- unique(c(needed_genes, genes_in_sets))
}

# Subset TPM
TPM_subset <- TPM[TPM$gene_name %in% needed_genes & 
                    TPM$depmap_id %in% needed_cell_lines, ]

# Keep only essential columns
TPM_subset <- TPM_subset[, c("depmap_id", "cell_line", "gene_name", "rna_expression")]

saveRDS(TPM_subset, "app/data/TPM.rds")

saveRDS(correlation_dfs, "data/correlation_dfs.rds")

saveRDS(gene_sets, "app/data/gene_sets.rds")