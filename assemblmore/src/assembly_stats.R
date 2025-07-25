#!/usr/bin/env Rscript

# assembly_stats.R - Custom assembly statistics script
# Computes contiguity, completeness, and correctness metrics for genome assemblies

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(scales)
  library(dbscan)  # For DBSCAN clustering
})

# Function to read FASTA file and get sequence lengths
read_fasta_lengths <- function(fasta_file) {
  if (!file.exists(fasta_file)) {
    stop(paste("FASTA file does not exist:", fasta_file))
  }
  
  # Read FASTA file line by line
  lines <- readLines(fasta_file)
  
  # Initialize variables
  lengths <- c()
  current_length <- 0
  
  for (line in lines) {
    if (startsWith(line, ">")) {
      # Header line - save previous sequence length if any
      if (current_length > 0) {
        lengths <- c(lengths, current_length)
        current_length <- 0
      }
    } else {
      # Sequence line - add to current length
      current_length <- current_length + nchar(gsub("[^ACGTNacgtn]", "", line))
    }
  }
  
  # Don't forget the last sequence
  if (current_length > 0) {
    lengths <- c(lengths, current_length)
  }
  
  return(lengths)
}

# Function to calculate NX values
calculate_nx <- function(lengths, x_values = seq(1, 100, length.out = 10000)) {
  # Sort lengths in descending order
  sorted_lengths <- sort(lengths, decreasing = TRUE)
  total_length <- sum(sorted_lengths)
  
  nx_values <- data.frame(
    X = x_values,
    NX = numeric(length(x_values))
  )
  
  cumulative_length <- 0
  length_index <- 1
  
  for (i in seq_along(x_values)) {
    target_length <- total_length * (x_values[i] / 100)
    
    while (cumulative_length < target_length && length_index <= length(sorted_lengths)) {
      cumulative_length <- cumulative_length + sorted_lengths[length_index]
      length_index <- length_index + 1
    }
    
    nx_values$NX[i] <- if (length_index > 1) sorted_lengths[length_index - 1] else 0
  }
  
  return(nx_values)
}

# Function to calculate auN (area under the NX curve)
calculate_aun <- function(lengths) {
  sum_li_squared <- sum(lengths^2)
  sum_li <- sum(lengths)
  aun <- sum_li_squared / sum_li
  return(aun)
}

# Function to create NX plot (single assembly)
create_nx_plot <- function(nx_data, assembly_name, output_file) {
  p <- ggplot(nx_data, aes(x = X, y = NX)) +
    geom_line(color = "blue", linewidth = 1.5) +
    scale_y_continuous(labels = comma_format()) +
    labs(
      title = paste("NX Plot for", assembly_name),
      x = "X (%)",
      y = "NX (bp)",
      subtitle = paste("Assembly contiguity metrics")
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    ) +
    geom_hline(yintercept = nx_data$NX[nx_data$X == 50], 
               linetype = "dashed", color = "red", alpha = 0.7) +
    annotate("text", 
             x = 25, 
             y = nx_data$NX[nx_data$X == 50] * 1.1, 
             label = paste("N50 =", comma(nx_data$NX[nx_data$X == 50])), 
             color = "red")
  
  ggsave(output_file, plot = p, width = 10, height = 6, dpi = 300)
  cat("NX plot saved to:", output_file, "\n")
}

# Function to create comparative NX plot (multiple assemblies)
create_comparative_nx_plot <- function(nx_data_list, output_file) {
  # Combine all NX data with assembly names
  combined_data <- data.frame()
  assembly_names <- names(nx_data_list)
  colors <- rainbow(length(nx_data_list))
  names(colors) <- assembly_names  # Name the colors vector to ensure proper mapping
  
  for (i in seq_along(nx_data_list)) {
    assembly_name <- names(nx_data_list)[i]
    nx_data <- nx_data_list[[i]]
    nx_data$Assembly <- assembly_name
    combined_data <- rbind(combined_data, nx_data)
  }
  
  p <- ggplot(combined_data, aes(x = X, y = NX, color = Assembly)) +
    geom_line(linewidth = 1.2) +
    scale_y_continuous(labels = comma_format()) +
    scale_color_manual(values = colors) +
    labs(
      title = "Comparative NX Plot",
      x = "X (%)",
      y = "NX (bp)",
      subtitle = "Assembly contiguity comparison"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.position = "bottom"
    ) +
    guides(color = guide_legend(title = "Assembly"))
  
  # Add N50 lines for each assembly with matching colors
  for (assembly_name in assembly_names) {
    nx_data <- nx_data_list[[assembly_name]]
    n50_value <- nx_data$NX[nx_data$X == 50]
    assembly_color <- colors[assembly_name]  # Get color by name to ensure match
    
    p <- p + geom_hline(yintercept = n50_value, 
                       linetype = "dashed", color = assembly_color, alpha = 0.5)
  }
  
  ggsave(output_file, plot = p, width = 12, height = 8, dpi = 300)
  cat("Comparative NX plot saved to:", output_file, "\n")
}

# Function to calculate basic assembly statistics
calculate_basic_stats <- function(lengths) {
  stats <- list(
    num_contigs = length(lengths),
    total_length = sum(lengths),
    mean_length = mean(lengths),
    median_length = median(lengths),
    min_length = min(lengths),
    max_length = max(lengths),
    std_length = sd(lengths)
  )
  return(stats)
}

# Function to create comparative statistics table
create_comparative_table <- function(results_list, output_file) {
  # Extract key metrics for comparison
  comparison_data <- data.frame()
  
  for (assembly_name in names(results_list)) {
    result <- results_list[[assembly_name]]
    basic_stats <- result$contiguity_results$basic_stats
    nx_data <- result$contiguity_results$nx_data
    aun_value <- result$contiguity_results$aun
    issues_result <- result$contiguity_results$issues
    
    # Extract issues count (handle case where issues calculation failed)
    issues_count <- if (!is.null(issues_result) && !any(is.na(issues_result)) && !is.null(issues_result$total_issues)) {
      issues_result$total_issues
    } else {
      NA
    }
    
    row_data <- data.frame(
      Assembly = assembly_name,
      Num_Contigs = basic_stats$num_contigs,
      Total_Length = basic_stats$total_length,
      Mean_Length = round(basic_stats$mean_length),
      Median_Length = basic_stats$median_length,
      N10 = nx_data$NX[which.min(abs(nx_data$X - 10))],
      N50 = nx_data$NX[which.min(abs(nx_data$X - 50))],
      N90 = nx_data$NX[which.min(abs(nx_data$X - 90))],
      auN = round(aun_value),
      Max_Contig = basic_stats$max_length,
      Issues = issues_count
    )
    
    comparison_data <- rbind(comparison_data, row_data)
  }
  
  # Save comparison table
  write_csv(comparison_data, output_file)
  cat("Comparative statistics saved to:", output_file, "\n")
  
  # Print comparison table to console
  cat("\n=== ASSEMBLY COMPARISON TABLE ===\n")
  print(comparison_data)
  
  return(comparison_data)
}

# Function to create length distribution comparison plot
create_length_distribution_plot <- function(lengths_list, output_file) {
  # Prepare data for plotting
  combined_data <- data.frame()
  assembly_names <- names(lengths_list)
  colors <- rainbow(length(lengths_list))
  names(colors) <- assembly_names  # Name the colors vector to ensure proper mapping

  for (assembly_name in assembly_names) {
    lengths <- lengths_list[[assembly_name]]
    length_data <- data.frame(
      Length = lengths,
      Assembly = assembly_name
    )
    combined_data <- rbind(combined_data, length_data)
  }

  # Create log-scale histogram
  p <- ggplot(combined_data, aes(x = Length, fill = Assembly)) +
    geom_histogram(alpha = 0.7, bins = 50, position = "identity") +
    scale_x_log10(labels = comma_format()) +
    scale_y_continuous(labels = comma_format()) +
    scale_fill_manual(values = colors) +  # Add this line to use consistent colors
    labs(
      title = "Contig Length Distribution Comparison",
      x = "Contig Length (bp, log scale)",
      y = "Count",
      subtitle = "Distribution of contig lengths across assemblies"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.position = "bottom"
    ) +
    facet_wrap(~Assembly, scales = "free_y")
  
  ggsave(output_file, plot = p, width = 12, height = 8, dpi = 300)
  cat("Length distribution plot saved to:", output_file, "\n")
}

# CONTIGUITY METRICS
calculate_contiguity_metrics <- function(lengths, assembly_name, output_dir, coverage_file = NULL) {
  cat("\n=== CONTIGUITY METRICS FOR", assembly_name, "===\n")
  
  # Basic statistics
  basic_stats <- calculate_basic_stats(lengths)
  
  cat("Basic Assembly Statistics:\n")
  cat(sprintf("  Number of contigs: %s\n", comma(basic_stats$num_contigs)))
  cat(sprintf("  Total length: %s bp\n", comma(basic_stats$total_length)))
  cat(sprintf("  Mean length: %s bp\n", comma(round(basic_stats$mean_length))))
  cat(sprintf("  Median length: %s bp\n", comma(basic_stats$median_length)))
  cat(sprintf("  Min length: %s bp\n", comma(basic_stats$min_length)))
  cat(sprintf("  Max length: %s bp\n", comma(basic_stats$max_length)))
  cat(sprintf("  Std deviation: %s bp\n", comma(round(basic_stats$std_length))))
  
  # Calculate NX values
  nx_data <- calculate_nx(lengths)
  
  cat("\nNX Values:\n")
  step_size <- 1000
  for (i in 1:nrow(nx_data)) {
    if (i %% step_size == 0){
      cat(sprintf("  N%f: %s bp\n", nx_data$X[i], comma(nx_data$NX[i])))
    }
  }
  
  # Calculate auN
  aun_value <- calculate_aun(lengths)
  cat(sprintf("\nauN (Area under NX curve): %s bp\n", comma(round(aun_value))))
  
  # Calculate Issues metric if coverage file is provided
  issues_result <- calculate_issues_metric(coverage_file, eps = 7500, output_dir, assembly_name)
  if (!is.null(issues_result) && !any(is.na(issues_result)) && !is.null(issues_result$total_issues)) {
    cat(sprintf("\nIssues metric: %s\n", comma(issues_result$total_issues)))
  } else {
    cat("\nIssues metric: Not available (no coverage file or calculation failed)\n")
  }
  
  # Create NX plot
  nx_plot_file <- file.path(output_dir, paste0(assembly_name, "_nx_plot.png"))
  create_nx_plot(nx_data, assembly_name, nx_plot_file)
  
  # Save NX data to CSV
  nx_csv_file <- file.path(output_dir, paste0(assembly_name, "_nx_values.csv"))
  write_csv(nx_data, nx_csv_file)
  cat("NX values saved to:", nx_csv_file, "\n")
  
  return(list(
    basic_stats = basic_stats,
    nx_data = nx_data,
    aun = aun_value,
    issues = issues_result
  ))
}

# Function to calculate Issues metric from coverage data
# This function identifies problematic regions in genome assemblies based on coverage discontinuities
# 
# Output files generated:
# 1. {assembly_name}_issues_clusters_detailed.csv - Complete list of all positions in issue clusters
#    Columns: chromosome, cluster_id, position, coverage, coverage_diff, cluster_size, cluster_span
# 2. {assembly_name}_issues_clusters_summary.csv - Summary statistics for each cluster
#    Columns: chromosome, cluster_id, cluster_size, cluster_span, start_position, end_position, 
#             min_coverage, max_coverage, mean_coverage, max_coverage_diff
#
# Users can use these files to:
# - Identify exact genomic coordinates of problematic regions
# - Prioritize clusters by size or coverage characteristics
# - Manual inspection and potential correction of assembly issues
calculate_issues_metric <- function(coverage_file, eps = 7500, output_dir = ".", assembly_name = "assembly") {
  if (is.null(coverage_file) || !file.exists(coverage_file)) {
    cat("Coverage file not provided or does not exist. Skipping Issues metric.\n")
    return(NA)
  }
  
  tryCatch({
    cat("Calculating Issues metric from coverage file:", coverage_file, "\n")
    
    # Read coverage data
    coverage_data <- read_tsv(coverage_file, col_names = c("chr_name", "position", "coverage"), 
                             show_col_types = FALSE)
    
    cat("Read", nrow(coverage_data), "coverage records\n")
    cat("Unique chromosomes:", length(unique(coverage_data$chr_name)), "\n")
    
    # Process each chromosome separately using traditional approach
    chromosomes <- unique(coverage_data$chr_name)
    issues_by_chr <- data.frame(
      chr_name = character(0),
      issues_count = numeric(0),
      stringsAsFactors = FALSE
    )
    
    # Store detailed cluster information for output file
    all_clusters <- data.frame(
      chromosome = character(0),
      cluster_id = character(0),
      position = numeric(0),
      coverage = numeric(0),
      coverage_diff = numeric(0),
      cluster_size = numeric(0),
      cluster_span = numeric(0),
      stringsAsFactors = FALSE
    )
    
    for (chr in chromosomes) {
      # Filter data for current chromosome
      chr_data <- coverage_data[coverage_data$chr_name == chr, ]
      
      cat("Processing chromosome:", chr, "- positions:", nrow(chr_data), "\n")
      
      # Calculate coverage differences
      coverage_diff <- c(0, diff(chr_data$coverage))
      
      # Calculate percentiles for outlier detection
      percentiles <- quantile(coverage_diff, probs = c(0.000001, 0.999999), na.rm = TRUE)
      
      cat("  Coverage diff range:", min(coverage_diff, na.rm = TRUE), "to", max(coverage_diff, na.rm = TRUE), "\n")
      cat("  Percentiles - 0.0001%:", percentiles[1], ", 99.9999%:", percentiles[2], "\n")
      
      # Find positions with extreme coverage differences
      outlier_mask <- coverage_diff > percentiles[2] | coverage_diff < percentiles[1]
      outlier_positions <- chr_data$position[outlier_mask]
      
      cat("  Found", length(outlier_positions), "outlier positions\n")
      
      if (length(outlier_positions) < 2) {
        cat("  Not enough outliers for clustering\n")
        issues_count <- 0
      } else {
        # Apply DBSCAN clustering
        clusters <- dbscan(matrix(outlier_positions, ncol = 1), eps = eps, minPts = 2)
        
        # Count non-outlier clusters (exclude noise points with cluster = 0)
        unique_clusters <- unique(clusters$cluster[clusters$cluster > 0])
        issues_count <- length(unique_clusters)
        
        cat("  DBSCAN found", length(unique(clusters$cluster)), "total clusters,", issues_count, "non-noise clusters\n")
        
        # Store detailed cluster information
        if (issues_count > 0) {
          for (cluster_label in unique_clusters) {
            cluster_positions <- outlier_positions[clusters$cluster == cluster_label]
            cluster_indices <- which(chr_data$position %in% cluster_positions)
            
            # Calculate cluster statistics
            cluster_span <- max(cluster_positions) - min(cluster_positions)
            cluster_size <- length(cluster_positions)
            
            # Store each position in this cluster
            for (idx in cluster_indices) {
              all_clusters <- rbind(all_clusters, data.frame(
                chromosome = chr,
                cluster_id = paste0(chr, "_cluster_", cluster_label),
                position = chr_data$position[idx],
                coverage = chr_data$coverage[idx],
                coverage_diff = coverage_diff[idx],
                cluster_size = cluster_size,
                cluster_span = cluster_span,
                stringsAsFactors = FALSE
              ))
            }
          }
        }
      }
      
      # Add to results
      issues_by_chr <- rbind(issues_by_chr, data.frame(
        chr_name = chr,
        issues_count = issues_count,
        stringsAsFactors = FALSE
      ))
    }
    
    # Sum up issues across all chromosomes
    total_issues <- sum(issues_by_chr$issues_count, na.rm = TRUE)
    
    cat("Issues metric calculation completed. Total issues:", total_issues, "\n")
    cat("Issues by chromosome:\n")
    print(issues_by_chr)
    
    # Save detailed cluster information to file
    if (nrow(all_clusters) > 0) {
      # Sort clusters by chromosome and position
      all_clusters <- all_clusters[order(all_clusters$chromosome, all_clusters$position), ]
      
      # Create cluster summary using base R
      cluster_ids <- unique(all_clusters$cluster_id)
      cluster_summary <- data.frame(
        chromosome = character(0),
        cluster_id = character(0),
        cluster_size = numeric(0),
        cluster_span = numeric(0),
        start_position = numeric(0),
        end_position = numeric(0),
        min_coverage = numeric(0),
        max_coverage = numeric(0),
        mean_coverage = numeric(0),
        max_coverage_diff = numeric(0),
        stringsAsFactors = FALSE
      )
      
      for (cid in cluster_ids) {
        cluster_data <- all_clusters[all_clusters$cluster_id == cid, ]
        summary_row <- data.frame(
          chromosome = cluster_data$chromosome[1],
          cluster_id = cid,
          cluster_size = cluster_data$cluster_size[1],
          cluster_span = cluster_data$cluster_span[1],
          start_position = min(cluster_data$position),
          end_position = max(cluster_data$position),
          min_coverage = min(cluster_data$coverage),
          max_coverage = max(cluster_data$coverage),
          mean_coverage = round(mean(cluster_data$coverage), 2),
          max_coverage_diff = max(abs(cluster_data$coverage_diff)),
          stringsAsFactors = FALSE
        )
        cluster_summary <- rbind(cluster_summary, summary_row)
      }
      
      # Save detailed clusters file
      clusters_file <- file.path(output_dir, paste0(assembly_name, "_issues_clusters_detailed.csv"))
      write_csv(all_clusters, clusters_file)
      cat("Detailed cluster information saved to:", clusters_file, "\n")
      
      # Save cluster summary file
      summary_file <- file.path(output_dir, paste0(assembly_name, "_issues_clusters_summary.csv"))
      write_csv(cluster_summary, summary_file)
      cat("Cluster summary saved to:", summary_file, "\n")
    } else {
      cat("No issue clusters found to save.\n")
    }
    
    return(list(
      total_issues = total_issues,
      issues_by_chr = issues_by_chr,
      cluster_details = all_clusters
    ))
    
  }, error = function(e) {
    cat("Error calculating Issues metric:", e$message, "\n")
    cat("Note: This function requires the 'dbscan' package. Install with: install.packages('dbscan')\n")
    cat("Alternative clustering methods available:\n")
    cat("  - k-means clustering: kmeans()\n")
    cat("  - hierarchical clustering: hclust()\n")
    cat("  - mixture models: mixtools package\n")
    return(NA)
  })
}

# COMPLETENESS METRICS (placeholder)
calculate_completeness_metrics <- function(lengths, assembly_name, output_dir, reference_file = NULL) {
  cat("\n=== COMPLETENESS METRICS ===\n")
  cat("Placeholder for completeness analysis\n")
  
  # TODO: Implement completeness metrics such as:
  # - BUSCO analysis for gene completeness
  # - Comparison with reference genome (if provided)
  # - Coverage analysis
  # - Gap analysis
  
  if (!is.null(reference_file) && file.exists(reference_file)) {
    cat("Reference genome provided:", reference_file, "\n")
    # TODO: Compare assembly to reference
    # - Calculate genome coverage
    # - Identify missing regions
    # - Calculate sequence identity
  } else {
    cat("No reference genome provided for completeness analysis\n")
  }
  
  cat("Completeness analysis functions to be implemented here\n")
  
  return(list(
    completeness_score = NA,
    coverage_stats = NA,
    missing_regions = NA
  ))
}

# CORRECTNESS METRICS (placeholder)
calculate_correctness_metrics <- function(lengths, assembly_name, output_dir, reference_file = NULL, reads_file = NULL) {
  cat("\n=== CORRECTNESS METRICS ===\n")
  cat("Placeholder for correctness analysis\n")
  
  # TODO: Implement correctness metrics such as:
  # - Misassembly detection
  # - Structural variant analysis
  # - Read mapping quality analysis
  # - Consensus accuracy
  
  if (!is.null(reference_file) && file.exists(reference_file)) {
    cat("Reference genome provided for correctness analysis:", reference_file, "\n")
    # TODO: Compare assembly to reference for correctness
    # - Detect inversions, translocations, duplications
    # - Calculate error rates
    # - Identify problematic regions
  }
  
  if (!is.null(reads_file) && file.exists(reads_file)) {
    cat("Reads file provided for correctness analysis:", reads_file, "\n")
    # TODO: Use reads for correctness assessment
    # - Map reads back to assembly
    # - Calculate mapping statistics
    # - Identify regions with poor read support
  }
  
  cat("Correctness analysis functions to be implemented here\n")
  
  return(list(
    error_rate = NA,
    misassemblies = NA,
    structural_variants = NA
  ))
}

# Main function
main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) < 1) {
    cat("Usage: Rscript assembly_stats.R <assembly1.fasta:name1[:coverage1.txt]> [assembly2.fasta:name2[:coverage2.txt]] ... [output_dir]\n")
    cat("\nArguments:\n")
    cat("  assemblyX.fasta:nameX[:coverageX.txt] - Path to assembly FASTA file with optional name and coverage file\n")
    cat("                                          Format: path:name or path:name:coverage_file\n")
    cat("  output_dir                            - Optional: Output directory (default: current directory)\n")
    cat("\nExamples:\n")
    cat("  # Single assembly\n")
    cat("  Rscript assembly_stats.R assembly.fasta:MyAssembly\n")
    cat("  # Multiple assemblies for comparison\n")
    cat("  Rscript assembly_stats.R original.fasta:Original improved.fasta:Improved output_dir\n")
    cat("  # With coverage file for Issues metric\n")
    cat("  Rscript assembly_stats.R assembly.fasta:MyAssembly:coverage.txt output_dir\n")
    cat("  # Mixed - some with coverage, some without\n")
    cat("  Rscript assembly_stats.R original.fasta:Original improved.fasta:Improved:coverage.txt output_dir\n")
    cat("  # Using default names\n")
    cat("  Rscript assembly_stats.R assembly1.fasta assembly2.fasta\n")
    quit(status = 1)
  }
  
  # Parse arguments
  output_dir <- "."
  assembly_inputs <- args
  
  # Check if last argument is output directory
  # An argument is considered an output directory if:
  # 1. It doesn't contain a colon (assembly:name format)
  # 2. AND it doesn't end with fasta/fa/fna extensions
  # 3. OR it exists as a directory
  last_arg <- args[length(args)]
  if (length(args) > 1 && 
      (!grepl(":", last_arg) && !grepl("\\.(fasta|fa|fna)$", last_arg)) ||
      (dir.exists(last_arg) && !grepl(":", last_arg))) {
    output_dir <- last_arg
    assembly_inputs <- args[-length(args)]
  }
  
  cat("Parsed arguments:\n")
  cat("  Assembly inputs:", paste(assembly_inputs, collapse = ", "), "\n")
  cat("  Output directory:", output_dir, "\n")
  cat("\n")
  
  # Parse assembly files, names, and coverage files
  assemblies <- list()
  assembly_coverage <- list()  # Store coverage file for each assembly
  
  for (input in assembly_inputs) {
    assembly_file <- NULL
    assembly_name <- NULL
    coverage_file <- NULL
    
    if (grepl(":", input)) {
      # Format: path:name or path:name:coverage
      parts <- strsplit(input, ":")[[1]]
      assembly_file <- parts[1]
      
      if (length(parts) >= 2) {
        assembly_name <- parts[2]
      } else {
        assembly_name <- tools::file_path_sans_ext(basename(assembly_file))
      }
      
      if (length(parts) >= 3) {
        coverage_file <- parts[3]
        # Check if coverage file exists
        if (!file.exists(coverage_file)) {
          cat("WARNING: Coverage file does not exist:", coverage_file, "\n")
          coverage_file <- NULL
        }
      }
    } else {
      # Format: just path, derive name from filename
      assembly_file <- input
      assembly_name <- tools::file_path_sans_ext(basename(assembly_file))
    }
    
    # Check if assembly file exists
    if (!file.exists(assembly_file)) {
      cat("WARNING: Assembly file does not exist:", assembly_file, "\n")
      next
    }
    
    assemblies[[assembly_name]] <- assembly_file
    assembly_coverage[[assembly_name]] <- coverage_file
    
    cat("Added assembly:", assembly_name, "->", assembly_file)
    if (!is.null(coverage_file)) {
      cat(" (coverage:", coverage_file, ")")
    }
    cat("\n")
  }
  
  if (length(assemblies) == 0) {
    cat("ERROR: No valid assembly files found\n")
    quit(status = 1)
  }
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  cat("Assembly Statistics Analysis\n")
  cat("============================\n")
  cat("Number of assemblies:", length(assemblies), "\n")
  cat("Output directory:", output_dir, "\n\n")
  
  # Process each assembly
  all_results <- list()
  all_lengths <- list()
  all_nx_data <- list()
  
  for (assembly_name in names(assemblies)) {
    assembly_file <- assemblies[[assembly_name]]
    coverage_file <- assembly_coverage[[assembly_name]]
    
    cat("Processing assembly:", assembly_name, "(", assembly_file, ")")
    if (!is.null(coverage_file)) {
      cat(" with coverage file:", coverage_file)
    }
    cat("\n")
    
    # Read assembly lengths
    lengths <- read_fasta_lengths(assembly_file)
    cat("Found", length(lengths), "sequences\n")
    
    # Store lengths for comparison
    all_lengths[[assembly_name]] <- lengths
    
    # Calculate metrics (pass assembly-specific coverage file)
    contiguity_results <- calculate_contiguity_metrics(lengths, assembly_name, output_dir, coverage_file)
    completeness_results <- calculate_completeness_metrics(lengths, assembly_name, output_dir, NULL)
    correctness_results <- calculate_correctness_metrics(lengths, assembly_name, output_dir, NULL, NULL)
    
    # Store results
    all_results[[assembly_name]] <- list(
      contiguity_results = contiguity_results,
      completeness_results = completeness_results,
      correctness_results = correctness_results
    )
    
    # Store NX data for comparison
    all_nx_data[[assembly_name]] <- contiguity_results$nx_data
    
    cat("\n" , rep("=", 50), "\n")
  }
  
  # Generate comparative analyses if multiple assemblies
  if (length(assemblies) > 1) {
    cat("\n=== GENERATING COMPARATIVE ANALYSES ===\n")
    
    # Create comparative NX plot
    comparative_nx_file <- file.path(output_dir, "comparative_nx_plot.png")
    create_comparative_nx_plot(all_nx_data, comparative_nx_file)
    
    # Create comparative statistics table
    comparative_table_file <- file.path(output_dir, "assembly_comparison.csv")
    comparison_data <- create_comparative_table(all_results, comparative_table_file)
    
    # Create length distribution comparison
    length_dist_file <- file.path(output_dir, "length_distribution_comparison.png")
    create_length_distribution_plot(all_lengths, length_dist_file)
    
    cat("\nComparative analysis complete!\n")
  }
  
  # Save comprehensive summary
  summary_file <- file.path(output_dir, "assembly_analysis_summary.txt")
  sink(summary_file)
  cat("Assembly Statistics Summary\n")
  cat("==========================\n")
  cat("Generated on:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("Number of assemblies analyzed:", length(assemblies), "\n\n")
  
  for (assembly_name in names(all_results)) {
    result <- all_results[[assembly_name]]
    basic_stats <- result$contiguity_results$basic_stats
    nx_data <- result$contiguity_results$nx_data
    aun_value <- result$contiguity_results$aun
    issues_result <- result$contiguity_results$issues
    
    # Extract issues count for display
    issues_count <- if (!is.null(issues_result) && !any(is.na(issues_result)) && !is.null(issues_result$total_issues)) {
      issues_result$total_issues
    } else {
      "N/A"
    }
    
    cat("ASSEMBLY:", assembly_name, "\n")
    cat("Number of contigs:", comma(basic_stats$num_contigs), "\n")
    cat("Total length:", comma(basic_stats$total_length), "bp\n")
    cat("N50:", comma(nx_data$NX[which.min(abs(nx_data$X - 50))]), "bp\n")
    cat("auN:", comma(round(aun_value)), "bp\n")
    cat("Issues:", issues_count, "\n\n")
  }
  
  if (length(assemblies) > 1) {
    cat("COMPARATIVE FILES GENERATED:\n")
    cat("- Comparative NX plot: comparative_nx_plot.png\n")
    cat("- Assembly comparison table: assembly_comparison.csv\n")
    cat("- Length distribution comparison: length_distribution_comparison.png\n")
  }
  
  sink()
  
  cat("\nAnalysis complete!\n")
  cat("Summary saved to:", summary_file, "\n")
  
  if (length(assemblies) > 1) {
    cat("\nComparative analyses generated for", length(assemblies), "assemblies\n")
    cat("Check the output directory for comparative plots and statistics\n")
  }
}

# Run main function
if (!interactive()) {
  main()
}
