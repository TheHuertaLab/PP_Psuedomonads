# Co-occurrence Analysis of Prophages and Defense Systems
# Author: Prasanna Joglekar
# Date: 2025-09-03
# Fixed version with proper fitPagel handling

# Load required libraries
library(phytools)
library(ape)
library(geiger)
library(foreach)
library(doParallel)

# Set up warning logging
warning_log <- "cooccurrence_warnings.log"
if (file.exists(warning_log)) file.remove(warning_log)
log_warning <- function(w, step) {
    msg <- paste0("[", step, "] ", conditionMessage(w), "\n")
    cat(msg, file = warning_log, append = TRUE)
}

# Load tree and data with warning handling
withCallingHandlers(
    {
        tree <- read.newick("core_gene_alignment_filtered_align2.newick")
        data <- read.csv("Converged_presence_absence_2.txt", row.names = 1, sep = "\t", header = TRUE, check.names = FALSE)
    },
    warning = function(w) log_warning(w, "Load tree and data")
)

# Standardize genome names in data and tree
rownames(data) <- gsub(" ", "_", rownames(data))
tree$tip.label <- gsub(" ", "_", tree$tip.label)
tree$tip.label <- gsub("^'|'$", "", tree$tip.label)

# Ensure common genomes
common_genomes <- intersect(rownames(data), tree$tip.label)
data <- data[common_genomes, , drop = FALSE]
tree <- drop.tip(tree, setdiff(tree$tip.label, common_genomes))

# Standardize tree branch lengths to mean 0.1
tree$edge.length <- tree$edge.length / mean(tree$edge.length) * 0.1
tree <- multi2di(tree)
tree <- collapse.singles(tree)

# Filter rare systems: present in <0.5% of genomes
min_prevalence <- 0.005
prevalence <- colMeans(data)
data_filtered <- data[, prevalence >= min_prevalence, drop = FALSE]
if (ncol(data_filtered) == 0) {
    stop("No systems remain after prevalence filtering. Lower the threshold.")
}
cat("Systems after prevalence filtering:", ncol(data_filtered), "\n")

# Prepare traits
systems <- colnames(data_filtered)
td <- treedata(tree, data_filtered)
tree_td <- td$phy
data_td <- td$data

# Fix order mismatch after treedata
if (!identical(rownames(data_td), tree_td$tip.label)) {
    cat("Order mismatch detected. Reordering data to match tree...\n")
    # Check if all genomes are present
    missing_in_tree <- setdiff(rownames(data_td), tree_td$tip.label)
    missing_in_data <- setdiff(tree_td$tip.label, rownames(data_td))

    if (length(missing_in_tree) > 0 || length(missing_in_data) > 0) {
        cat("Content mismatch after treedata:\n")
        cat("In data_td but not tree_td:", missing_in_tree, "\n")
        cat("In tree_td but not data_td:", missing_in_data, "\n")
        stop("Some genomes in data_td not found in tree_td after treedata!")
    }

    # Reorder data to match tree tip labels
    data_td <- data_td[tree_td$tip.label, , drop = FALSE]
}

# Final alignment check
if (!identical(rownames(data_td), tree_td$tip.label)) {
    cat("Final alignment check failed:\n")
    cat("Data rownames:", head(rownames(data_td), 10), "\n")
    cat("Tree tip labels:", head(tree_td$tip.label, 10), "\n")
    stop("Tree tip labels and data rownames do not match after reordering!")
}

# Convert traits to 1/2 coding for fitPagel and fitDiscrete (required by geiger package)
data_td <- data_td + 1

# Remove invariant systems (all 1 or all 2)
variable_systems <- systems[apply(data_td, 2, function(x) length(unique(x)) > 1)]
data_td <- data_td[, variable_systems, drop = FALSE]
systems <- variable_systems
cat("Variable systems after removing invariants:", length(systems), "\n")

# Parallelized pairwise comparisons
rds_file <- "cooccurrence_results.rds"
progress_log <- "progress.log"
warning_log <- "cooccurrence_warnings.log"
if (file.exists(progress_log)) file.remove(progress_log)
if (file.exists(warning_log)) file.remove(warning_log)

if (file.exists(rds_file)) {
    cat("Loading previous results from RDS file...\n")
    results <- readRDS(rds_file)
    cat("Loaded", nrow(results), "pairs from RDS.\n")
} else {
    total <- (length(systems) - 1) * length(systems) / 2
    cat("Starting co-occurrence analysis for", total, "system pairs...\n")

    n_cores <- parallel::detectCores() - 1
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)

    results <- foreach(i = 1:(length(systems) - 1), .combine = rbind, .packages = c("phytools", "ape", "geiger")) %:%
        foreach(j = (i + 1):length(systems), .combine = rbind) %dopar% {
            step <- sprintf("Co-occurrence analysis: %s vs %s", systems[i], systems[j])
            cat(sprintf("Processing pair: %s vs %s\n", systems[i], systems[j]), file = progress_log, append = TRUE)

            tryCatch(
                withCallingHandlers(
                    {
                        A <- data_td[, systems[i]]
                        B <- data_td[, systems[j]]

                        # Check if we have sufficient states for analysis
                        combined <- paste(A, B, sep = ":")
                        state_counts <- table(combined)

                        # Log state information
                        state_msg <- sprintf(
                            "[%s vs %s] Has %d states: %s",
                            systems[i], systems[j],
                            length(state_counts),
                            paste(names(state_counts), collapse = ", ")
                        )
                        cat(state_msg, "\n", file = warning_log, append = TRUE)

                        if (length(state_counts) < 3) {
                            # Skip pairs with too few states
                            return(data.frame(
                                systemA = systems[i],
                                systemB = systems[j],
                                p_indep_vs_dep = NA,
                                best_model = "insufficient_states",
                                association = NA,
                                p_adj_bonf = NA,
                                p_adj_bh = NA
                            ))
                        }

                        # Create multistate character for independent model
                        unique_states <- sort(unique(combined))
                        combined_states <- as.numeric(factor(combined, levels = unique_states))
                        names(combined_states) <- rownames(data_td)

                        # Fit independent model (multistate character)
                        fit_indep <- tryCatch(
                            fitDiscrete(tree_td, combined_states, model = "ER"),
                            error = function(e) NULL
                        )

                        # Fit dependent model (Pagel's method)
                        fit_dep <- tryCatch(
                            fitPagel(tree_td, A, B, model = "ARD"),
                            error = function(e) NULL
                        )

                        # Check if models fitted successfully
                        if (is.null(fit_indep) || is.null(fit_dep)) {
                            return(data.frame(
                                systemA = systems[i],
                                systemB = systems[j],
                                p_indep_vs_dep = NA,
                                best_model = "model_failed",
                                association = NA,
                                p_adj_bonf = NA,
                                p_adj_bh = NA
                            ))
                        }

                        # Extract p-value (fitPagel provides this directly)
                        p_indep_vs_dep <- NA
                        if (!is.null(fit_dep$P)) {
                            p_indep_vs_dep <- fit_dep$P
                        }

                        # Determine best model and association
                        best_model <- "independent"
                        association <- "not_significant"

                        if (!is.na(p_indep_vs_dep) && p_indep_vs_dep < 0.05) {
                            best_model <- "dependent"
                            association <- "significant"
                        }

                        data.frame(
                            systemA = systems[i],
                            systemB = systems[j],
                            p_indep_vs_dep = p_indep_vs_dep,
                            best_model = best_model,
                            association = association,
                            p_adj_bonf = NA,
                            p_adj_bh = NA
                        )
                    },
                    warning = function(w) log_warning(w, step)
                ),
                error = function(e) {
                    msg <- paste0("[", step, "] ERROR: ", conditionMessage(e), "\n")
                    cat(msg, file = warning_log, append = TRUE)
                    data.frame(
                        systemA = systems[i],
                        systemB = systems[j],
                        p_indep_vs_dep = NA,
                        best_model = "error",
                        association = NA,
                        p_adj_bonf = NA,
                        p_adj_bh = NA
                    )
                }
            )
        }

    stopCluster(cl)

    cat("Analysis complete. Processed", nrow(results), "pairs.\n")
    saveRDS(results, rds_file)
    cat("Results saved to", rds_file, "\n")
}

# Apply multiple testing corrections and export data
if (exists("results") && nrow(results) > 0) {
    # Only apply corrections to non-NA p-values
    valid_p <- !is.na(results$p_indep_vs_dep)
    if (sum(valid_p) > 0) {
        results$p_adj_bonf[valid_p] <- p.adjust(results$p_indep_vs_dep[valid_p], method = "bonferroni")
        results$p_adj_bh[valid_p] <- p.adjust(results$p_indep_vs_dep[valid_p], method = "BH")
    }

    write.csv(results, "cooccurrence_results.csv", row.names = FALSE)
    cat("Results saved to cooccurrence_results.csv\n")

    # Print summary
    cat("\nSummary of results:\n")
    cat("Total pairs analyzed:", nrow(results), "\n")
    cat("Successful analyses:", sum(!is.na(results$p_indep_vs_dep)), "\n")
    cat("Significant pairs (p < 0.05):", sum(results$p_indep_vs_dep < 0.05, na.rm = TRUE), "\n")
    cat("Significant pairs after Bonferroni correction:", sum(results$p_adj_bonf < 0.05, na.rm = TRUE), "\n")
    cat("Significant pairs after BH correction:", sum(results$p_adj_bh < 0.05, na.rm = TRUE), "\n")
} else {
    cat("No results to export.\n")
}
