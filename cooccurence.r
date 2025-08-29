# Co-occurrence Analysis of Prophages and Defense Systems
# Author: Prasanna Joglekar
# Date: 2025-08-19

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
    invokeRestart("muffleWarning")
}

# Load tree and data with warning handling
withCallingHandlers(
    {
        tree <- read.newick("/Users/pjoglekar/work/pseudomonas/Pseudomonas_data_from_selected_genomes/pacbio_genomes/statistical_analysis/core_gene_alignment_filtered_align2.newick")
        data <- read.csv("/Users/pjoglekar/work/pseudomonas/Pseudomonas_data_from_selected_genomes/pacbio_genomes/statistical_analysis/Converged_presence_absence_2.txt", row.names = 1, sep = "\t", header = TRUE, check.names = FALSE)
    },
    warning = function(w) log_warning(w, "Load tree and data")
)

# Standardize genome names in data and tree
rownames(data) <- gsub(" ", "_", rownames(data))
tree$tip.label <- gsub(" ", "_", tree$tip.label)
tree$tip.label <- gsub("^'|'$", "", tree$tip.label)

# Print names in data not in tree
not_in_tree <- setdiff(rownames(data), tree$tip.label)
if (length(not_in_tree) > 0) {
    cat("Genomes in data not in tree:\n")
    print(not_in_tree)
}

# Print names in tree not in data
not_in_data <- setdiff(tree$tip.label, rownames(data))
if (length(not_in_data) > 0) {
    cat("Genomes in tree not in data:\n")
    print(not_in_data)
}

# Standardize tree branch lengths to mean 0.1
tree$edge.length <- tree$edge.length / mean(tree$edge.length) * 0.1
tree <- multi2di(tree)
tree <- collapse.singles(tree)

# Filter rare systems
min_prevalence <- 0.005
prevalence <- colMeans(data)
data_filtered <- data[, prevalence >= min_prevalence, drop = FALSE]
if (ncol(data_filtered) == 0) {
    stop("No systems remain after filtering. Lower the threshold.")
}
cat("Remaining systems:", ncol(data_filtered), "\n")

# Prepare traits
systems <- colnames(data_filtered)
td <- treedata(tree, data_filtered)
tree_td <- td$phy
data_td <- td$data

# Parallelized pairwise comparisons
rds_file <- "cooccurrence_results.rds"

progress_log <- "progress.log"
if (file.exists(progress_log)) file.remove(progress_log)

if (file.exists(rds_file)) {
    cat("Loading previous results from RDS file...\n")
    results <- readRDS(rds_file)
    cat("Loaded", nrow(results), "pairs from RDS.\n")
} else {
    data_td_12 <- data_td + 1
    total <- (length(systems) - 1) * length(systems) / 2
    cat("Starting co-occurrence analysis for", total, "system pairs...\n")

    n_cores <- parallel::detectCores() - 1
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)

    results <- foreach(i = 1:(length(systems) - 1), .combine = rbind, .packages = c("phytools", "ape", "geiger")) %:%
        foreach(j = (i + 1):length(systems), .combine = rbind) %dopar% {
            step <- sprintf("Co-occurrence analysis: %s vs %s", systems[i], systems[j])
            # Log progress for each pair
            cat(sprintf("Processing pair: %s vs %s\n", systems[i], systems[j]), file = progress_log, append = TRUE)
            tryCatch(
                {
                    A <- data_td_12[, systems[i]]
                    B <- data_td_12[, systems[j]]

                    fit_indep <- fitDiscrete(tree_td, cbind(A, B), model = "ARD")
                    fit_dep <- fitPagel(tree_td, A, B, model = "ARD")
                    fit_A_on_B <- fitPagel(tree_td, A, B, dep.var = "x")
                    fit_B_on_A <- fitPagel(tree_td, A, B, dep.var = "y")

                    p_indep_vs_dep <- NA
                    if (!is.null(fit_indep$AIC) && !is.null(fit_dep$AIC)) {
                        p_indep_vs_dep <- fit_dep$AIC - fit_indep$AIC
                    }

                    best_model <- NA
                    flux_cooccur <- NA
                    flux_neg <- NA

                    if (!is.na(p_indep_vs_dep) && p_indep_vs_dep < 0) {
                        aic_values <- c(fit_dep$AIC, fit_A_on_B$AIC, fit_B_on_A$AIC)
                        best_model_idx <- which.min(aic_values)
                        best_model <- c("dependent", "A_on_B", "B_on_A")[best_model_idx]
                        best_fit <- list(fit_dep, fit_A_on_B, fit_B_on_A)[[best_model_idx]]

                        if (!is.null(best_fit$fit) && !is.null(best_fit$fit$rates)) {
                            rates <- best_fit$fit$ratestop
                            q01 <- rates["q13"]
                            q10 <- rates["q12"]
                            q11_from01 <- rates["q34"]
                            q11_from10 <- rates["q24"]
                            q00_from01 <- rates["q31"]
                            q00_from10 <- rates["q21"]
                            q01_from11 <- rates["q43"]
                            q10_from11 <- rates["q42"]

                            flux_cooccur <- (q11_from01 / q01_from11) + (q11_from10 / q10_from11)
                            flux_neg <- (q01 / q00_from01) + (q10 / q00_from10)
                        }
                    }

                    data.frame(
                        systemA = systems[i],
                        systemB = systems[j],
                        p_indep_vs_dep = p_indep_vs_dep,
                        best_model = best_model,
                        flux_cooccur = flux_cooccur,
                        flux_neg = flux_neg,
                        p_adj_bonf = NA,
                        p_adj_bh = NA
                    )
                },
                warning = function(w) {
                    msg <- paste0("[", step, "] ", conditionMessage(w), "\n")
                    cat(msg, file = warning_log, append = TRUE)
                    invokeRestart("muffleWarning")
                    data.frame(
                        systemA = systems[i],
                        systemB = systems[j],
                        p_indep_vs_dep = NA,
                        best_model = NA,
                        flux_cooccur = NA,
                        flux_neg = NA,
                        p_adj_bonf = NA,
                        p_adj_bh = NA
                    )
                },
                error = function(e) {
                    msg <- paste0("[", step, "] ERROR: ", conditionMessage(e), "\n")
                    cat(msg, file = warning_log, append = TRUE)
                    data.frame(
                        systemA = systems[i],
                        systemB = systems[j],
                        p_indep_vs_dep = NA,
                        best_model = NA,
                        flux_cooccur = NA,
                        flux_neg = NA,
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

# Export data
if (exists("results") && nrow(results) > 0) {
    results$p_adj_bonf <- p.adjust(results$p_indep_vs_dep, method = "bonferroni")
    results$p_adj_bh <- p.adjust(results$p_indep_vs_dep, method = "BH")
    write.csv(results, "cooccurrence_results.csv", row.names = FALSE)
    cat("Results saved to cooccurrence_results.csv\n")
} else {
    cat("No significant pairs found.\n")
}
