# Co-occurrence Analysis of Prophages and Defense Systems
# Author: Prasanna Joglekar
# Date: 2025-08-19

# Load libraries
library(phytools)
library(ape)
library(geiger)
library(tidyverse)
library(doParallel)
library(foreach)

# Load tree and data
tree <- read.newick("/Users/pjoglekar/work/pseudomonas/Pseudomonas_data_from_selected_genomes/pacbio_genomes/statistical_analysis/core_gene_alignment_filtered_align2.newick")
data <- read.csv("/Users/pjoglekar/work/pseudomonas/Pseudomonas_data_from_selected_genomes/pacbio_genomes/statistical_analysis/Converged_presence_absence_2.txt", row.names = 1, sep = "\t", header = TRUE, check.names = FALSE)

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
tree <- multi2di(tree) # remove polytomies
tree <- collapse.singles(tree) # Collapse unbranched nodes

# Filter rare systems
min_prevalence <- 0.005 # 0.5%; adjust to 0.01 for 1%
prevalence <- colMeans(data)
data_filtered <- data[, prevalence >= min_prevalence, drop = FALSE] # Only filter columns, keep all rows
if (ncol(data_filtered) == 0) {
    stop("No systems remain after filtering. Lower the threshold.")
}
cat("Remaining systems:", ncol(data_filtered), "\n")

# Prepare traits (binary 0/1 for each pair of systems)
systems <- colnames(data_filtered)
td <- treedata(tree, data_filtered)
tree_td <- td$phy
data_td <- td$data

results <- data.frame(systemA = character(), systemB = character(), p_indep_vs_dep = numeric(), best_model = character(), flux_cooccur = numeric(), flux_neg = numeric(), p_adj_bonf = numeric(), p_adj_bh = numeric())
p_values <- numeric()

data_td_12 <- data_td + 1

total <- (length(systems) - 1) * length(systems) / 2
counter <- 0

cat("Starting co-occurrence analysis for", total, "system pairs...\n")

for (i in 1:(length(systems) - 1)) {
    for (j in (i + 1):length(systems)) {
        counter <- counter + 1
        cat(sprintf("Processing pair %d/%d: %s vs %s\n", counter, total, systems[i], systems[j]))

        A <- data_td_12[, systems[i]]
        B <- data_td_12[, systems[j]]

        # Fit models using treedata-matched tree and data
        fit_indep <- fitDiscrete(tree_td, cbind(A, B), model = "ARD", ncores = 16)
        fit_dep <- fitPagel(tree_td, A, B, model = "ARD")
        fit_A_on_B <- fitPagel(tree_td, A, B, dep.var = "x")
        fit_B_on_A <- fitPagel(tree_td, A, B, dep.var = "y")

        # Compare independent vs. dependent using AIC difference
        p_indep_vs_dep <- NA
        if (!is.null(fit_indep$AIC) && !is.null(fit_dep$AIC)) {
            p_indep_vs_dep <- fit_dep$AIC - fit_indep$AIC
        }
        p_values <- c(p_values, p_indep_vs_dep)

        if (!is.na(p_indep_vs_dep) && p_indep_vs_dep < -2) {
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

                results <- rbind(results, data.frame(systemA = systems[i], systemB = systems[j], p_indep_vs_dep = p_indep_vs_dep, best_model = best_model, flux_cooccur = flux_cooccur, flux_neg = flux_neg, p_adj_bonf = NA, p_adj_bh = NA))
            } else {
                cat(sprintf("Rates not found for %s vs %s - skipping.\n", systems[i], systems[j]))
            }
        }

        # Print a heartbeat every 10 pairs
        if (counter %% 10 == 0) {
            cat(sprintf("Heartbeat: %d pairs processed...\n", counter))
        }
    }
}

cat("Analysis complete. Processed", counter, "pairs.\n")

# Export data
if (length(p_values) > 0 && nrow(results) > 0) {
    results$p_adj_bonf <- p.adjust(p_values, method = "bonferroni")
    results$p_adj_bh <- p.adjust(p_values, method = "BH")
    write.csv(results, "cooccurrence_results.csv", row.names = FALSE)
    cat("Results saved to cooccurrence_results.csv\n")
} else {
    cat("No significant pairs found.\n")
}
