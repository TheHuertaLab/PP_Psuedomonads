# Co-occurrence Analysis of Prophages and Defense Systems
# Author: Prasanna Joglekar
# Date: 2025-09-02

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
        tree <- read.newick("/Users/pjoglekar/work/pseudomonas/Pseudomonas_data_from_selected_genomes/pacbio_genomes/statistical_analysis/core_gene_alignment_filtered_align2.newick")
        data <- read.csv("/Users/pjoglekar/work/pseudomonas/Pseudomonas_data_from_selected_genomes/pacbio_genomes/statistical_analysis/Converged_presence_absence_2.txt", row.names = 1, sep = "\t", header = TRUE, check.names = FALSE)
    },
    warning = function(w) log_warning(w, "Load tree and data")
)

# Standardize genome names in data and tree (original)
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

# Filter rare systems: present in <0.5% of genomes (E. coli dataset)
min_prevalence <- 0.005 # Adjust to 0.01 for order-level datasets
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

# Check for near-invariant pairs
for (i in 1:(length(systems) - 1)) {
    for (j in (i + 1):length(systems)) {
        A <- data_td[, systems[i]]
        B <- data_td[, systems[j]]
        combined <- paste(A, B, sep = ":")
        state_counts <- table(combined)
        if (length(state_counts) < 4) {
            cat(sprintf("Pair %s vs %s has only %d states: %s\n", systems[i], systems[j], length(state_counts), paste(names(state_counts), collapse = ", ")))
        }
    }
}

# Parallelized pairwise comparisons
rds_file <- "cooccurrence_results.rds"
progress_log <- "progress.log"
if (file.exists(progress_log)) file.remove(progress_log)

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
                        # Use data_td with 1/2 coding (required by fitPagel and fitDiscrete)
                        A <- data_td[, systems[i]]
                        B <- data_td[, systems[j]]

                        fit_indep <- fitDiscrete(tree_td, cbind(A, B), model = "ARD", control = list(niter = 1000, fail = 1e-8))
                        fit_dep <- fitPagel(tree_td, A, B, model = "ARD", control = list(niter = 1000, fail = 1e-8))
                        fit_A_on_B <- fitPagel(tree_td, A, B, dep.var = "x", control = list(niter = 1000, fail = 1e-8))
                        fit_B_on_A <- fitPagel(tree_td, A, B, dep.var = "y", control = list(niter = 1000, fail = 1e-8))

                        p_indep_vs_dep <- NA
                        if (!is.null(fit_indep$opt$lnL) && !is.null(fit_dep$opt$lnL)) {
                            lr_stat <- 2 * (fit_dep$opt$lnL - fit_indep$opt$lnL)
                            p_indep_vs_dep <- pchisq(lr_stat, df = 4, lower.tail = FALSE)
                        }

                        best_model <- NA
                        flux_cooccur <- NA
                        flux_neg <- NA
                        association <- NA

                        if (!is.na(p_indep_vs_dep) && p_indep_vs_dep < 0.01) {
                            aic_values <- c(fit_dep$opt$aic, fit_A_on_B$opt$aic, fit_B_on_A$opt$aic)
                            best_model_idx <- which.min(aic_values)
                            best_model <- c("dependent", "A_on_B", "B_on_A")[best_model_idx]
                            best_fit <- list(fit_dep, fit_A_on_B, fit_B_on_A)[[best_model_idx]]

                            if (!is.null(best_fit$fit) && !is.null(best_fit$fit$rates)) {
                                rates <- best_fit$fit$rates
                                required_rates <- c("q12", "q21", "q13", "q24", "q31", "q42", "q43", "q14")
                                missing_rates <- setdiff(required_rates, names(rates))
                                if (length(missing_rates) > 0) {
                                    cat(sprintf("[%s] Missing rates: %s\n", step, paste(missing_rates, collapse = ", ")), file = warning_log, append = TRUE)
                                }
                                # Adjusted for fitPagel rate naming (0/1 coding)
                                q01 <- ifelse(exists("q13", rates), rates["q13"], NA) # (0,0) -> (0,1)
                                q10 <- ifelse(exists("q12", rates), rates["q12"], NA) # (0,0) -> (1,0)
                                q11_from01 <- ifelse(exists("q34", rates), rates["q34"], NA) # (0,1) -> (1,1)
                                q11_from10 <- ifelse(exists("q24", rates), rates["q24"], NA) # (1,0) -> (1,1)
                                q00_from01 <- ifelse(exists("q31", rates), rates["q31"], NA) # (0,1) -> (0,0)
                                q00_from10 <- ifelse(exists("q21", rates), rates["q21"], NA) # (1,0) -> (0,0)
                                q01_from11 <- ifelse(exists("q43", rates), rates["q43"], NA) # (1,1) -> (0,1)
                                q10_from11 <- ifelse(exists("q42", rates), rates["q42"], NA) # (1,1) -> (1,0)

                                if (all(!is.na(c(q11_from01, q01_from11, q11_from10, q10_from11, q01, q00_from01, q10, q00_from10)))) {
                                    flux_cooccur <- (q11_from01 / q01_from11) + (q11_from10 / q10_from11)
                                    flux_neg <- (q01 / q00_from01) + (q10 / q00_from10)
                                    association <- ifelse(flux_cooccur > flux_neg, "cooccur", "negative")
                                }
                            }
                        }

                        data.frame(
                            systemA = systems[i],
                            systemB = systems[j],
                            p_indep_vs_dep = p_indep_vs_dep,
                            best_model = best_model,
                            flux_cooccur = flux_cooccur,
                            flux_neg = flux_neg,
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
                        best_model = NA,
                        flux_cooccur = NA,
                        flux_neg = NA,
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

# Export data
if (exists("results") && nrow(results) > 0) {
    results$p_adj_bonf <- p.adjust(results$p_indep_vs_dep, method = "bonferroni")
    results$p_adj_bh <- p.adjust(results$p_indep_vs_dep, method = "BH")
    write.csv(results, "cooccurrence_results.csv", row.names = FALSE)
    cat("Results saved to cooccurrence_results.csv\n")
} else {
    cat("No significant pairs found.\n")
}
