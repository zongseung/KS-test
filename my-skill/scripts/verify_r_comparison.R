#!/usr/bin/env Rscript
#
# R Cross-Validation Script for Kruskal-Wallis Exact Distribution
#
# This script computes exact Kruskal-Wallis critical values and probabilities
# using R's coin package for comparison with Python implementation.
#
# Reference: Murakami, Lee & Ha (2009) - Higher Order Asymptotic
# Approximations of Kruskal-Wallis Statistics
#

# Install coin package if not available
if (!require("coin", quietly = TRUE)) {
    install.packages("coin", repos = "https://cloud.r-project.org/")
    library(coin)
}

# Function to compute H statistic from rank sums
compute_H <- function(rank_sums, sample_sizes) {
    N <- sum(sample_sizes)
    H <- 12 / (N * (N + 1)) * sum(rank_sums^2 / sample_sizes) - 3 * (N + 1)
    return(H)
}

# Function to enumerate all permutations (brute force for small N)
# This is the ground truth reference
enumerate_h_distribution <- function(sample_sizes) {
    N <- sum(sample_sizes)
    k <- length(sample_sizes)

    # All permutations of ranks
    ranks <- 1:N

    # For small N, enumerate all permutations
    library(gtools)
    perms <- permutations(N, N, ranks)

    h_values <- numeric(nrow(perms))

    for (i in 1:nrow(perms)) {
        perm <- perms[i,]
        idx <- 1
        rank_sums <- numeric(k)
        for (j in 1:k) {
            rank_sums[j] <- sum(perm[idx:(idx + sample_sizes[j] - 1)])
            idx <- idx + sample_sizes[j]
        }
        h_values[i] <- compute_H(rank_sums, sample_sizes)
    }

    return(h_values)
}

# Compute exact tail probability
exact_tail_prob <- function(h_values, h) {
    return(mean(h_values >= h))
}

# Compute exact critical value for given alpha
exact_critical_value <- function(h_values, alpha) {
    quantile(h_values, probs = 1 - alpha, type = 1)
}

# Main verification
cat("==========================================================\n")
cat("R VERIFICATION OF KRUSKAL-WALLIS EXACT DISTRIBUTION\n")
cat("==========================================================\n\n")

# Test case: (3,3,3)
cat("Configuration: (3,3,3)\n")
cat("--------------------------\n")

sample_sizes <- c(3, 3, 3)
N <- sum(sample_sizes)

# Load gtools for permutations
if (!require("gtools", quietly = TRUE)) {
    install.packages("gtools", repos = "https://cloud.r-project.org/")
    library(gtools)
}

# Enumerate all H values
h_values <- enumerate_h_distribution(sample_sizes)

# Unique H values and their frequencies
h_table <- table(round(h_values, 6))
cat("Distribution of H values:\n")
print(h_table / length(h_values))

# Critical value for alpha = 0.10
alpha <- 0.10
cv <- exact_critical_value(h_values, alpha)
ep <- exact_tail_prob(h_values, cv)

cat("\nFor alpha = 0.10:\n")
cat(sprintf("  Critical Value: %.6f\n", cv))
cat(sprintf("  P(H >= cv): %.6f\n", ep))

# Paper's CV
paper_cv <- 4.622222
paper_ep <- exact_tail_prob(h_values, paper_cv)
cat(sprintf("\nAt paper's CV (%.6f):\n", paper_cv))
cat(sprintf("  P(H >= %.6f) = %.6f\n", paper_cv, paper_ep))

# Using coin package for comparison (if available)
cat("\n==========================================================\n")
cat("USING COIN PACKAGE (kruskal_test)\n")
cat("==========================================================\n\n")

# Create sample data
set.seed(42)
n1 <- 3; n2 <- 3; n3 <- 3
data <- data.frame(
    value = c(rnorm(n1), rnorm(n2), rnorm(n3)),
    group = factor(rep(1:3, c(n1, n2, n3)))
)

# Exact test using coin
cat("Running exact Kruskal-Wallis test with coin package...\n")
result <- kruskal_test(value ~ group, data = data, distribution = "exact")
print(result)

# Monte Carlo approximation
cat("\nRunning Monte Carlo approximation (99999 replicates)...\n")
result_mc <- kruskal_test(value ~ group, data = data,
                          distribution = approximate(nresample = 99999))
print(result_mc)

cat("\n==========================================================\n")
cat("OUTPUT FOR PYTHON COMPARISON\n")
cat("==========================================================\n")

# Output in format easy to parse by Python
output <- list(
    config = "3,3,3",
    N = N,
    unique_h_values = as.numeric(names(h_table)),
    probabilities = as.numeric(h_table / length(h_values)),
    cv_alpha_0.10 = cv,
    ep_at_cv = ep,
    paper_cv = paper_cv,
    ep_at_paper_cv = paper_ep
)

# Save to JSON
if (require("jsonlite", quietly = TRUE)) {
    json_output <- toJSON(output, auto_unbox = TRUE, pretty = TRUE)
    cat("\nJSON Output:\n")
    cat(json_output)
    cat("\n")

    # Save to file
    write(json_output, file = "r_verification_output.json")
    cat("\nSaved to r_verification_output.json\n")
}

cat("\n==========================================================\n")
cat("ADDITIONAL CASES\n")
cat("==========================================================\n")

# Test other cases
test_cases <- list(
    c(2, 2, 2),
    c(4, 4, 4),
    c(2, 2, 2, 2),
    c(3, 3, 3, 3)
)

for (sample_sizes in test_cases) {
    config <- paste(sample_sizes, collapse = ",")
    N <- sum(sample_sizes)

    cat(sprintf("\nConfiguration: (%s), N=%d\n", config, N))

    if (N <= 10) {  # Only for small N due to computational limits
        h_values <- enumerate_h_distribution(sample_sizes)
        cv <- exact_critical_value(h_values, 0.10)
        ep <- exact_tail_prob(h_values, cv)
        cat(sprintf("  CV (alpha=0.10): %.6f\n", cv))
        cat(sprintf("  P(H >= cv): %.6f\n", ep))
    } else {
        cat("  Skipped (N too large for brute force)\n")
    }
}

cat("\n==========================================================\n")
cat("VERIFICATION COMPLETE\n")
cat("==========================================================\n")
