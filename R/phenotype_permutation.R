# Optional phenotype-permutation check, run locally only. It shuffles the
# phenotype labels, reruns the GWAS (GAPIT) on each shuffled version, and
# compares the observed number of significant markers with the resulting null
# distribution, which tests whether the association signal itself exceeds
# what an unrelated phenotype would produce. It assesses the raw GWAS signal,
# before COJO.
#
# It needs phenotype data, which the rest of the tool does not, and each
# replicate is a full GWAS, so it is not part of the app. Only the phenotype
# vector changes between replicates.
library(data.table)

#' Count markers reaching a significance threshold in one GAPIT run.
run_one_gwas_replicate = function(taxa, pheno_vec, GD, GM, model, pca_total,
                                    out_dir, sig_threshold) {
  y = data.frame(Taxa = taxa, Trait = pheno_vec)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  old_wd = getwd()
  on.exit(setwd(old_wd))
  setwd(out_dir)
  invisible(GAPIT(Y = y, GD = GD, GM = GM, model = model, PCA.total = pca_total, Model.selection = FALSE))
  result_file = list.files(out_dir, pattern = "GWAS_Results.*\\.csv$", full.names = TRUE)
  if (length(result_file) == 0) stop("GAPIT did not produce a GWAS_Results file in ", out_dir)
  r = fread(result_file[1])
  setnames(r, "P.value", "P", skip_absent = TRUE)
  sum(r$P < sig_threshold, na.rm = TRUE)
}

#' Run the phenotype-label permutation check.
#'
#' @param taxa character vector of sample IDs, in the same order as
#'   pheno_vec and matching GD's own taxa column.
#' @param pheno_vec numeric vector, the real (unpermuted) phenotype values.
#' @param GD,GM the GAPIT-format genotype dosage and map data.frames.
#' @param model,pca_total passed to GAPIT; the defaults (FarmCPU with three
#'   principal components) match the maize analysis.
#' @param sig_threshold NULL (default) uses a genome-wide Bonferroni threshold
#'   (0.05 divided by the number of markers in GD), as COJO does elsewhere in
#'   the tool.
#' @param n_perm number of permutation replicates. It has no default, because
#'   each replicate is a full GWAS and the caller should see the estimated
#'   runtime, printed after the observed run, before choosing a large value.
#' @param work_dir directory for GAPIT's per-replicate output subfolders.
#' @param seed for reproducibility.
#' @param progress_callback optional function(i, n_perm, elapsed) called
#'   after each replicate, for a UI progress indicator.
#'
#' @return a list: $observed_n_sig, $perm_counts (one count per replicate,
#'   for a histogram), $empirical_p, $sig_threshold, $n_perm and
#'   $observed_runtime_sec (the time taken by the unpermuted run).
run_phenotype_permutation = function(taxa, pheno_vec, GD, GM, model = "FarmCPU",
                                       pca_total = 3, sig_threshold = NULL,
                                       n_perm, work_dir, seed = 20260918,
                                       progress_callback = NULL) {
  if (is.null(sig_threshold)) sig_threshold = 0.05 / (ncol(GD) - 1)

  cat("Running observed (unpermuted) GWAS, for both the real result and a runtime estimate...\n")
  t0 = Sys.time()
  observed_n_sig = run_one_gwas_replicate(taxa, pheno_vec, GD, GM, model, pca_total,
                                            file.path(work_dir, "observed"), sig_threshold)
  observed_runtime_sec = as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat("Observed run took", round(observed_runtime_sec / 60, 1), "min;",
      "estimated total for", n_perm, "permutations: ~",
      round(observed_runtime_sec * n_perm / 60, 1), "min\n")

  set.seed(seed)
  perm_counts = integer(n_perm)
  for (i in seq_len(n_perm)) {
    shuffled = sample(pheno_vec)
    t_i = Sys.time()
    perm_counts[i] = run_one_gwas_replicate(taxa, shuffled, GD, GM, model, pca_total,
                                              file.path(work_dir, paste0("perm_", i)), sig_threshold)
    if (!is.null(progress_callback)) {
      progress_callback(i, n_perm, as.numeric(difftime(Sys.time(), t_i, units = "secs")))
    }
  }

  empirical_p = (sum(perm_counts >= observed_n_sig) + 1) / (n_perm + 1)

  list(observed_n_sig = observed_n_sig, perm_counts = perm_counts, empirical_p = empirical_p,
       sig_threshold = sig_threshold, n_perm = n_perm, observed_runtime_sec = observed_runtime_sec)
}
