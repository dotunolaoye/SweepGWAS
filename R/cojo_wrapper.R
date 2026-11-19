# Resolves independent signals from GWAS summary statistics: builds the COJO
# input file from the user's columns, runs GCTA-COJO, applies the signal
# rules, and falls back to LD clumping when the joint model is unstable.
#
# The allele check (R/allele_check.R, which must be sourced first) runs
# before the input file is built, so that COJO never sees misoriented
# alleles or a frequency belonging to the other allele.
#
# plink2 and gcta64 are always called from inside work_dir with bare file
# names rather than full paths. plink2 on Windows cannot read WSL paths, and
# gcta64 splits a --bfile path at any space, even when the path is quoted,
# which breaks on a folder name that contains a space.
library(data.table)

run_plink2_relative = function(plink2_path, args, work_dir) {
  old_wd = getwd()
  on.exit(setwd(old_wd))
  setwd(work_dir)
  res = system2(plink2_path, args, stdout = TRUE, stderr = TRUE)
  res
}

#' Run COJO on user-supplied GWAS summary statistics and genotype data.
#'
#' @param gwas_raw data.frame or data.table of the user's GWAS output, with
#'   any column names.
#' @param column_map named character vector mapping standard names to the
#'   user's column names, e.g. c(SNP = "rsid", CHR = "chrom", POS = "bp",
#'   A1 = "effect_allele", A2 = "other_allele", FREQ = "eaf", BETA = "beta",
#'   P = "pval", N = "n_obs"). SE is optional; without it the standard error
#'   is calculated from BETA and P.
#' @param bfile_prefix path prefix (no extension) of the genotypes in PLINK
#'   bed/bim/fam format, which must already be in work_dir under this name.
#'   VCF input is converted first with convert_vcf_to_bed().
#' @param plink2_path,gcta_path paths to the plink2 and gcta64 programs;
#'   plink2 computes allele frequencies and, when needed, clumps.
#' @param work_dir a writable directory where both programs run.
#' @param cojo_p COJO significance threshold. NULL (default) uses a
#'   genome-wide Bonferroni threshold, 0.05 divided by the number of markers
#'   that pass the allele check.
#' @param out_prefix output file prefix (no extension).
#' @param few_signals_threshold below this number of signals the result
#'   carries a low-power warning.
#' @param same_source passed to the allele check; see
#'   verify_and_orient_alleles().
#' @param max_median_ratio,clump_p2,clump_r2,clump_kb see apply_cojo_rules()
#'   and clump_signals().
#'
#' @return a list: $signals (the independent signals the convergence test
#'   uses, either the COJO signals kept by apply_cojo_rules() or LD-clumped
#'   loci), $n_signals, $status ("ok", "no_signals" or "few_signals"),
#'   $message, $method ("cojo" or "clumping"), $cojo_output (GCTA's full
#'   output with the ratio and retained columns added), $rules,
#'   $clump_settings, $allele_check, $threshold_used, and $timing (seconds
#'   per stage).
run_cojo = function(gwas_raw, column_map, bfile_prefix, plink2_path, gcta_path,
                     work_dir, cojo_p = NULL, out_prefix = "cojo_result",
                     few_signals_threshold = 10, same_source = FALSE,
                     max_median_ratio = 2, clump_p2 = 1e-4, clump_r2 = 0.1, clump_kb = 250) {
  # seconds per stage, for the app's timing summary
  timing = c(allele_check = NA_real_, cojo = NA_real_, clumping = NA_real_)
  t_stage = Sys.time()
  lap = function() { el = as.numeric(difftime(Sys.time(), t_stage, units = "secs")); t_stage <<- Sys.time(); el }

  gwas_raw = as.data.table(gwas_raw)
  std = data.table(
    SNP = gwas_raw[[column_map["SNP"]]],
    CHR = gwas_raw[[column_map["CHR"]]],
    POS = gwas_raw[[column_map["POS"]]],
    A1 = gwas_raw[[column_map["A1"]]],
    A2 = gwas_raw[[column_map["A2"]]],
    FREQ = as.numeric(gwas_raw[[column_map["FREQ"]]]),
    BETA = as.numeric(gwas_raw[[column_map["BETA"]]]),
    P = as.numeric(gwas_raw[[column_map["P"]]]),
    N = as.numeric(gwas_raw[[column_map["N"]]])
  )
  if ("SE" %in% names(column_map) && !is.na(column_map["SE"])) {
    std[, SE := as.numeric(gwas_raw[[column_map["SE"]]])]
  } else {
    z = qnorm(std$P / 2, lower.tail = FALSE)
    std[, SE := abs(BETA) / z]
  }
  std = std[is.finite(SE) & SE > 0 & !is.na(P) & P > 0]

  cat("Computing genotype allele frequencies (plink2 --freq)...\n")
  freq_out = run_plink2_relative(plink2_path,
    c("--bfile", basename(bfile_prefix), "--freq", "--out", paste0(out_prefix, "_freq"), "--allow-extra-chr"),
    work_dir)
  freq_path = file.path(work_dir, paste0(out_prefix, "_freq.afreq"))
  if (!file.exists(freq_path)) stop("plink2 --freq failed:\n", paste(freq_out, collapse = "\n"))

  cat("Running mandatory allele-orientation check...\n")
  check = verify_and_orient_alleles(std, bim_path = paste0(bfile_prefix, ".bim"), freq_path = freq_path,
                                     same_source = same_source)
  print(check$summary)
  if (check$summary$match_rate < 0.20) {
    stop("Only ", round(check$summary$match_rate * 100, 1), "% of markers passed the allele-orientation ",
         "check (minimum 20% required). This usually means the genotype data does not correspond to ",
         "this GWAS summary statistics file. Check that the two files use the same marker IDs and the same ",
         "genome build, and that they come from the same dataset.")
  }

  timing["allele_check"] = lap()
  n_markers = nrow(check$matched)
  if (is.null(cojo_p)) cojo_p = 0.05 / n_markers
  cat("Building .ma file (", n_markers, "markers) and running COJO at p <", format(cojo_p, scientific = TRUE), "...\n")

  ma_path = file.path(work_dir, paste0(out_prefix, ".ma"))
  # GCTA expects exactly these eight columns in this order; the extra
  # freq_corrected column in check$matched would break its parser.
  fwrite(check$matched[, .(SNP, A1, A2, freq = FREQ, b = BETA, se = SE, p = P, N)],
         ma_path, sep = "\t", quote = FALSE)

  old_wd = getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(work_dir)
  cojo_args = c("--bfile", basename(bfile_prefix), "--cojo-file", paste0(out_prefix, ".ma"),
                "--cojo-slct", "--cojo-p", format(cojo_p, scientific = TRUE),
                "--out", out_prefix, "--thread-num", "3")
  cojo_res = system2(gcta_path, cojo_args, stdout = TRUE, stderr = TRUE)
  cojo_status = attr(cojo_res, "status")
  setwd(old_wd)
  if (!is.null(cojo_status) && cojo_status != 0) {
    stop("GCTA COJO failed (status ", cojo_status, "). Command: ", gcta_path, " ",
         paste(cojo_args, collapse = " "), "\nOutput:\n", paste(cojo_res, collapse = "\n"))
  }

  timing["cojo"] = lap()
  jma_path = file.path(work_dir, paste0(out_prefix, ".jma.cojo"))
  if (!file.exists(jma_path)) {
    return(list(signals = data.table(), n_signals = 0L, status = "no_signals", method = "cojo",
                message = "No independent signals passed the significance threshold. There is nothing to test. A less stringent threshold or a larger sample may help.",
                allele_check = check, threshold_used = cojo_p, timing = timing))
  }

  cojo_output = fread(jma_path)
  rules = apply_cojo_rules(cojo_output, cojo_p, max_median_ratio = max_median_ratio)
  cojo_output = rules$flagged

  if (rules$use_clumping) {
    cat("COJO joint model", if (rules$n_retained == 0) "retained no signal" else "unstable",
        "; taking independent loci by LD clumping instead...\n")
    signals = clump_signals(out_prefix, bfile_prefix, plink2_path, work_dir, cojo_p,
                            clump_p2 = clump_p2, clump_r2 = clump_r2, clump_kb = clump_kb)
    timing["clumping"] = lap()
    method = "clumping"
  } else {
    signals = cojo_output[retained == TRUE]
    method = "cojo"
  }

  n_signals = nrow(signals)
  status = if (n_signals == 0) "no_signals" else if (n_signals < few_signals_threshold) "few_signals" else "ok"
  message = switch(status,
    no_signals = "No independent signals passed the significance threshold. There is nothing to test. A less stringent threshold or a larger sample may help.",
    few_signals = sprintf("Only %d independent signal(s) were found. With fewer than %d signals, the test has little power to detect convergence. Interpret the result with caution.", n_signals, few_signals_threshold),
    ok = NA_character_
  )

  list(signals = signals, n_signals = n_signals, status = status, message = message,
       method = method, cojo_output = cojo_output, rules = rules[names(rules) != "flagged"],
       clump_settings = c(p2 = clump_p2, r2 = clump_r2, kb = clump_kb),
       allele_check = check, threshold_used = cojo_p, timing = timing)
}

#' Decide which COJO signals count as independent, and whether to trust the
#' joint model at all.
#'
#' Two rules are applied in order.
#'
#' 1. Marginal and joint (Rowan et al. 2024): a signal is kept only if it is
#'    significant both in the joint model (pJ) and on its own (p). A signal
#'    significant only jointly is usually one of a pair of correlated markers
#'    given large offsetting effects rather than a separate association. On
#'    the maize example 22 of 46 COJO signals were joint-only, with joint
#'    effects a median 9 times their marginal effects.
#'
#' 2. Stability: the joint model can be unstable even for signals that pass
#'    rule 1. On the Arabidopsis clim-bio6 GWAS (selfing, long-range LD and a
#'    trait that follows ancestry) the 19 kept signals still had joint
#'    effects a median 6.69 times their marginal effects and fell into only
#'    16 distinct LD clumps, whereas the 24 kept maize signals had a median
#'    ratio of 1.00, with one outlier at 10.1. The check therefore uses the
#'    median of |bJ|/|b| over the kept signals, which one outlier cannot
#'    trigger, and the default limit of 2 lies between the two values. When
#'    it is exceeded, or when rule 1 keeps nothing, independent loci are taken
#'    by LD clumping, which fits no joint model and so cannot inflate effects.
#'
#' @param cojo_output parsed .jma.cojo (columns p, pJ, b, bJ).
#' @param cojo_p the significance threshold COJO was run at.
#' @param max_median_ratio stability cut on the median |bJ|/|b|.
#' @return list: flagged (cojo_output plus ratio and retained columns),
#'   n_cojo, n_retained, n_joint_only, median_ratio, max_median_ratio,
#'   stable, use_clumping.
apply_cojo_rules = function(cojo_output, cojo_p, max_median_ratio = 2) {
  d = copy(as.data.table(cojo_output))
  d[, joint_marginal_ratio := abs(bJ) / abs(b)]
  d[, retained := p < cojo_p & pJ < cojo_p]
  n_retained = sum(d$retained)
  med = if (n_retained > 0) median(d[retained == TRUE]$joint_marginal_ratio) else NA_real_
  stable = n_retained > 0 && med <= max_median_ratio
  list(flagged = d, n_cojo = nrow(d), n_retained = n_retained,
       n_joint_only = sum(!d$retained), median_ratio = med,
       max_median_ratio = max_median_ratio, stable = stable,
       use_clumping = !stable)
}

#' Independent loci by LD clumping (plink2 --clump), used when the COJO
#' joint model cannot be trusted.
#'
#' Reads the .ma file already written by run_cojo(), so the alleles have
#' passed the orientation check. Index markers must pass the COJO threshold,
#' and the defaults (r2 0.1 within 250 kb, secondary threshold 1e-4) are
#' those of the Arabidopsis case study, which gave 110 loci from 1,502
#' significant markers. The output has the same core columns as the COJO
#' output (Chr, SNP, bp, p, and freq, b and se from the .ma file), so the
#' rest of the pipeline treats both alike.
clump_signals = function(out_prefix, bfile_prefix, plink2_path, work_dir, cojo_p,
                         clump_p2 = 1e-4, clump_r2 = 0.1, clump_kb = 250) {
  clump_out = paste0(out_prefix, "_clump")
  res = run_plink2_relative(plink2_path,
    c("--bfile", basename(bfile_prefix), "--clump", paste0(out_prefix, ".ma"),
      "--clump-id-field", "SNP", "--clump-p-field", "p",
      "--clump-p1", format(cojo_p, scientific = TRUE), "--clump-p2", format(clump_p2, scientific = TRUE),
      "--clump-r2", clump_r2, "--clump-kb", clump_kb,
      "--out", clump_out, "--allow-extra-chr"),
    work_dir)
  clumps_path = file.path(work_dir, paste0(clump_out, ".clumps"))
  if (!file.exists(clumps_path)) {
    if (any(grepl("No significant --clump results", res))) return(data.table())
    stop("plink2 --clump failed:\n", paste(res, collapse = "\n"))
  }
  cl = fread(clumps_path)
  setnames(cl, "#CHROM", "CHROM", skip_absent = TRUE)
  ma = fread(file.path(work_dir, paste0(out_prefix, ".ma")), select = c("SNP", "freq", "b", "se"))
  out = cl[, .(Chr = as.integer(CHROM), SNP = ID, bp = as.integer(POS), p = P,
               markers_in_clump = TOTAL + 1L)]
  out = merge(out, ma, by = "SNP", all.x = TRUE, sort = FALSE)
  setcolorder(out, c("Chr", "SNP", "bp", "freq", "b", "se", "p", "markers_in_clump"))
  setorder(out, Chr, bp)
  out
}
