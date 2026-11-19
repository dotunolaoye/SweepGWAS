# Convergence test: counts the independent signals that lie within a
# distance threshold of a selection-scan region, and compares the count with a
# permutation null in which each signal is redrawn at a tested marker on its
# own chromosome. The distance threshold must come from the user, from the reported
# linkage-disequilibrium decay of the population, because an arbitrary wide
# distance threshold places most of the genome near a region and removes the test's
# power (at 500 kb, 85.3% of the maize genome).
library(data.table)

# Chromosome names are compared after dropping a "chr" prefix, surrounding
# spaces and letter case, so that 1, chr1 and Chr1 match across input files
# written by different programs.
std_chr = function(x) sub("^chr", "", trimws(as.character(x)), ignore.case = TRUE)

#' Build a window of the distance threshold on each side of every independent signal.
#'
#' @param cojo_signals data.table with at least Chr, SNP and bp, as returned
#'   by run_cojo()$signals.
#' @param tolerance half-width of the window in bp.
build_signal_regions = function(cojo_signals, tolerance) {
  data.table(
    chr = cojo_signals$Chr,
    lead_snp = cojo_signals$SNP,
    lead_pos = cojo_signals$bp,
    start = pmax(0, cojo_signals$bp - tolerance),
    end = cojo_signals$bp + tolerance
  )
}

#' Split the selection-scan regions by chromosome once, so the permutation
#' loop does not filter the whole table on every draw. With the 6,862 maize
#' regions and 1,000 permutations this cut the runtime from about 2.5 minutes
#' to 2 seconds.
split_selscan_by_chr = function(selscan_regions) {
  s = as.data.table(selscan_regions)
  split(s[, .(start, end)], as.character(s$chr))
}

#' Whether each position lies within `window` of at least one selection-scan
#' region on the same chromosome, from any method.
#'
#' @param selscan_by_chr the lookup from split_selscan_by_chr().
near_selscan = function(chrs, mids, selscan_by_chr, window) {
  chrs = as.character(chrs)
  vapply(seq_along(chrs), function(i) {
    sub = selscan_by_chr[[chrs[i]]]
    if (is.null(sub) || nrow(sub) == 0) return(FALSE)
    any(sub$start <= mids[i] + window & sub$end >= mids[i] - window)
  }, logical(1))
}

#' Run the convergence permutation test.
#'
#' @param cojo_signals data.table of independent signals (Chr, SNP, bp).
#' @param selscan_regions data.table with chr, start, end and an optional
#'   method. A signal counts as convergent if it lies within `tolerance` of
#'   any one region, from any method, rather than requiring every method to
#'   agree.
#' @param tolerance required, with no default: the reported
#'   linkage-disequilibrium decay distance for the population, for example
#'   10000 for the maize Ames panel within breeding groups (Romay et al.
#'   2013). It is not estimated automatically, because naive LD estimates are
#'   inflated by population structure; users without a published value can
#'   estimate one with PopLDdecay (Zhang et al. 2019).
#' @param marker_positions data.table with chr and pos for every marker
#'   tested in the GWAS, which are the only places a signal could have been
#'   detected and are what permuted positions are drawn from by default.
#' @param null how permuted positions are drawn. "markers" (default) draws
#'   from marker_positions. "sequence" draws uniformly along each chromosome,
#'   which misstates chance whenever markers and regions are unevenly
#'   distributed relative to each other. On the maize genotyping-by-sequencing
#'   data, where only 5 of 2,108 one-megabase bins hold a region without
#'   markers and bins with regions carry twice the marker density, the
#'   uniform null gives 10.0% chance convergence at 10 kb against 36.4% from
#'   tested markers, so the observed 37.5% appears significant (p = 0.000999)
#'   when it is not (p = 0.524). "sequence" is kept only for comparison and
#'   warns when used.
#' @param n_perm number of permutations. Each permutation only redraws
#'   positions, so 1,000 take seconds.
#' @param seed for reproducibility.
#'
#' @return a list: $regions (per signal: position, whether convergent, and
#'   distance to the nearest region, used by the plots), $observed,
#'   $n_regions, $observed_pct, $null_distribution (one count per
#'   permutation), $null_mean_pct, $empirical_p, $null_model and
#'   $tolerance_used.
run_convergence_test = function(cojo_signals, selscan_regions, tolerance,
                                 marker_positions, n_perm = 1000, seed = 20260918,
                                 null = c("markers", "sequence")) {
  null = match.arg(null)
  if (is.null(tolerance)) {
    stop("A distance threshold must be supplied: the linkage-disequilibrium decay distance reported for this ",
         "population. There is no default. An arbitrary distance threshold can leave the test unable to ",
         "detect anything.")
  }
  if (nrow(cojo_signals) == 0) {
    stop("There are no independent signals. The convergence test was not run.")
  }

  # Both lookups are built once and reused by every permutation.
  selscan_by_chr = split_selscan_by_chr(selscan_regions)
  mp = as.data.table(marker_positions)
  stopifnot(all(c("chr", "pos") %in% names(mp)))
  # The permutations draw markers and signals by row, so both are put in
  # genome order first; otherwise the same data in a different row order
  # would give a different p-value under the same seed.
  mp = mp[order(pos)]
  pos_by_chr = split(mp$pos, as.character(mp$chr))
  chr_lengths = mp[, .(max_pos = max(pos, na.rm = TRUE)), by = .(chr)]
  max_pos_by_chr = setNames(chr_lengths$max_pos, as.character(chr_lengths$chr))
  if (null == "sequence") {
    warning("null = \"sequence\" draws permuted positions uniformly along each chromosome, ",
            "including stretches where no marker was assayed and therefore where neither a GWAS ",
            "signal nor a selection-scan region could ever have been called. This can misstate chance ",
            "convergence in either direction when markers and regions are unevenly distributed. ",
            "Use null = \"markers\" unless reproducing a published analysis that used the uniform null.")
  }

  regions = build_signal_regions(cojo_signals, tolerance)
  regions = regions[order(suppressWarnings(as.numeric(as.character(chr))), as.character(chr), lead_pos)]
  regions[, mid := round((start + end) / 2)]
  regions[, convergent := near_selscan(chr, mid, selscan_by_chr, tolerance)]
  # distance to the nearest region, shown in the plots; 0 means inside one
  regions[, nearest_selscan_dist := vapply(seq_len(.N), function(i) {
    sub = selscan_by_chr[[as.character(regions$chr[i])]]
    if (is.null(sub) || nrow(sub) == 0) return(NA_real_)
    d = ifelse(regions$mid[i] >= sub$start & regions$mid[i] <= sub$end, 0,
               pmin(abs(regions$mid[i] - sub$start), abs(regions$mid[i] - sub$end)))
    min(d)
  }, numeric(1))]

  observed = sum(regions$convergent)
  n_regions = nrow(regions)

  set.seed(seed)
  rand_chrs = regions$chr
  rand_chr_chr = as.character(rand_chrs)
  chr_bounds = max_pos_by_chr[rand_chr_chr]
  draw = if (null == "markers") {
    # a signal can only arise where a marker was tested
    function() vapply(rand_chr_chr, function(cc) {
      v = pos_by_chr[[cc]]
      if (is.null(v) || length(v) == 0) return(NA_real_)
      v[sample.int(length(v), 1)]
    }, numeric(1))
  } else {
    function() vapply(chr_bounds, function(mx) {
      if (is.na(mx)) return(NA_real_)
      as.numeric(sample.int(mx, 1))
    }, numeric(1))
  }
  null_distribution = vapply(seq_len(n_perm), function(p) {
    sum(near_selscan(rand_chrs, draw(), selscan_by_chr, tolerance))
  }, integer(1))

  empirical_p = (sum(null_distribution >= observed) + 1) / (n_perm + 1)

  list(
    regions = regions,
    observed = observed,
    n_regions = n_regions,
    observed_pct = round(100 * observed / n_regions, 1),
    null_distribution = null_distribution,
    null_mean_pct = round(100 * mean(null_distribution) / n_regions, 1),
    empirical_p = empirical_p,
    null_model = null,
    tolerance_used = tolerance
  )
}

#' Run the convergence test under both null models and compare them.
#'
#' The reported result always uses the marker-based null. The uniform null is
#' run beside it as a check, because the difference between the two shows how unevenly
#' markers and regions are distributed, and because its error is silent: on the
#' maize example it gives p = 0.000999 at 10 kb against 0.524 for the
#' marker-based null, and nothing about its result looks wrong.
#'
#' The flag is raised when the two nulls lead to different conclusions at
#' `alpha`, with no other threshold, and the ratio of their chance baselines is
#' reported either way.
#'
#' @param alpha significance level used only to compare conclusions.
#' @inheritParams run_convergence_test
#' @return a list: $markers (the full marker-null result, the one to report),
#'   $sequence (the full uniform-null result), $table (one row per null),
#'   $conclusions_differ, $baseline_ratio (marker-null chance rate divided by
#'   the uniform-null chance rate), $status ("agree", "uniform_only" or
#'   "markers_only") and $message (plain-language summary for display).
compare_null_models = function(cojo_signals, selscan_regions, tolerance,
                                marker_positions, n_perm = 1000, seed = 20260918,
                                alpha = 0.05) {
  markers = run_convergence_test(cojo_signals, selscan_regions, tolerance,
                                 marker_positions, n_perm = n_perm, seed = seed,
                                 null = "markers")
  # the uniform null's warning is suppressed here, since it is run on purpose
  # as a comparison and its result is never the one reported
  sequence = suppressWarnings(
    run_convergence_test(cojo_signals, selscan_regions, tolerance,
                         marker_positions, n_perm = n_perm, seed = seed,
                         null = "sequence"))

  sig_m = markers$empirical_p < alpha
  sig_s = sequence$empirical_p < alpha
  status = if (sig_m == sig_s) "agree" else if (sig_s) "uniform_only" else "markers_only"
  ratio = if (sequence$null_mean_pct > 0) markers$null_mean_pct / sequence$null_mean_pct else NA_real_
  ratio_txt = if (is.na(ratio)) "substantially" else sprintf("%.1f times", ratio)

  fmt_p = function(p) format(signif(p, 3), scientific = FALSE)
  message = switch(status,
    uniform_only = sprintf(paste0(
      "A uniform null would report significant convergence (p = %s) that the ",
      "marker-based null does not support (p = %s). Chance convergence is %s higher once ",
      "permuted signals are restricted to tested markers (%.1f%% against %.1f%%). This ",
      "indicates that markers are concentrated where the selection-scan regions lie. ",
      "The marker-based result is the one reported."),
      fmt_p(sequence$empirical_p), fmt_p(markers$empirical_p), ratio_txt,
      markers$null_mean_pct, sequence$null_mean_pct),
    markers_only = sprintf(paste0(
      "Convergence is significant under the marker-based null (p = %s) but not under a ",
      "uniform null (p = %s). Chance convergence is %s lower once permuted signals ",
      "are restricted to tested markers (%.1f%% against %.1f%%). This happens when the ",
      "selection-scan regions lie where markers are sparse. The uniform null then places ",
      "signals where no marker could have been significant, and so overstates chance. The ",
      "marker-based result is the one reported."),
      fmt_p(markers$empirical_p), fmt_p(sequence$empirical_p),
      if (is.na(ratio) || ratio <= 0) "substantially" else sprintf("%.1f times", 1 / ratio),
      markers$null_mean_pct, sequence$null_mean_pct),
    agree = sprintf(paste0(
      "Both null models lead to the same conclusion (marker-based p = %s, uniform p = %s). ",
      "Chance convergence under the marker-based null is %s the uniform-null rate ",
      "(%.1f%% against %.1f%%)."),
      fmt_p(markers$empirical_p), fmt_p(sequence$empirical_p),
      if (is.na(ratio)) "not comparable to" else sprintf("%.1f times", ratio),
      markers$null_mean_pct, sequence$null_mean_pct))

  table = data.table(
    null_model = c("markers (reported)", "sequence (diagnostic)"),
    observed = c(markers$observed, sequence$observed),
    n_signals = c(markers$n_regions, sequence$n_regions),
    observed_pct = c(markers$observed_pct, sequence$observed_pct),
    null_mean_pct = c(markers$null_mean_pct, sequence$null_mean_pct),
    empirical_p = c(markers$empirical_p, sequence$empirical_p))

  list(markers = markers, sequence = sequence, table = table,
       conclusions_differ = status != "agree", baseline_ratio = ratio,
       status = status, message = message, alpha = alpha)
}

#' Largest marker position on each chromosome, from the GWAS summary
#' statistics alone, used to draw the genome view.
#' @param gwas_std data.table with CHR and POS columns.
derive_chr_lengths = function(gwas_std) {
  gwas_std[, .(max_pos = max(POS, na.rm = TRUE)), by = .(chr = CHR)]
}
