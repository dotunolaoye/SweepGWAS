# Region placement diagnostic: a marker-level logistic model of whether each
# tested marker lies within the distance threshold of a selection-scan region, from
# local marker density, the gap to the nearest marker, allele frequency,
# effect size and association strength, with the local recombination rate
# from a genetic map and gene density from an annotation when the user
# supplies them. It explains a disagreement between the two nulls:
# regions in marker-dense windows make a uniform null understate chance, and
# regions in marker-sparse windows make it overstate chance.
library(data.table)

placement_labels = c(density = "Markers per 100 kb", nn_gap = "Gap to nearest marker (bp)",
                     maf = "Minor allele frequency", effect = "Absolute effect size",
                     logp = "Association strength (-log10 p)", recombination = "Recombination rate (cM/Mb)",
                     genes = "Genes per 100 kb", crossovers = "Crossovers per 100 kb")

# AUC by the Mann-Whitney identity, with counts as doubles so that large
# panels do not overflow R's integer range
auc_mw = function(score, y) {
  r = rank(score); n1 = as.numeric(sum(y)); n0 = as.numeric(sum(!y))
  (sum(r[y]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# Average precision, the area under the precision-recall curve: the mean
# precision at the rank of each marker that is near a region. Unlike the AUC
# it is compared with a baseline, the share of markers near a region, which
# is what a model without information would score; it is more informative
# than the AUC when that share is small.
average_precision = function(score, y) {
  o = order(score, decreasing = TRUE)
  yy = y[o]
  sum((cumsum(yy) / seq_along(yy))[yy]) / sum(yy)
}

# counts of events (crossovers or gene starts) per 100 kb window, joined to
# each marker's window; windows without events get zero
window_counts = function(mk, events, win, name) {
  ev = as.data.table(events)
  setnames(ev, 1:2, c("chr", "pos"))
  ev = ev[, .(n = .N), by = .(chr = as.character(chr), wbin = floor(as.numeric(pos) / win))]
  mk[ev, (name) := i.n, on = .(chr, wbin)]
  mk[is.na(get(name)), (name) := 0L]
  invisible(mk)
}

# Recombination rate of each marker's 100 kb window, as the slope of the
# genetic map in cM per Mb: genetic positions are interpolated linearly at
# the window's ends. Genetic positions must not fall along a chromosome, so
# small reversals left by map estimation are flattened. A window outside the
# span of the map takes its chromosome's average rate, and a chromosome
# missing from the map takes the average rate of the whole map.
window_recombination = function(mk, gmap, win) {
  g = as.data.table(gmap)[, .(chr = as.character(chr), bp = as.numeric(bp), cm = as.numeric(cm))]
  setorder(g, chr, bp)
  g[, cm := cummax(cm), by = chr]
  g = g[, .(cm = mean(cm)), by = .(chr, bp)]
  span = g[, .(avg = if (.N > 1 && max(bp) > min(bp)) (max(cm) - min(cm)) / ((max(bp) - min(bp)) / 1e6) else NA_real_,
               lo = min(bp), hi = max(bp)), by = chr]
  overall = sum(span$avg * (span$hi - span$lo), na.rm = TRUE) / sum((span$hi - span$lo)[!is.na(span$avg)])
  bins = unique(mk[, .(chr, wbin)])
  bins[, rate := {
    gg = g[chr == .BY$chr]
    if (nrow(gg) < 2) rep(NA_real_, .N) else {
      a = pmax(wbin * win, min(gg$bp)); b = pmin((wbin + 1) * win, max(gg$bp))
      ok = b > a
      r = rep(NA_real_, .N)
      r[ok] = (approx(gg$bp, gg$cm, b[ok])$y - approx(gg$bp, gg$cm, a[ok])$y) / ((b[ok] - a[ok]) / 1e6)
      r
    }
  }, by = chr]
  bins = merge(bins, span[, .(chr, avg)], by = "chr", all.x = TRUE)
  bins[is.na(rate), rate := avg]
  bins[is.na(rate), rate := overall]
  mk[bins, recombination := i.rate, on = .(chr, wbin)]
  invisible(mk)
}

#' Fit the region placement diagnostic.
#'
#' @param markers data.table with chr, pos, freq (effect-allele frequency),
#'   beta and p for every tested marker.
#' @param selscan_regions data.table with chr, start and end.
#' @param tolerance the distance threshold in bp.
#' @param genetic_map optional data.frame with chr, bp and cm, from
#'   read_genetic_map(), giving the recombination-rate feature.
#' @param genes optional data.frame whose first two columns are chromosome
#'   and gene start, from read_gene_annotation().
#' @param crossovers optional data.frame of crossover positions (chromosome,
#'   bp), counted per window. The app uses a genetic map instead; this is kept
#'   for the case-study scripts, which used a crossover map.
#' @param win window size in bp for the density features.
#'
#' @return a list: $coefficients (standardised, largest first, with
#'   chromosome-jackknife standard errors and 95% intervals), $feature_means
#'   (mean of each feature near and not near a region, and their ratio),
#'   $auc_by_chr (AUC, average precision and baseline for each withheld
#'   chromosome), $auc, $avg_precision and $ap_baseline (their means),
#'   $n_chr, $pct_near, $n_markers, $features, and $message.
fit_region_placement = function(markers, selscan_regions, tolerance, genetic_map = NULL,
                                genes = NULL, crossovers = NULL, win = 1e5) {
  mk = as.data.table(markers)[, .(chr = as.character(chr), pos = as.numeric(pos),
                                  maf = pmin(as.numeric(freq), 1 - as.numeric(freq)),
                                  effect = abs(as.numeric(beta)),
                                  logp = -log10(pmax(as.numeric(p), 1e-300)))]
  mk = mk[!is.na(pos) & !is.na(maf) & is.finite(logp) & !is.na(effect)]

  sel = as.data.table(selscan_regions)[, .(chr = as.character(chr), start = as.numeric(start),
                                           end = as.numeric(end))]
  setkey(sel, chr, start, end)
  q = mk[, .(chr, lo = pos - tolerance, hi = pos + tolerance)]
  ov = foverlaps(q, sel, by.x = c("chr", "lo", "hi"), by.y = c("chr", "start", "end"),
                 nomatch = NULL, which = TRUE)
  mk[, near := FALSE]
  mk[unique(ov$xid), near := TRUE]

  mk[, wbin := floor(pos / win)]
  mk[, density := .N, by = .(chr, wbin)]
  setorder(mk, chr, pos)
  mk[, gap_prev := pos - shift(pos), by = chr]
  mk[, gap_next := shift(pos, type = "lead") - pos, by = chr]
  mk[, nn_gap := pmin(gap_prev, gap_next, na.rm = TRUE)]
  mk[is.na(nn_gap), nn_gap := median(mk$nn_gap, na.rm = TRUE)]

  feats = c("density", "nn_gap", "maf", "effect", "logp")
  if (!is.null(genetic_map) && nrow(genetic_map) > 1) { window_recombination(mk, genetic_map, win); feats = c(feats, "recombination") }
  if (!is.null(crossovers) && nrow(crossovers) > 0) { window_counts(mk, crossovers, win, "crossovers"); feats = c(feats, "crossovers") }
  if (!is.null(genes) && nrow(genes) > 0) { window_counts(mk, genes, win, "genes"); feats = c(feats, "genes") }

  if (sum(mk$near) == 0 || all(mk$near)) {
    return(list(coefficients = data.table(), feature_means = data.table(), auc_by_chr = data.table(),
                auc = NA_real_, pct_near = 100 * mean(mk$near), n_markers = nrow(mk), features = feats,
                message = "Every marker, or none, lies within the distance threshold of a region. There is nothing for the model to distinguish."))
  }

  fm = rbindlist(lapply(feats, function(f) {
    a = mean(mk[near == TRUE][[f]]); b = mean(mk[near == FALSE][[f]])
    data.table(feature = f, label = placement_labels[[f]], mean_near = a, mean_not_near = b, ratio = a / b)
  }))

  d = mk[, c("chr", "near", feats), with = FALSE]
  for (f in feats) set(d, j = f, value = as.numeric(scale(d[[f]])))
  form = reformulate(feats, "near")
  full = coef(glm(form, data = d, family = binomial()))[feats]

  # The model is refitted with each chromosome withheld in turn. The withheld
  # chromosome measures discrimination, since the window features are shared
  # by neighbouring markers and a random split would overstate it, and the
  # refits give a delete-one-chromosome jackknife for the coefficients. The
  # usual glm standard errors are not used: they treat every marker as
  # independent, which neighbouring markers are not, so they are far too small.
  chrs = sort(unique(d$chr))
  loo = lapply(chrs, function(cc) {
    tr = d$chr != cc
    te = !tr
    m = glm(form, data = d[tr], family = binomial())
    held = if (sum(d$near[te]) == 0 || all(d$near[te])) NULL else {
      sc = predict(m, d[te])
      data.table(chr = cc, n_test = sum(te), test_near = sum(d$near[te]),
                 auc = auc_mw(sc, d$near[te]), avg_precision = average_precision(sc, d$near[te]),
                 baseline = mean(d$near[te]))
    }
    list(coef = coef(m)[feats], held = held)
  })
  auc_by_chr = rbindlist(lapply(loo, `[[`, "held"))
  auc = if (nrow(auc_by_chr)) mean(auc_by_chr$auc) else NA_real_
  ap = if (nrow(auc_by_chr)) mean(auc_by_chr$avg_precision) else NA_real_
  ap_base = if (nrow(auc_by_chr)) mean(auc_by_chr$baseline) else NA_real_

  g = length(chrs)
  jk = do.call(rbind, lapply(loo, `[[`, "coef"))
  jk_se = sqrt((g - 1) / g * colSums(sweep(jk, 2, colMeans(jk))^2))
  tq = qt(0.975, df = g - 1)
  co = data.table(feature = feats, estimate = unname(full), jackknife_se = unname(jk_se[feats]))
  co[, `:=`(ci_low = estimate - tq * jackknife_se, ci_high = estimate + tq * jackknife_se)]
  co[, label := placement_labels[feature]]
  co = co[order(-abs(estimate))]

  list(coefficients = co, feature_means = fm, auc_by_chr = auc_by_chr, auc = auc,
       avg_precision = ap, ap_baseline = ap_base, n_chr = g,
       pct_near = 100 * mean(mk$near), n_markers = nrow(mk), features = feats,
       message = placement_message(co, fm, auc, tolerance))
}

# One or two plain sentences, generated from the result. The wording is
# descriptive: the model shows where regions sit, not why.
placement_message = function(co, fm, auc, tolerance) {
  tol_txt = if (tolerance %% 1000 == 0) paste(tolerance / 1000, "kb") else paste(format(tolerance, big.mark = ","), "bp")
  auc_txt = sprintf("%.3f", auc)
  if (is.na(auc)) {
    return("The model could not be evaluated with chromosomes withheld, because no withheld chromosome had markers both near and away from a region.")
  }
  if (auc < 0.6) {
    return(sprintf(paste0("Whether a marker lies within %s of a selection-scan region is only weakly predictable from ",
                          "these features (AUC %s with each chromosome withheld, where 0.5 is no better than chance). ",
                          "Thus, the regions are not concentrated where markers are unusually dense or sparse. The ",
                          "two null models should give similar chance baselines."), tol_txt, auc_txt))
  }
  top = co$feature[1]
  dens = fm[feature == "density"]
  lead = sprintf(paste0("Whether a marker lies within %s of a selection-scan region is predictable from these features ",
                        "(AUC %s with each chromosome withheld, where 0.5 is no better than chance). The strongest ",
                        "predictor is %s."), tol_txt, auc_txt, tolower(placement_labels[[top]]))
  if (top %in% c("density", "nn_gap") && nrow(dens)) {
    pct = round(100 * abs(dens$ratio - 1))
    if (dens$ratio > 1) {
      return(paste(lead, sprintf(paste0("Markers near a region lie in windows with %d%% more markers than markers ",
                                        "elsewhere. The regions lie where markers are dense. A null that places ",
                                        "signals uniformly along each chromosome therefore tends to understate chance convergence."), pct)))
    } else {
      return(paste(lead, sprintf(paste0("Markers near a region lie in windows with %d%% fewer markers than markers ",
                                        "elsewhere. The regions lie where markers are sparse. A null that places ",
                                        "signals uniformly along each chromosome therefore tends to overstate chance convergence."), pct)))
    }
  }
  lead
}
