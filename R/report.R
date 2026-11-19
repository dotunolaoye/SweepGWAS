# Reports of a finished run, as a PDF and as a single HTML file, so that the
# result can be kept or shared without the app. Both are built in R alone: the
# PDF with the Cairo graphics device and the HTML with the figures embedded as
# images, because the hosted server has neither LaTeX nor a guaranteed pandoc.
# The interactive plots stay in the app.
library(data.table)
library(ggplot2)

# formatted one value at a time, so that a small p-value does not pad the others
fmt_p = function(x) vapply(x, function(v) format(signif(v, 3), scientific = FALSE), "")

# The sentence that opens the Convergence tab, shared with the reports so that
# the app and both reports say exactly the same thing.
convergence_verdict = function(cv, nc = NULL) {
  tol = cv$tolerance_used
  tol_txt = if (tol %% 1000 == 0) paste(tol / 1000, "kb") else paste(format(tol, big.mark = ","), "bp")
  p_low = (sum(cv$null_distribution <= cv$observed) + 1) / (length(cv$null_distribution) + 1)
  lead = sprintf(paste0("%d of %d independent signals (%.1f%%) lie within %s of a selection-scan region, against %.1f%% ",
                        "expected by chance under the marker-based null (p = %s). "),
                 cv$observed, cv$n_regions, cv$observed_pct, tol_txt, cv$null_mean_pct, fmt_p(cv$empirical_p))
  verdict = if (cv$empirical_p < 0.05) {
    "This is more than chance alone produces. It is evidence of convergence at this distance threshold."
  } else if (p_low < 0.05) {
    sprintf(paste0("Fewer signals lie near a region than chance produces (lower-tail p = %s). The result indicates ",
                   "a shortfall rather than convergence."), fmt_p(p_low))
  } else {
    "This is consistent with chance. There is no evidence of convergence beyond what the marker layout alone produces."
  }
  caveat = c(
    if (cv$n_regions < 10) "With fewer than 10 signals, the test has little power. A result consistent with chance is therefore not evidence against convergence.",
    if (!is.null(nc) && isTRUE(nc$conclusions_differ)) "A uniform null would reach a different conclusion. See the null-model check.")
  list(text = paste0(lead, verdict), caveat = paste(caveat, collapse = " "))
}

gg_permutation = function(cv, nc) {
  ink = okabe_ito$text
  h = data.table(count = cv$null_distribution)[, .N, by = count][, null := "Marker-based null (reported)"]
  u = if (!is.null(nc)) data.table(count = nc$sequence$null_distribution)[, .N, by = count] else NULL
  p = ggplot() + geom_col(data = h, aes(count, N), fill = okabe_ito$non_convergent, width = 0.9)
  if (!is.null(u)) {
    full = data.table(count = seq(min(u$count), max(u$count)))
    u = merge(full, u, by = "count", all.x = TRUE)[is.na(N), N := 0]
    p = p + geom_step(data = u[, count := count - 0.5], aes(count, N), colour = okabe_ito$selscan_track, linewidth = 0.6)
  }
  p + geom_vline(xintercept = cv$observed, colour = okabe_ito$convergent, linewidth = 0.9) +
    annotate("label", x = cv$observed + 0.3, y = Inf, vjust = 1.2, hjust = 0, size = 3, colour = ink, fill = "white",
             label = sprintf("Observed: %d of %d\nmarker-based null p = %s%s", cv$observed, cv$n_regions, fmt_p(cv$empirical_p),
                             if (!is.null(nc)) sprintf("\nuniform null p = %s", fmt_p(nc$sequence$empirical_p)) else "")) +
    labs(x = sprintf("Signals within the distance threshold of a region (of %d)", cv$n_regions), y = "Permutations",
         caption = "Bars: marker-based null (reported). Grey outline: uniform null (check). Orange line: observed count.") +
    theme_minimal(base_size = 10) + theme(panel.grid.minor = element_blank(), text = element_text(colour = ink),
                                          plot.caption = element_text(hjust = 0, colour = "#5B6573"))
}

gg_placement = function(co) {
  ink = okabe_ito$text
  co = copy(co)[, label := factor(label, levels = rev(label))]
  co[, dir := ifelse(estimate > 0, "Higher near regions", "Lower near regions")]
  ggplot(co, aes(estimate, label)) +
    geom_vline(xintercept = 0, colour = ink, linewidth = 0.3) +
    geom_col(aes(fill = dir), width = 0.6) +
    geom_errorbar(aes(xmin = ci_low, xmax = ci_high), width = 0.25, orientation = "y", linewidth = 0.35, colour = ink) +
    scale_fill_manual(values = c("Higher near regions" = okabe_ito$non_convergent, "Lower near regions" = okabe_ito$selscan_track)) +
    labs(x = "Standardised coefficient (log-odds per SD), with 95% chromosome-jackknife interval", y = NULL, fill = NULL) +
    theme_minimal(base_size = 10) + theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
                                          legend.position = "top", legend.justification = "left", text = element_text(colour = ink))
}

# Program versions, read from the binaries themselves.
tool_versions = function(plink2_path, gcta_path) {
  first = function(x, pat) { v = grep(pat, x, value = TRUE, ignore.case = TRUE); if (length(v)) trimws(sub("^\\*\\s*", "", v[1])) else "not found" }
  p2 = tryCatch(system2(plink2_path, "--version", stdout = TRUE, stderr = TRUE), error = function(e) "")
  # GCTA prints its version banner and exits with a non-zero status when run without options
  gc = tryCatch(suppressWarnings(system2(gcta_path, stdout = TRUE, stderr = TRUE)), error = function(e) "")
  c(PLINK2 = first(p2, "PLINK"), GCTA = first(gc, "version"))
}

# Everything the reports say, gathered once from the app state.
report_content = function(st, plink2_path, gcta_path) {
  cv = st$convergence; nc = st$null_check; pl = st$placement; cj = st$cojo
  v = convergence_verdict(cv, nc)
  method_txt = if (identical(cj$method, "clumping")) {
    sprintf(paste0("The COJO joint model was judged unstable (median joint-to-marginal effect ratio %.2f, above the tool's stability limit of %s). ",
                   "Independent loci (%d) were therefore taken by LD clumping."), cj$rules$median_ratio, cj$rules$max_median_ratio, cj$n_signals)
  } else {
    sprintf(paste0("COJO selected %d signals, of which %d were significant both jointly and marginally and were retained ",
                   "(median joint-to-marginal effect ratio %.2f, within the tool's stability limit of %s)."), cj$rules$n_cojo, cj$rules$n_retained,
            cj$rules$median_ratio, cj$rules$max_median_ratio)
  }
  placement_txt = if (!is.null(pl) && nrow(pl$coefficients) > 0) paste(pl$message, sprintf(paste0(
    "Average precision %.3f, which is %.1f times its baseline of %.3f. The baseline is the share of markers near a region."),
    pl$avg_precision, pl$avg_precision / pl$ap_baseline, pl$ap_baseline))
    else if (!is.null(pl)) pl$message else ""
  regions_by_method = if ("method" %in% names(st$selscan)) st$selscan[, .N, by = method] else data.table(method = "all", N = nrow(st$selscan))
  sig = cv$regions[, .(Chromosome = chr, Signal = lead_snp, Position = format(lead_pos, big.mark = ",", trim = TRUE),
                       `Within distance threshold` = ifelse(convergent, "yes", "no"),
                       `Distance to nearest region (bp)` = ifelse(is.na(nearest_selscan_dist), "none on chromosome",
                                                                  format(round(nearest_selscan_dist), big.mark = ",", trim = TRUE)))]
  pkgs = c("shiny", "bslib", "shinyjs", "plotly", "data.table", "DT", "ggplot2", "patchwork")
  versions = c(R = paste(R.version$major, R.version$minor, sep = "."),
               vapply(pkgs, function(p) as.character(packageVersion(p)), ""), tool_versions(plink2_path, gcta_path))
  list(
    title = "SweepGWAS report",
    date = format(Sys.time(), "%Y-%m-%d %H:%M %Z"),
    verdict = v$text, caveat = v$caveat, method = method_txt,
    null_check = if (!is.null(nc)) nc$message else "",
    placement = placement_txt,
    settings = data.table(Setting = c("Distance threshold", "Permutations", "COJO significance threshold", "Markers supplied",
                                      "Markers after the allele check", "Independent signals tested", "Selection-scan regions"),
                          Value = c(paste(format(cv$tolerance_used, big.mark = ","), "bp"), format(length(cv$null_distribution), big.mark = ","),
                                    format(cj$threshold_used, digits = 3), format(nrow(st$markers), big.mark = ","),
                                    format(nrow(cj$allele_check$matched), big.mark = ","), cv$n_regions,
                                    paste0(format(nrow(st$selscan), big.mark = ","), " (",
                                           paste(sprintf("%s %s", regions_by_method$method, format(regions_by_method$N, big.mark = ",", trim = TRUE)), collapse = "; "), ")"))),
    null_table = if (!is.null(nc)) data.table(`Null model` = c("Drawn from tested markers (reported)", "Drawn uniformly along each chromosome (check)"),
                                              `Expected by chance` = sprintf("%.1f%%", nc$table$null_mean_pct),
                                              `p-value` = fmt_p(nc$table$empirical_p)) else NULL,
    signals = sig,
    stages = st$stages %||% "",
    versions = data.table(Software = names(versions), Version = unname(versions)))
}

report_figures = function(st) {
  chrs = chr_order(st$markers$chr)
  figs = list()
  g = gg_convergence_manhattan(st$markers, st$convergence$regions, st$selscan, st$cojo$threshold_used)
  figs[["Manhattan plot, all chromosomes"]] = g$plot
  thinned = g$thinned
  for (ch in chrs) {
    gc = gg_convergence_manhattan(st$markers, st$convergence$regions, st$selscan, st$cojo$threshold_used, chr = ch)
    figs[[paste("Chromosome", ch)]] = gc$plot
    thinned = thinned || gc$thinned
  }
  figs[["Permutation distributions under both null models"]] = gg_permutation(st$convergence, st$null_check)
  if (!is.null(st$placement) && nrow(st$placement$coefficients) > 0) figs[["Region placement diagnostic"]] = gg_placement(st$placement$coefficients)
  list(figs = figs, thinned = thinned)
}

thinning_note = "For speed, markers with p above 0.001 were randomly thinned in some panels. Every marker below 0.001 is shown."

write_report_pdf = function(file, st, plink2_path, gcta_path) {
  rc = report_content(st, plink2_path, gcta_path)
  rf = report_figures(st)
  grDevices::cairo_pdf(file, width = 11, height = 8.5, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)

  # text pages: each paragraph is wrapped and placed down the page, and a new
  # page starts when the space runs out
  text_pages = function(blocks) {
    y = 0.95; grid::grid.newpage()
    for (bl in blocks) {
      lines = if (isTRUE(bl$mono)) bl$text else unlist(lapply(bl$text, strwrap, width = if (isTRUE(bl$head)) 90 else 125))
      step = if (isTRUE(bl$head)) 0.045 else 0.028
      for (ln in lines) {
        if (y < 0.06) { grid::grid.newpage(); y = 0.95 }
        grid::grid.text(ln, x = 0.05, y = y, just = c("left", "top"),
                        gp = grid::gpar(fontsize = if (isTRUE(bl$head)) 15 else 10.5, fontface = if (isTRUE(bl$head)) "bold" else "plain",
                                        fontfamily = if (isTRUE(bl$mono)) "mono" else "sans", col = okabe_ito$text))
        y = y - step
      }
      y = y - 0.015
    }
  }
  tab_lines = function(d) {
    d = as.data.frame(lapply(d, as.character), check.names = FALSE)
    w = pmax(nchar(names(d)), vapply(d, function(x) max(nchar(x), 0L), 0L))
    row = function(v) paste(mapply(function(x, ww) formatC(x, width = -ww), v, w), collapse = "  ")
    c(row(names(d)), row(strrep("-", w)), apply(d, 1, row))
  }
  text_pages(list(
    list(text = rc$title, head = TRUE), list(text = rc$date),
    list(text = "Result", head = TRUE), list(text = c(rc$verdict, rc$caveat)),
    list(text = "Independent signals", head = TRUE), list(text = rc$method),
    list(text = "Null-model check", head = TRUE), list(text = rc$null_check),
    if (!is.null(rc$null_table)) list(text = tab_lines(rc$null_table), mono = TRUE),
    list(text = "Region placement diagnostic", head = TRUE), list(text = rc$placement),
    list(text = "Settings", head = TRUE), list(text = tab_lines(rc$settings), mono = TRUE)))
  for (nm in names(rf$figs)) print(rf$figs[[nm]] + patchwork::plot_annotation(title = nm,
    caption = if (rf$thinned && grepl("Manhattan|Chromosome", nm)) thinning_note else NULL))
  text_pages(list(list(text = "Independent signals and their nearest region", head = TRUE),
                  list(text = tab_lines(rc$signals), mono = TRUE),
                  list(text = "Run time and software", head = TRUE), list(text = rc$stages),
                  list(text = tab_lines(rc$versions), mono = TRUE)))
  invisible(file)
}

write_report_html = function(file, st, plink2_path, gcta_path) {
  rc = report_content(st, plink2_path, gcta_path)
  rf = report_figures(st)
  esc = function(x) htmltools::htmlEscape(x)
  img = function(p, w = 11, h = 6.5) {
    f = tempfile(fileext = ".png"); on.exit(unlink(f))
    ggsave(f, p, width = w, height = h, dpi = 110, bg = "white")
    sprintf('<img alt="" src="data:image/png;base64,%s">', jsonlite::base64_enc(readBin(f, "raw", file.info(f)$size)))
  }
  tbl = function(d) paste0("<table><thead><tr>", paste0("<th>", esc(names(d)), "</th>", collapse = ""), "</tr></thead><tbody>",
    paste0(apply(as.data.frame(lapply(d, as.character), check.names = FALSE), 1, function(r) paste0("<tr>", paste0("<td>", esc(r), "</td>", collapse = ""), "</tr>")), collapse = ""),
    "</tbody></table>")
  para = function(x) paste0("<p>", esc(x), "</p>", collapse = "")
  figs = paste0(vapply(names(rf$figs), function(nm) {
    h = if (grepl("placement", nm)) 4.5 else if (grepl("Permutation", nm)) 4.5 else 6.5
    paste0("<h3>", esc(nm), "</h3>", img(rf$figs[[nm]], h = h))
  }, ""), collapse = "\n")
  html = paste0('<!DOCTYPE html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">',
    "<title>", esc(rc$title), "</title><style>",
    "body{font-family:-apple-system,Segoe UI,Helvetica,Arial,sans-serif;color:#1B2430;background:#fff;max-width:1100px;margin:0 auto;padding:24px 16px;line-height:1.5}",
    "h1{font-size:1.6rem;margin:0 0 .2rem}h2{font-size:1.2rem;margin-top:2rem;border-bottom:1px solid #E4E1DB;padding-bottom:.2rem}h3{font-size:1rem;margin-top:1.4rem}",
    ".date{color:#5B6573;margin:0 0 1rem}.verdict{background:#F5F4F0;border:1px solid #E4E1DB;border-radius:6px;padding:.8rem 1rem}",
    "img{max-width:100%;height:auto}table{border-collapse:collapse;font-size:.9rem;margin:.5rem 0}th,td{text-align:left;padding:.25rem .8rem .25rem 0;border-bottom:1px solid #E8EAED}",
    ".note{color:#5B6573;font-size:.9rem}</style></head><body>",
    "<h1>", esc(rc$title), '</h1><p class="date">', esc(rc$date), "</p>",
    '<h2>Result</h2><div class="verdict">', para(rc$verdict), if (nzchar(rc$caveat)) para(rc$caveat) else "", "</div>",
    "<h2>Independent signals</h2>", para(rc$method),
    "<h2>Null-model check</h2>", para(rc$null_check), if (!is.null(rc$null_table)) tbl(rc$null_table) else "",
    "<h2>Region placement diagnostic</h2>", para(rc$placement),
    "<h2>Settings</h2>", tbl(rc$settings),
    "<h2>Figures</h2>", if (rf$thinned) paste0('<p class="note">', esc(thinning_note), "</p>") else "", figs,
    "<h2>Independent signals and their nearest region</h2>", tbl(rc$signals),
    "<h2>Run time and software</h2>", para(rc$stages), tbl(rc$versions),
    "</body></html>")
  writeLines(html, file, useBytes = TRUE)
  invisible(file)
}
