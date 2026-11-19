# Manhattan view of the GWAS with the selection-scan regions beneath it, one
# track per method, as a genome browser shows annotation tracks. The
# independent signals are marked on the Manhattan plot, orange if they lie
# within the distance threshold of a region and blue otherwise, with a faint line from
# each signal down through the tracks. plot_convergence_manhattan() draws the
# interactive version for the app and gg_convergence_manhattan() the static
# version for the reports; both take the same inputs.
library(data.table)
library(plotly)
library(ggplot2)
library(patchwork)

# Chromosome names sorted as a reader expects: 1, 2, ..., 10, then X and Y,
# with a "chr" prefix ignored. Sorted as plain text, chr10 would come before chr2.
chr_order = function(x) {
  x = unique(as.character(x))
  core = sub("^chr", "", x, ignore.case = TRUE)
  num = suppressWarnings(as.numeric(core))
  x[order(is.na(num), num, core)]
}

# Keeps every marker with p below keep_p and a random sample of the rest, so a
# panel of millions of markers still draws quickly. The associated markers,
# which are the ones a reader looks at, are never dropped.
thin_markers = function(m, max_points = 250000, keep_p = 1e-3, seed = 1) {
  if (nrow(m) <= max_points) return(list(markers = m, thinned = FALSE))
  strong = m[p < keep_p]
  rest = m[!(p < keep_p)]
  set.seed(seed)
  n_rest = max(0, max_points - nrow(strong))
  if (n_rest < nrow(rest)) rest = rest[sample.int(nrow(rest), n_rest)]
  list(markers = rbind(strong, rest), thinned = TRUE)
}

region_hover = function(chr, start, end, method) {
  span = ifelse(start == end, sprintf("%s bp", format(start, big.mark = ",", trim = TRUE)),
                sprintf("%s to %s bp", format(start, big.mark = ",", trim = TRUE), format(end, big.mark = ",", trim = TRUE)))
  sprintf("<b>Selection-scan region</b><br>Chr %s : %s%s", chr, span,
          ifelse(is.na(method), "", paste0("<br>Method: ", method)))
}

signal_hover = function(d) {
  sprintf("<b>%s</b><br>Chr %s : %s bp<br>p = %s<br>%s", d$lead_snp, d$chr, format(d$lead_pos, big.mark = ","),
          format(signif(d$p, 3)),
          ifelse(d$convergent, "Within the distance threshold of a selection-scan region",
                 paste0("Nearest selection-scan region: ",
                        ifelse(is.na(d$nearest_selscan_dist), "none on this chromosome",
                               paste(format(round(d$nearest_selscan_dist), big.mark = ","), "bp")))))
}

# Shared by both versions: positions on one x axis, in Mb for a single
# chromosome or laid end to end for the whole genome, and the tracks.
manhattan_data = function(markers, convergence_regions, selscan_regions, chr = NULL, max_points = 250000) {
  m = as.data.table(markers)[, .(chr = as.character(chr), pos = as.numeric(pos), p = as.numeric(p))]
  m = m[!is.na(pos) & !is.na(p) & p > 0]
  sig = as.data.table(convergence_regions)[, chr := as.character(chr)]
  # one p-value per position, in case two markers share a position
  sig = merge(sig, m[, .(p = min(p)), by = .(chr, lead_pos = pos)], by = c("chr", "lead_pos"), all.x = TRUE)
  sel = as.data.table(selscan_regions)[, chr := as.character(chr)]
  if (!"method" %in% names(sel)) sel[, method := "Selection-scan regions"]
  sel[, method := as.character(method)]
  methods = sort(unique(sel$method))
  sel[, track := length(methods) + 1 - match(method, methods)]

  chrs = chr_order(m$chr)
  if (is.null(chr)) {
    len = m[, .(len = max(pos)), by = chr][match(chrs, chr)]
    gap = max(len$len) * 0.02
    len[, offset := c(0, cumsum(len + gap))[seq_len(.N)]]
    off = setNames(len$offset, len$chr)
    m[, x := (pos + off[chr]) / 1e6]
    sig = sig[chr %in% chrs][, x := (lead_pos + off[chr]) / 1e6]
    sel = sel[chr %in% chrs][, x := ((start + end) / 2 + off[chr]) / 1e6]
    m[, band := match(chr, chrs) %% 2 == 0]
    axis = list(ticks = (len$offset + len$len / 2) / 1e6, labels = len$chr, title = "Chromosome",
                range = c(0, (max(len$offset + len$len)) / 1e6))
  } else {
    this = as.character(chr)
    m = m[chr == this][, x := pos / 1e6][, band := FALSE]
    sig = sig[chr == this][, x := lead_pos / 1e6]
    sel = sel[chr == this][, x := (start + end) / 2e6]
    axis = list(ticks = NULL, labels = NULL, title = sprintf("Position on chromosome %s (Mb)", this),
                range = c(0, max(m$pos) / 1e6))
  }
  th = thin_markers(m, max_points)
  sel[, hover := region_hover(chr, start, end, method)]
  list(markers = th$markers, thinned = th$thinned, n_markers = nrow(m), signals = sig, regions = sel,
       methods = methods, axis = axis)
}

plot_convergence_manhattan = function(markers, convergence_regions, selscan_regions, threshold, chr = NULL) {
  d = manhattan_data(markers, convergence_regions, selscan_regions, chr)
  m = d$markers; sig = d$signals; sel = d$regions
  nc = sig[convergent == FALSE]; cv = sig[convergent == TRUE]

  top = plot_ly()
  for (bd in c(FALSE, TRUE)) {
    mm = m[band == bd]
    if (nrow(mm)) top = top |> add_trace(type = "scattergl", mode = "markers", x = mm$x, y = -log10(mm$p),
      marker = list(color = if (bd) "#A7AFB9" else "#C9CED6", size = 3), hoverinfo = "skip", showlegend = FALSE)
  }
  if (nrow(nc)) top = top |> add_trace(type = "scatter", mode = "markers", x = nc$x, y = -log10(nc$p),
    marker = list(color = okabe_ito$non_convergent, size = 10, symbol = "circle", line = list(color = "white", width = 1.2)),
    name = b("Independent signal, no region within the distance threshold"), text = signal_hover(nc), hoverinfo = "text")
  if (nrow(cv)) top = top |> add_trace(type = "scatter", mode = "markers", x = cv$x, y = -log10(cv$p),
    marker = list(color = okabe_ito$convergent, size = 13, symbol = "diamond", line = list(color = "white", width = 1.2)),
    name = b("Independent signal within the distance threshold of a region"), text = signal_hover(cv), hoverinfo = "text")

  bottom = plot_ly() |> add_trace(type = "scatter", mode = "markers", x = sel$x, y = sel$track,
    marker = list(color = okabe_ito$selscan_track, size = 14, opacity = 0.45, symbol = "line-ns-open",
                  line = list(color = okabe_ito$selscan_track, width = 1.1)),
    name = b("Selection-scan region"), text = sel$hover, hoverinfo = "text")

  guides = lapply(seq_len(nrow(sig)), function(i) list(type = "line", xref = "x", yref = "paper", layer = "below",
    x0 = sig$x[i], x1 = sig$x[i], y0 = 0, y1 = 1, opacity = 0.35,
    line = list(color = if (sig$convergent[i]) okabe_ito$convergent else okabe_ito$non_convergent, width = 1, dash = "dot")))
  thr_line = list(list(type = "line", xref = "paper", yref = "y", x0 = 0, x1 = 1,
                       y0 = -log10(threshold), y1 = -log10(threshold), line = list(color = okabe_ito$text, width = 1, dash = "dash")))

  xaxis = list(title = list(text = b(d$axis$title), font = axis_title_font), range = d$axis$range,
               showgrid = FALSE, zeroline = FALSE, tickfont = tick_font)
  if (!is.null(d$axis$ticks)) xaxis = c(xaxis, list(tickvals = d$axis$ticks, ticktext = b(d$axis$labels), tickangle = 0))
  else xaxis = c(xaxis, list(tickformat = ",~f"))

  p = subplot(top, bottom, nrows = 2, heights = c(0.7, 0.3), shareX = TRUE, margin = 0.03) |> layout(
    shapes = c(guides, thr_line),
    xaxis = xaxis,
    yaxis = list(title = list(text = b("-log<sub>10</sub>(p)"), font = axis_title_font), zeroline = FALSE, tickfont = tick_font),
    yaxis2 = list(title = "", tickvals = seq_along(d$methods), ticktext = b(rev(d$methods)), range = c(0.4, length(d$methods) + 0.6),
                  fixedrange = TRUE, showgrid = FALSE, zeroline = FALSE, automargin = TRUE, tickfont = tick_font),
    legend = list(orientation = "h", x = 0, y = 1.1, xanchor = "left"),
    margin = list(l = 10, r = 20, t = 50, b = 55), hovermode = "closest")
  tidy_config(base_plot_theme(p), if (is.null(chr)) "manhattan_all" else sprintf("manhattan_chr%s", chr))
}

gg_convergence_manhattan = function(markers, convergence_regions, selscan_regions, threshold, chr = NULL, max_points = 150000) {
  d = manhattan_data(markers, convergence_regions, selscan_regions, chr, max_points)
  m = d$markers; sig = d$signals; sel = d$regions
  ink = okabe_ito$text
  sig[, status := factor(ifelse(convergent, "Signal within the distance threshold of a region", "Signal, no region within the distance threshold"),
                         levels = c("Signal within the distance threshold of a region", "Signal, no region within the distance threshold"))]
  xs = if (!is.null(d$axis$ticks)) scale_x_continuous(breaks = d$axis$ticks, labels = d$axis$labels, limits = d$axis$range, expand = c(0.005, 0))
       else scale_x_continuous(limits = d$axis$range, expand = c(0.005, 0))
  base = theme_minimal(base_size = 10) + theme(panel.grid.minor = element_blank(), text = element_text(colour = ink),
    axis.text = element_text(colour = ink), legend.position = "top", legend.justification = "left", legend.title = element_blank())
  sig_col = ifelse(sig$convergent, okabe_ito$convergent, okabe_ito$non_convergent)
  top = ggplot(m, aes(x, -log10(p))) +
    geom_point(aes(colour = band), size = 0.3, show.legend = FALSE) +
    scale_colour_manual(values = c(`FALSE` = "#C9CED6", `TRUE` = "#A7AFB9")) +
    geom_hline(yintercept = -log10(threshold), colour = ink, linewidth = 0.3, linetype = "dashed") +
    geom_vline(xintercept = sig$x, colour = sig_col, linewidth = 0.3, alpha = 0.5, linetype = "dotted") +
    geom_point(data = sig, aes(x, -log10(p), fill = status, shape = status), size = 2.6, stroke = 0.5, colour = "white") +
    scale_fill_manual(values = c(okabe_ito$convergent, okabe_ito$non_convergent), drop = FALSE) +
    scale_shape_manual(values = c(23, 21), drop = FALSE) +
    xs + labs(x = NULL, y = expression(-log[10](italic(p)))) + base +
    theme(axis.text.x = element_blank(), panel.grid.major.x = element_blank())
  bottom = ggplot(sel) +
    geom_linerange(aes(x = x, ymin = track - 0.38, ymax = track + 0.38), colour = okabe_ito$selscan_track, linewidth = 0.25) +
    geom_vline(xintercept = sig$x, colour = sig_col, linewidth = 0.35, alpha = 0.7) +
    scale_y_continuous(breaks = seq_along(d$methods), labels = rev(d$methods), limits = c(0.4, length(d$methods) + 0.6)) +
    xs + labs(x = d$axis$title, y = NULL) + base + theme(panel.grid.major = element_blank())
  list(plot = top / bottom + plot_layout(heights = c(2.4, 1)), thinned = d$thinned, n_markers = d$n_markers)
}
