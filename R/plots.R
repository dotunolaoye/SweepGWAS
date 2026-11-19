# The app's interactive plots, besides the Manhattan view in manhattan.R:
# 1. plot_cojo_scatter(): marginal against joint significance of the COJO
#    signals, which shows what conditioning did and which signals were
#    removed.
# 2. plot_permutation_null(): the permutation distribution with the
#    observed value marked, so the p-value can be read as a distribution.
#
# Colours are fixed to colour-blind-safe Okabe-Ito values rather than
# plotly's defaults: orange #E69F00 for convergent signals, blue #0072B2 for
# the others, and grey #7A828C for selection-scan regions. Other distinctions
# use marker shape, so the legend needs only two colours.
library(data.table)
library(plotly)

okabe_ito = list(
  convergent = "#E69F00",
  non_convergent = "#0072B2",
  selscan_track = "#7A828C",
  band = "rgba(27,36,48,0.04)",
  grid = "#E8EAED",
  text = "#1B2430",
  muted = "#6B7280"
)

plot_font = list(family = "Inter, -apple-system, Segoe UI, Helvetica, Arial, sans-serif",
                 color = okabe_ito$text, size = 13)
# Axis titles, tick labels and legend text are bold, so the figures stay
# legible when reduced for a manuscript or slide.
#
# Bold is applied with <b> tags because the plotly.js bundled with the R
# plotly package ignores the font weight setting. Numeric axes therefore get
# explicit tick values and labels, since automatic tick labels cannot carry
# markup.
axis_title_font = list(size = 13, color = okabe_ito$text)
tick_font       = list(size = 12, color = okabe_ito$text)
legend_font     = list(size = 12, color = okabe_ito$text)

b = function(x) paste0("<b>", x, "</b>")

#' Explicit, bold tick labels for a numeric axis over `rng`.
#' @param fmt formatter applied to each tick value before bolding.
bold_numeric_ticks = function(rng, n = 6, fmt = function(v) format(v, big.mark = ",", trim = TRUE)) {
  v = pretty(rng, n = n)
  v = v[v >= min(rng) & v <= max(rng)]
  list(tickmode = "array", tickvals = v, ticktext = b(fmt(v)))
}

# Keep only the plot toolbar buttons that are useful for these figures.
tidy_config = function(p, filename) {
  config(p, displaylogo = FALSE, responsive = TRUE, toImageButtonOptions = list(format = "png", scale = 3, filename = filename),
         modeBarButtonsToRemove = list("select2d", "lasso2d", "autoScale2d",
                                       "hoverClosestCartesian", "hoverCompareCartesian",
                                       "toggleSpikelines"))
}

base_plot_theme = function(p) {
  p |> layout(
    font = plot_font,
    paper_bgcolor = "white",
    plot_bgcolor = "white",
    xaxis = list(gridcolor = okabe_ito$grid, zerolinecolor = okabe_ito$grid,
                 linecolor = okabe_ito$grid, ticks = "outside", tickcolor = okabe_ito$grid,
                 tickfont = tick_font),
    yaxis = list(gridcolor = okabe_ito$grid, zerolinecolor = okabe_ito$grid,
                 linecolor = okabe_ito$grid, ticks = "outside", tickcolor = okabe_ito$grid,
                 tickfont = tick_font),
    legend = list(bgcolor = "rgba(255,255,255,0.92)", bordercolor = okabe_ito$grid,
                  borderwidth = 1, font = legend_font),
    hoverlabel = list(font = list(family = plot_font$family, size = 12),
                      bordercolor = okabe_ito$grid)
  )
}

#' Marginal (single-marker GWAS) against joint (COJO) -log10(p) for each COJO
#' signal. Points near the diagonal kept their significance after
#' conditioning. Points far above it gained significance only in the joint
#' model, which is usually a sign of correlated markers given large
#' offsetting effects rather than of a masked real effect, and such signals
#' are removed by the marginal-and-joint rule.
#'
#' @param cojo_signals GCTA's COJO output with p and pJ columns, normally
#'   run_cojo()$cojo_output.
#' @param convergence_regions optional $regions from run_convergence_test();
#'   when supplied, tested signals are coloured by whether they converge.
#' @param method "cojo" when the retained COJO signals were tested, or
#'   "clumping" when LD-clumped loci replaced them; with a `retained` column
#'   present, signals removed by the marginal-and-joint rule are drawn apart
#'   so they are not mistaken for tested, non-convergent signals.
plot_cojo_scatter = function(cojo_signals, convergence_regions = NULL, method = "cojo") {
  d = copy(cojo_signals)
  if (!"retained" %in% names(d)) d[, retained := TRUE]
  d[, marginal_log10p := -log10(p)]
  d[, joint_log10p := -log10(pJ)]

  if (!is.null(convergence_regions)) {
    cr = unique(as.data.table(convergence_regions)[, .(SNP = lead_snp, convergent)])
    d = merge(d, cr, by = "SNP", all.x = TRUE)
    d[is.na(convergent), convergent := FALSE]
  } else {
    d[, convergent := NA]
  }

  lim = max(c(d$marginal_log10p, d$joint_log10p), na.rm = TRUE) * 1.05
  hover = function(x) sprintf("<b>%s</b><br>Chr %s : %s bp<br>Marginal p = %.2e<br>Joint p = %.2e",
                               x$SNP, x$Chr, format(x$bp, big.mark = ","), x$p, x$pJ)

  p = plot_ly() |>
    add_segments(x = 0, xend = lim, y = 0, yend = lim, inherit = FALSE,
                 line = list(color = okabe_ito$grid, width = 1.5, dash = "dot"),
                 name = b("Marginal = joint"), hoverinfo = "skip")

  # Signals the test did not use: joint-only signals removed by the rule, or
  # every COJO signal when LD-clumped loci were tested instead
  unused = if (method == "clumping") d else d[retained == FALSE]
  if (nrow(unused) > 0) {
    p = p |> add_markers(data = unused, x = ~marginal_log10p, y = ~joint_log10p,
      # open symbols take their outline colour from marker$color
      marker = list(color = okabe_ito$selscan_track, size = 9, symbol = "circle-open",
                    line = list(width = 1.6)),
      text = hover(unused), hoverinfo = "text",
      name = if (method == "clumping") b("COJO signal (not used; LD-clumped loci tested)") else b("Removed (significant only jointly)"))
  }
  d = if (method == "clumping") d[0] else d[retained == TRUE]

  if (nrow(d) == 0) {
    # nothing further to draw
  } else if (all(is.na(d$convergent))) {
    p = p |> add_markers(data = d, x = ~marginal_log10p, y = ~joint_log10p,
      marker = list(color = okabe_ito$non_convergent, size = 10, opacity = 0.85,
                    line = list(color = "white", width = 1.2)),
      text = hover(d), hoverinfo = "text", name = b("COJO signal"))
  } else {
    for (grp in c(FALSE, TRUE)) {
      sub = d[convergent == grp]
      if (nrow(sub) == 0) next
      p = p |> add_markers(data = sub, x = ~marginal_log10p, y = ~joint_log10p,
        marker = list(color = if (grp) okabe_ito$convergent else okabe_ito$non_convergent,
                      size = if (grp) 12 else 10, opacity = 0.9,
                      symbol = if (grp) "diamond" else "circle",
                      line = list(color = "white", width = 1.2)),
        text = hover(sub), hoverinfo = "text",
        name = if (grp) b("Convergent") else b("Non-convergent"))
    }
  }

  tk = bold_numeric_ticks(c(0, lim), n = 6, fmt = function(v) format(v, trim = TRUE))
  p = p |> layout(
    xaxis = c(list(title = list(text = b("-log<sub>10</sub>(p), marginal (single-SNP GWAS)"),
                                font = axis_title_font),
                   range = c(0, lim), zeroline = FALSE, tickfont = tick_font), tk),
    yaxis = c(list(title = list(text = b("-log<sub>10</sub>(p), joint (COJO-conditional)"),
                                font = axis_title_font),
                   range = c(0, lim), zeroline = FALSE, tickfont = tick_font), tk),
    legend = list(orientation = "h", x = 0, y = 1.1, xanchor = "left"),
    margin = list(l = 70, r = 20, t = 50, b = 60),
    hovermode = "closest"
  )
  tidy_config(base_plot_theme(p), "cojo_marginal_vs_joint")
}

#' Histogram of the permutation distribution with the observed value marked,
#' for the output of either run_convergence_test() or
#' run_phenotype_permutation(), whose fields are passed explicitly.
#'
#' @param null_distribution numeric vector, the permutation replicates.
#' @param observed the real (unpermuted) value.
#' @param empirical_p the already-computed empirical p-value, shown as an
#'   annotation.
#' @param x_label axis label, since the two callers count different things
#'   (convergent signals or significant markers).
plot_permutation_null = function(null_distribution, observed, empirical_p,
                                   x_label = "Value under permutation",
                                   uniform_null = NULL, uniform_p = NA_real_) {
  d = data.table(x = null_distribution)
  tail_n = sum(null_distribution >= observed)
  span = range(c(null_distribution, uniform_null, observed))
  # tallest bar, used to place explicit (bold) y ticks
  ymax = max(hist(null_distribution, plot = FALSE,
                  breaks = seq(min(null_distribution) - 0.5,
                               max(null_distribution) + 0.5, by = 1))$counts)
  if (!is.null(uniform_null)) {
    # the uniform null as a grey outline of counts per value, drawn over the
    # reported bars so the two distributions can be compared on one axis
    u = data.table(x = uniform_null)[, .(n = .N), by = x][order(x)]
    u = merge(data.table(x = seq(min(u$x) - 1, max(u$x) + 1)), u, by = "x", all.x = TRUE)
    u[is.na(n), n := 0L]
    ymax = max(ymax, u$n)
  }

  p = plot_ly(d, x = ~x, type = "histogram",
              marker = list(color = okabe_ito$non_convergent, opacity = 0.8,
                            line = list(color = "white", width = 0.6)),
              name = b(if (is.null(uniform_null))
                sprintf("Permutation null (n=%s)", format(length(null_distribution), big.mark = ","))
                else sprintf("Drawn from tested markers (reported), p = %s", format(signif(empirical_p, 3), scientific = FALSE))),
              hovertemplate = "%{y} replicates at %{x}<extra></extra>")
  if (!is.null(uniform_null)) {
    p = p |> add_trace(data = u, x = ~x, y = ~n, type = "scatter", mode = "lines", inherit = FALSE,
                       line = list(shape = "hvh", color = okabe_ito$selscan_track, width = 2),
                       name = b(sprintf("Drawn uniformly along chromosomes (check), p = %s",
                                        format(signif(uniform_p, 3), scientific = FALSE))),
                       hovertemplate = "Uniform null: %{y} replicates at %{x}<extra></extra>")
  }
  p = p |>
    layout(
      xaxis = c(list(title = list(text = b(x_label), font = axis_title_font),
                     zeroline = FALSE, tickfont = tick_font),
                bold_numeric_ticks(span, n = 8, fmt = function(v) format(v, trim = TRUE))),
      # explicit y ticks need the bar heights, so the binning is computed here
      yaxis = c(list(title = list(text = b("Permutation replicates"),
                                  font = axis_title_font),
                     zeroline = FALSE, tickfont = tick_font,
                     range = c(0, ymax * 1.08)),
                bold_numeric_ticks(c(0, ymax), n = 5,
                                   fmt = function(v) format(round(v), trim = TRUE))),
      shapes = list(
        # shade the values at or above the observed one, which the p-value counts
        list(type = "rect", xref = "x", yref = "paper", layer = "below",
             x0 = observed, x1 = span[2] + diff(span) * 0.04, y0 = 0, y1 = 1,
             fillcolor = "rgba(230,159,0,0.10)", line = list(width = 0)),
        list(type = "line", x0 = observed, x1 = observed, y0 = 0, y1 = 1,
             xref = "x", yref = "paper",
             line = list(color = okabe_ito$convergent, width = 2.5))
      ),
      # placed just inside the top of the plot, since above it the label is clipped
      annotations = list(list(
        x = observed, y = 0.97, yref = "paper", yanchor = "top",
        xanchor = if (observed > mean(span)) "right" else "left",
        text = sprintf("<b>Observed = %.3g</b><br>%s of %s replicates at or above<br>p = %.3g",
                        observed, format(tail_n, big.mark = ","),
                        format(length(null_distribution), big.mark = ","), empirical_p),
        showarrow = FALSE, align = if (observed > mean(span)) "right" else "left",
        font = list(color = okabe_ito$text, size = 11.5),
        bgcolor = "rgba(255,255,255,0.88)", borderpad = 4
      )),
      legend = list(orientation = "h", x = 0, y = 1.02, xanchor = "left", yanchor = "bottom"),
      margin = list(l = 70, r = 25, t = if (is.null(uniform_null)) 55 else 75, b = 60),
      bargap = 0.06
    )
  tidy_config(base_plot_theme(p), "permutation_null")
}

#' Standardised coefficients of the region placement diagnostic, largest
#' first, as horizontal bars with 95% chromosome-jackknife intervals. Positive
#' values mean the feature is higher for markers near a selection-scan region.
#'
#' @param coefficients the $coefficients data.table from fit_region_placement().
plot_placement_coefficients = function(coefficients) {
  d = copy(as.data.table(coefficients))[order(abs(estimate))]
  d[, label := factor(label, levels = label)]
  d[, colour := ifelse(estimate >= 0, okabe_ito$non_convergent, okabe_ito$selscan_track)]
  has_ci = "ci_low" %in% names(d)
  hover = if (has_ci) sprintf("<b>%s</b><br>Standardised coefficient = %.3f<br>95%% interval %.3f to %.3f",
                              d$label, d$estimate, d$ci_low, d$ci_high)
          else sprintf("<b>%s</b><br>Standardised coefficient = %.3f", d$label, d$estimate)
  lim = max(abs(c(d$estimate, if (has_ci) c(d$ci_low, d$ci_high)))) * 1.1
  tk = bold_numeric_ticks(c(-lim, lim), n = 6, fmt = function(v) format(round(v, 2), trim = TRUE))
  p = plot_ly(d, x = ~estimate, y = ~label, type = "bar", orientation = "h",
              marker = list(color = ~colour), text = hover, hoverinfo = "text", textposition = "none",
              error_x = if (has_ci) list(type = "data", symmetric = FALSE, array = d$ci_high - d$estimate,
                                         arrayminus = d$estimate - d$ci_low, color = okabe_ito$text,
                                         thickness = 1.5, width = 5)) |>
    layout(
      xaxis = c(list(title = list(text = b("Standardised coefficient (log-odds per SD)"), font = axis_title_font),
                     range = c(-lim, lim), zeroline = TRUE, zerolinecolor = okabe_ito$text,
                     zerolinewidth = 1.2, tickfont = tick_font), tk),
      yaxis = list(title = "", tickfont = tick_font, ticktext = b(levels(d$label)),
                   tickvals = levels(d$label), automargin = TRUE),
      showlegend = FALSE, bargap = 0.35,
      margin = list(l = 10, r = 20, t = 10, b = 60))
  tidy_config(base_plot_theme(p), "region_placement_coefficients")
}
