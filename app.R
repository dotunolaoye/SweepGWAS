# SweepGWAS: the Shiny app that joins the functions in R/ (allele check, VCF
# conversion, COJO, convergence test and plots). The phenotype permutation
# in R/phenotype_permutation.R is not offered here, because each replicate
# reruns the GWAS; it is run from an R session instead.
library(shiny)
library(bslib)
library(DT)
library(shinyjs)
library(plotly)
library(data.table)

`%||%` = function(a, b) if (is.null(a)) b else a

# A rate just below 100% is shown with enough decimals not to read as 100%,
# which would contradict the dropped-marker table beneath it.
fmt_match_rate = function(n_matched, n_input) {
  if (n_input == 0) return("n/a")
  if (n_matched == n_input) return("100%")
  pct = 100 * n_matched / n_input
  digits = 1
  while (digits < 6 && round(pct, digits) >= 100) digits = digits + 1
  paste0(formatC(pct, format = "f", digits = digits), "%")
}

# Shiny runs this file from the app's own directory, so relative paths work.
for (f in c("allele_check.R", "convert_vcf_to_bed.R", "cojo_wrapper.R",
            "convergence_test.R", "region_placement.R", "plots.R", "manhattan.R", "report.R", "annotation_inputs.R")) {
  source(file.path("R", f))
}

# Shiny accepts only 5 MB uploads by default, far below a genotype panel or
# a genome-wide summary-statistic file, so the limit is raised here rather
# than in each launch script. SWEEPGWAS_MAX_UPLOAD_MB changes it.
options(shiny.maxRequestSize = as.numeric(Sys.getenv("SWEEPGWAS_MAX_UPLOAD_MB", "2000")) * 1024^2)

# PLINK 2 and GCTA are found through PLINK2_PATH and GCTA_PATH when set, then
# in bin/, where setup_binaries.sh puts them for the hosted app, and otherwise
# on the system path. The bundled copies are Linux builds, so they are used
# only on Linux. The execute permission is set again here because it is not
# always kept when the app is bundled for deployment.
find_tool = function(env_var, name) {
  set = Sys.getenv(env_var, "")
  if (nzchar(set)) return(set)
  bundled = file.path("bin", name)
  if (Sys.info()[["sysname"]] == "Linux" && file.exists(bundled)) {
    Sys.chmod(bundled, "755")
    return(normalizePath(bundled))
  }
  name
}

# The maize example in example_data/ (the genotypes and the summary
# statistics are stored with Git LFS). The genotype panel must be the one the
# GWAS was run on, filtered at MAF 0.05 (142,779 markers, matching the
# summary statistics one to one). The unfiltered 347,301-marker panel still
# runs, but the allele frequencies and COJO's LD reference then come from a
# different marker set and the signal count changes.
example_paths = list(
  gwas = file.path("example_data", "gwas_farmcpu_gdd_dts_kansas.csv"),
  bfile_prefix = file.path("example_data", "genotypes"),
  selscan = file.path("example_data", "selection_scan_regions.csv"),
  tolerance = 10000,          # Romay et al. 2013 LD decay for this panel
  same_source = TRUE          # the GWAS was run on this exact genotype panel
)

theme = bs_theme(
  version = 5,
  bg = "#FCFCFA",
  fg = "#20242B",
  primary = "#0072B2",
  secondary = "#E69F00",
  success = "#029E73",
  base_font = font_google("Inter"),
  heading_font = font_google("Source Serif 4"),
  code_font = font_google("IBM Plex Mono"),
  "border-radius" = "0.4rem",
  "body-line-height" = "1.55"
) |>
  bs_add_rules("
    .navbar { border-bottom: 1px solid #E4E1DB; }
    .sidebar { background-color: #F5F4F0 !important; }
    .card { border: 1px solid #E4E1DB; box-shadow: none; }
    .card-header { background-color: transparent; font-family: var(--bs-heading-font); font-weight: 600; }
    .btn-primary { font-weight: 500; }
    .step-label { font-size: 0.78rem; letter-spacing: 0.04em; text-transform: uppercase; color: #6B7280; font-weight: 600; margin-bottom: 0.15rem; }
    .metric-value { font-family: var(--bs-heading-font); font-size: 1.9rem; font-weight: 600; }
    .metric-label { font-size: 0.82rem; color: #6B7280; }
    .app-lede { color: #4B5157; font-size: 0.95rem; max-width: 58ch; }
    pre.diagnostic { background: #F5F4F0; border: 1px solid #E4E1DB; padding: 0.6rem 0.8rem; font-size: 0.82rem; border-radius: 0.4rem; }
    .null-check { border: 1px solid #E4E1DB; border-left: 4px solid #029E73; border-radius: 0.4rem; padding: 0.85rem 1rem; margin-bottom: 1rem; background: #FFFFFF; }
    .null-check.flag { border-left-color: #E69F00; background: #FFF9EE; }
    .null-check-head { font-family: var(--bs-heading-font); font-weight: 600; font-size: 1.02rem; margin-bottom: 0.3rem; }
    .null-check-msg { font-size: 0.9rem; color: #3A3F46; max-width: 80ch; margin-bottom: 0.6rem; }
    .null-check-table { font-family: var(--bs-code-font-family, 'IBM Plex Mono', monospace); font-size: 0.82rem; font-variant-numeric: tabular-nums; border-collapse: collapse; }
    .null-check-table th { font-weight: 600; color: #6B7280; text-align: left; padding: 0.2rem 1.2rem 0.2rem 0; border-bottom: 1px solid #E4E1DB; }
    .null-check-table td { padding: 0.2rem 1.2rem 0.2rem 0; }
    .null-check-table tr.reported td { font-weight: 600; }
    .result-verdict { font-size: 1rem; line-height: 1.55; max-width: 80ch; margin: 1rem 0 0.25rem; padding: 0.75rem 1rem; border-radius: 0.4rem; background: #F5F4F0; border: 1px solid #E4E1DB; }
    .result-verdict .caveat { display: block; margin-top: 0.35rem; font-size: 0.88rem; color: #4B5157; }
    .placement-table { font-family: var(--bs-code-font-family, 'IBM Plex Mono', monospace); font-size: 0.8rem; font-variant-numeric: tabular-nums; border-collapse: collapse; margin-top: 0.5rem; }
    .placement-table th { font-weight: 600; color: #6B7280; text-align: right; padding: 0.2rem 0 0.2rem 1.4rem; border-bottom: 1px solid #E4E1DB; }
    .placement-table td { text-align: right; padding: 0.2rem 0 0.2rem 1.4rem; }
    .placement-table th:first-child, .placement-table td:first-child { text-align: left; padding-left: 0; font-family: var(--bs-body-font-family); }
    .plot-note { font-size: 0.84rem; color: #4B5157; max-width: 85ch; margin: 0.5rem 0 0; line-height: 1.5; }
    .run-timer { position: fixed; top: 0.55rem; right: 1rem; z-index: 1050; max-width: 24rem; padding: 0.45rem 0.75rem; border: 1px solid #E4E1DB; border-radius: 0.4rem; background: #FFFFFF; font-size: 0.85rem; box-shadow: 0 2px 8px rgba(32, 36, 43, 0.08); }
    @media (max-width: 30rem) { .run-timer { left: 1rem; } }
    .run-timer.running { border-left: 4px solid #0072B2; }
    .run-timer.done { border-left: 4px solid #029E73; }
    .run-timer.stopped { border-left: 4px solid #E69F00; background: #FFF9EE; }
    .run-timer-label { color: #6B7280; margin-right: 0.35rem; }
    .run-timer-value { font-family: var(--bs-code-font-family, 'IBM Plex Mono', monospace); font-variant-numeric: tabular-nums; font-weight: 600; }
    .run-timer-stages { margin-top: 0.3rem; color: #4B5157; font-size: 0.78rem; line-height: 1.45; }
    .run-timer-stages:empty { display: none; }
    .run-timer-close { float: right; margin: -0.2rem -0.35rem 0 0.6rem; border: 0; background: transparent; color: #6B7280; font-size: 1.1rem; line-height: 1; cursor: pointer; padding: 0 0.2rem; }
    .run-timer-close:hover, .run-timer-close:focus-visible { color: #20242B; }
    .run-timer.running .run-timer-close { display: none; }
  ")

# Elapsed-time display. It runs in the browser, because R is blocked while
# GCTA or plink2 runs and a timer driven by the server would stop exactly
# when it is needed. It starts on the click and stops when the server sends
# run_timer_stop, which both run handlers send on exit, so it also stops
# when inputs are missing (it is hidden) and when a run fails.
run_timer_js = "
(function () {
  var t0 = null, iv = null;
  function el(id) { return document.getElementById(id); }
  function fmt(ms) { var s = Math.round(ms / 1000), m = Math.floor(s / 60); s = s % 60; return (m > 0 ? m + ' min ' : '') + s + ' s'; }
  function tick() { el('run_timer_value').textContent = fmt(Date.now() - t0); }
  function start() {
    t0 = Date.now();
    var box = el('run_timer');
    box.hidden = false; box.className = 'run-timer running';
    el('run_timer_label').textContent = 'Running';
    el('run_timer_stages').textContent = '';
    tick(); clearInterval(iv); iv = setInterval(tick, 1000);
  }
  document.addEventListener('click', function (e) {
    if (!e.target.closest) return;
    if (e.target.closest('#run_btn, #example_btn')) start();
    if (e.target.closest('#run_timer_close')) el('run_timer').hidden = true;
  });
  $(document).on('shown.bs.tab', function () { window.dispatchEvent(new Event('resize')); });
  $(document).on('shiny:connected', function () {
    Shiny.addCustomMessageHandler('run_timer_stop', function (msg) {
      clearInterval(iv);
      var box = el('run_timer');
      if (msg.status === 'cancelled' || t0 === null) { box.hidden = true; return; }
      el('run_timer_value').textContent = fmt(Date.now() - t0);
      el('run_timer_label').textContent = msg.status === 'done' ? 'Finished in' : 'Stopped after';
      box.className = 'run-timer ' + (msg.status === 'done' ? 'done' : 'stopped');
      el('run_timer_stages').textContent = msg.stages || '';
    });
  });
})();
"

# User interface

ui = page_sidebar(
  title = "SweepGWAS",
  theme = theme,
  useShinyjs(),
  tags$head(tags$script(HTML(run_timer_js))),
  # fixed to the top-right corner, over the empty end of the title bar
  div(id = "run_timer", class = "run-timer", hidden = NA,
    tags$button(id = "run_timer_close", class = "run-timer-close", type = "button", title = "Close", `aria-label` = "Close", HTML("&times;")),
    span(id = "run_timer_label", class = "run-timer-label", "Running"),
    span(id = "run_timer_value", class = "run-timer-value", "0 s"),
    div(id = "run_timer_stages", class = "run-timer-stages")),
  sidebar = sidebar(
    width = 340,
    div(class = "step-label", "Step 1 · GWAS summary statistics"),
    fileInput("gwas_file", NULL, accept = c(".csv", ".tsv", ".txt")),
    uiOutput("gwas_column_map_ui"),

    div(class = "step-label", "Step 2 · Genotype data"),
    radioButtons("geno_format", NULL, choices = c("PLINK bed/bim/fam" = "plink", "VCF" = "vcf"), inline = TRUE),
    conditionalPanel("input.geno_format == 'plink'",
      fileInput("bed_file", ".bed", accept = ".bed"),
      fileInput("bim_file", ".bim", accept = ".bim"),
      fileInput("fam_file", ".fam", accept = ".fam")
    ),
    conditionalPanel("input.geno_format == 'vcf'",
      fileInput("vcf_file", ".vcf / .vcf.gz", accept = c(".vcf", ".gz"))
    ),

    div(class = "step-label", "Step 3 · Selection-scan regions"),
    fileInput("selscan_file", "CSV with columns chr, start, end and method", accept = c(".csv")),

    div(class = "step-label", "Optional · Region placement inputs"),
    fileInput("genmap_file", "Genetic map: chromosome, position (bp) and genetic position (cM), as CSV or PLINK .map",
              accept = c(".csv", ".txt", ".tsv", ".map", ".gz")),
    fileInput("gene_file", "Gene annotation: GFF3 or GTF, compressed or not, or CSV of chromosome and gene start",
              accept = c(".gff", ".gff3", ".gtf", ".gz", ".csv")),
    helpText("Optional. A genetic map and a gene annotation let the region placement diagnostic test whether the ",
             "selection-scan regions sit in low-recombination or gene-rich parts of the genome, in 100 kb windows. The ",
             "convergence test itself does not use them. The recombination rate is the slope of the genetic map in cM per ",
             "Mb. Map columns are found by name (for example chr, pos and cM), or else the first three columns are taken ",
             "in that order. Gene density is counted from the gene rows of the annotation. Chromosomes may be written as ",
             "1 or chr1."),

    div(class = "step-label", "Step 4 · Parameters"),
    numericInput("tolerance", "Distance threshold in bp, from LD decay (required)", value = NA, min = 1),
    helpText("Use the linkage-disequilibrium decay distance reported for your population, or estimate it with a tool ",
             "such as PopLDdecay. There is no default. An arbitrary distance threshold can place most of the genome ",
             "within range of a region and leave the test unable to detect anything."),
    checkboxInput("use_default_cojo_p", "Use default COJO threshold (genome-wide Bonferroni)", value = TRUE),
    conditionalPanel("!input.use_default_cojo_p",
      numericInput("cojo_p", "COJO p-value threshold", value = 5e-8, min = 0, max = 1)
    ),
    numericInput("n_perm", "Number of permutations", value = 1000, min = 100, max = 10000, step = 100),
    checkboxInput("same_source", "The GWAS was run on these same genotypes", value = FALSE),
    helpText("Tick this only if the GWAS was run on these same genotype data, so that a strand flip is impossible. ",
             "Markers with A/T or C/G alleles are then matched by allele rather than by frequency. None are ",
             "dropped for having a frequency too close to 50% to resolve."),

    br(),
    actionButton("run_btn", "Run analysis", class = "btn-primary w-100"),
    actionButton("example_btn", "Load maize case-study example", class = "btn-outline-secondary w-100 mt-2"),
    br(), br(),
    p(class = "text-muted", style = "font-size: 0.78rem;",
      "Uploaded files are used only for this session and are deleted when it ends. They are not shared or sent anywhere else. ",
      "An optional check permutes phenotypes and reruns the GWAS. It runs only when the tool is installed locally.")
  ),

  navset_tab(
    nav_panel("Overview", icon = icon("circle-info"),
      div(class = "app-lede mt-3",
        p("This tool tests whether more independent GWAS signals lie near selection-scan regions than expected by",
          " chance. The GWAS markers are first checked against the genotypes for allele orientation and then resolved",
          " into independent signals by conditional and joint analysis (GCTA-COJO). The test counts signals within a",
          " distance threshold taken from linkage-disequilibrium decay, rather than an arbitrary window."),
        p("Upload your own GWAS summary statistics, genotype data and selection-scan regions. Alternatively, load the",
          " maize Ames-panel case study included with the tool to see the pipeline on real data first.")
      ),
      layout_column_wrap(
        width = "16rem",
        card(card_header("1. Allele check"), p(class = "app-lede",
          "Checks each marker's effect allele against your genotype data before any other step runs. It corrects the allele frequency where needed.")),
        card(card_header("2. COJO"), p(class = "app-lede",
          "Separates independent signals from markers that are significant only through linkage. Only signals ",
          "significant both jointly and marginally are kept. If the joint model is unstable, LD clumping is used instead.")),
        card(card_header("3. Convergence test"), p(class = "app-lede",
          "Tests whether more independent signals lie near selection-scan regions than expected by chance. Chance is ",
          "estimated by redrawing signal positions from the tested markers. A uniform null is run alongside as a ",
          "check, since it can misstate chance in either direction."))
      ),
      div(class = "plot-note mt-4", p("Contact: Dotun Olaoye, University of Missouri-Columbia. Questions and bug reports: dsolaoye@missouri.edu."))
    ),

    # The guide is the same text as the downloadable PDF, shown as HTML so it fits the
    # screen and can be searched; www/guide.html is built from the PDF's source.
    nav_panel("Guide", icon = icon("book-open"),
      div(class = "d-flex flex-wrap align-items-center gap-3 mt-3 mb-2",
        p(class = "app-lede mb-0", "How to read every step, number, plot and table in the app, including the region placement diagnostic in full."),
        downloadButton("dl_guide_pdf", "Guide (PDF)", class = "btn-sm btn-outline-primary")),
      includeHTML("www/guide.html")
    ),

    nav_panel("Allele check", icon = icon("check-double"),
      uiOutput("allele_check_summary_ui"),
      # Table cards keep their natural height (fill = FALSE), so they never scroll inside the
      # page, and the download button sits in the footer at its normal size.
      card(fill = FALSE, card_header("Dropped markers"), card_body(DTOutput("allele_check_dropped_table"), fillable = FALSE, fill = FALSE),
           card_footer(downloadButton("dl_allele_dropped", "Download dropped markers (CSV)", class = "btn-sm")))
    ),

    nav_panel("COJO signals", icon = icon("dna"),
      uiOutput("cojo_rule_ui"),
      card(card_header("Marginal vs. joint significance"), plotlyOutput("cojo_scatter_plot", height = "420px"),
           uiOutput("cojo_plot_note")),
      card(fill = FALSE, card_header("Independent signals used in the test"), card_body(DTOutput("cojo_signals_table"), fillable = FALSE, fill = FALSE),
           card_footer(class = "d-flex gap-2",
             downloadButton("dl_cojo_signals", "Download independent signals (CSV)", class = "btn-sm"),
             downloadButton("dl_cojo_output", "Download full COJO output (CSV)", class = "btn-sm")))
    ),

    nav_panel("Convergence", icon = icon("map-location-dot"),
      uiOutput("convergence_summary_ui"),
      uiOutput("null_check_ui"),
      uiOutput("placement_ui"),
      card(card_header("Manhattan plot and selection-scan regions"), uiOutput("genome_chr_ui"),
           plotlyOutput("convergence_genome_plot", height = "560px"),
           p(class = "plot-note", paste(
             "The upper panel is a Manhattan plot of every tested marker. The dashed line marks the significance",
             "threshold used for the independent signals. Signals within the distance threshold of a selection-scan",
             "region are drawn as orange diamonds. Other signals are blue circles. The lower panel shows the",
             "selection-scan regions in one track per method. A faint line runs from each signal down through the tracks.",
             "Choose a chromosome to see it alone, with coordinates in Mb. Drag across the plot to zoom in. Hover over",
             "a signal or region for its details. With very large panels, markers with p above 0.001 are thinned at random",
             "for speed. Every marker with p below 0.001 is always drawn."))),
      card(card_header("Permutation distributions under both null models"), plotlyOutput("permutation_null_plot", height = "380px"),
           p(class = "plot-note", "The blue bars show how many signals fell within the distance threshold of a region in each permutation. In these permutations, signal positions are redrawn at random from the tested markers on the same chromosome. This is the reported null. The grey outline shows the same count when positions are drawn uniformly along each chromosome instead, as a check. The observed count is marked. Each p-value is the proportion of permutations that reached it. The two distributions can overlap. How far apart their centres lie shows how much the choice of null changes the answer.")),
      card(fill = FALSE, card_header("Distance from each signal to the nearest region"), card_body(DTOutput("convergence_regions_table"), fillable = FALSE, fill = FALSE),
           card_footer(downloadButton("dl_convergence_regions", "Download convergence detail (CSV)", class = "btn-sm")))
    ),

    nav_panel("Download all", icon = icon("download"),
      p(class = "app-lede mt-3", paste("Download a report of the whole analysis as a single HTML file or as a PDF.",
                                       "The archive contains every table, a summary of the run parameters and both reports.")),
      div(class = "d-flex flex-wrap gap-2",
        downloadButton("dl_report_html", "Report (HTML)", class = "btn-primary"),
        downloadButton("dl_report_pdf", "Report (PDF)", class = "btn-primary"),
        downloadButton("dl_all_zip", "All results (.zip)", class = "btn-outline-primary")),
      p(class = "plot-note mt-3", paste("The reports give the result sentence and the rule used for the independent signals.",
        "They also give the null-model check, the region placement diagnostic and the settings. Manhattan plots of all",
        "chromosomes and of each chromosome follow, with the permutation distributions and the signal table. The run time",
        "and the software versions close the report. Its figures are static. The interactive plots are only in the app."))
    )
  )
)

# Server

server = function(input, output, session) {
  gwas_raw = reactiveVal(NULL)
  state = reactiveValues(allele_check = NULL, cojo = NULL, convergence = NULL, placement = NULL,
                          null_check = NULL, chr_lengths = NULL, selscan = NULL,
                          run_params = NULL)

  observeEvent(input$gwas_file, {
    ext = tolower(tools::file_ext(input$gwas_file$name))
    sep = if (ext == "csv") "," else "\t"
    gwas_raw(fread(input$gwas_file$datapath, sep = sep))
  })

  output$gwas_column_map_ui = renderUI({
    req(gwas_raw())
    cols = names(gwas_raw())
    guess = function(patterns) {
      hit = cols[grepl(paste(patterns, collapse = "|"), cols, ignore.case = TRUE)]
      if (length(hit) > 0) hit[1] else cols[1]
    }
    tagList(
      helpText("Choose the column in your file that holds each field. Likely matches are preselected."),
      selectInput("map_SNP", "Marker ID", cols, selected = guess(c("^snp$", "^id$", "^rsid$", "^marker"))),
      selectInput("map_CHR", "Chromosome", cols, selected = guess(c("^chr", "^chromosome"))),
      selectInput("map_POS", "Position", cols, selected = guess(c("^pos", "^bp$"))),
      selectInput("map_A1", "Effect allele", cols, selected = guess(c("^a1$", "effect.?allele", "^allele1"))),
      selectInput("map_A2", "Other allele", cols, selected = guess(c("^a2$", "other.?allele", "^allele0"))),
      selectInput("map_FREQ", "Effect-allele frequency", cols, selected = guess(c("^freq", "^maf$", "^eaf$"))),
      selectInput("map_BETA", "Effect size", cols, selected = guess(c("^beta$", "^effect$", "^b$"))),
      selectInput("map_P", "P-value", cols, selected = guess(c("^p$", "p.value", "^pval"))),
      selectInput("map_N", "Sample size", cols, selected = guess(c("^n$", "nobs", "n_obs")))
    )
  })

  build_column_map = function() {
    c(SNP = input$map_SNP, CHR = input$map_CHR, POS = input$map_POS,
      A1 = input$map_A1, A2 = input$map_A2, FREQ = input$map_FREQ,
      BETA = input$map_BETA, P = input$map_P, N = input$map_N)
  }

  work_dir = tempfile("convergence_tool_")
  dir.create(work_dir)
  # The sidebar promises that uploads are used only for this session, so the
  # session's working files and the uploaded files are deleted when it closes.
  upload_ids = c("gwas_file", "bed_file", "bim_file", "fam_file", "vcf_file", "selscan_file", "genmap_file", "gene_file")
  session$onSessionEnded(function() {
    unlink(work_dir, recursive = TRUE)
    for (id in upload_ids) {
      path = tryCatch(isolate(input[[id]]$datapath), error = function(e) NULL)
      # each upload sits in its own folder inside R's temporary directory; that directory itself is kept
      up = if (length(path)) dirname(path) else character(0)
      up = up[normalizePath(up, mustWork = FALSE) != normalizePath(tempdir(), mustWork = FALSE)]
      if (length(up)) unlink(up, recursive = TRUE)
    }
  })
  plink2_path = find_tool("PLINK2_PATH", "plink2")
  gcta_path = find_tool("GCTA_PATH", "gcta64")

  # same_source is passed in rather than read from input, because the
  # example sets the checkbox just before running and input$same_source
  # still holds the old value at that point.
  timer_stop = function(status, stages = "") {
    session$sendCustomMessage("run_timer_stop", list(status = status, stages = stages))
  }
  # stage times for the timer, e.g. "Allele check 2 s · COJO 1 min 28 s";
  # stages that did not run are left out
  format_stages = function(timing) {
    timing = timing[!is.na(timing)]
    fmt = function(x) { x = round(x); if (x >= 60) sprintf("%d min %d s", x %/% 60, x %% 60) else sprintf("%d s", x) }
    labels = c(prep = "Data preparation", allele_check = "Allele check", cojo = "COJO",
               clumping = "LD clumping", convergence = "Convergence test",
               placement = "Region placement")
    paste(paste(labels[names(timing)], vapply(timing, fmt, "")), collapse = " \u00b7 ")
  }

  run_pipeline = function(gwas_dt, column_map, bfile_prefix, selscan_dt, tolerance, cojo_p, n_perm,
                           same_source = isTRUE(input$same_source), prep_s = NA_real_) {
    timing = c(prep = prep_s)
    withProgress(message = "Running analysis", value = 0, {
      incProgress(0.1, detail = "Allele check + COJO")
      cojo_res = run_cojo(gwas_dt, column_map, bfile_prefix, plink2_path, gcta_path, work_dir,
                           cojo_p = cojo_p, same_source = same_source)
      state$allele_check = cojo_res$allele_check
      state$cojo = cojo_res

      timing = c(timing, cojo_res$timing)
      if (cojo_res$status == "no_signals") {
        showNotification(cojo_res$message, type = "error", duration = NULL)
        return(list(status = "stopped", stages = format_stages(timing)))
      }
      if (cojo_res$status == "few_signals") showNotification(cojo_res$message, type = "warning", duration = NULL)
      if (cojo_res$method == "clumping") {
        showNotification(paste0("The COJO joint model was unstable on these data. Independent loci were taken ",
                                "by LD clumping instead. See the COJO signals tab."), type = "warning", duration = 12)
      }

      incProgress(0.5, detail = "Convergence permutation test (both null models)")
      gwas_std = data.table(CHR = std_chr(gwas_dt[[column_map["CHR"]]]), POS = gwas_dt[[column_map["POS"]]])
      gwas_chrs = unique(gwas_std$CHR)
      example_chrs = function(x) paste(head(unique(x), 3), collapse = ", ")
      # A regions file whose chromosome names do not match the GWAS would
      # otherwise give no convergence at all without any warning.
      selscan_dt = copy(as.data.table(selscan_dt))[, chr := std_chr(chr)]
      unmatched = setdiff(unique(selscan_dt$chr), gwas_chrs)
      if (length(unmatched) == uniqueN(selscan_dt$chr)) {
        stop(sprintf(paste0("none of the chromosome names in the selection-scan regions file (for example %s) match those ",
                            "in the GWAS file (for example %s). Check that both use the same chromosome names and genome assembly."),
                     example_chrs(selscan_dt$chr), example_chrs(gwas_chrs)))
      }
      if (length(unmatched)) {
        showNotification(sprintf("Selection-scan regions on %d chromosome(s) absent from the GWAS file were ignored: %s.",
                                 length(unmatched), paste(head(unmatched, 10), collapse = ", ")), type = "warning", duration = 12)
      }
      signals = copy(cojo_res$signals)[, Chr := std_chr(Chr)]
      # every tested marker's p-value, for the Manhattan plot and the reports
      state$markers = data.table(chr = gwas_std$CHR, pos = gwas_std$POS, p = as.numeric(gwas_dt[[column_map["P"]]]))
      chr_lengths = derive_chr_lengths(gwas_std)
      # The reported result draws permuted positions from the tested
      # markers. The uniform null is run beside it only as a check, since a
      # disagreement between the two would otherwise go unnoticed.
      marker_positions = gwas_std[, .(chr = CHR, pos = POS)]
      t_conv = Sys.time()
      null_check = compare_null_models(signals, selscan_dt, tolerance,
                                       marker_positions, n_perm = n_perm)
      timing["convergence"] = as.numeric(difftime(Sys.time(), t_conv, units = "secs"))
      conv_res = null_check$markers

      incProgress(0.3, detail = "Region placement diagnostic")
      t_place = Sys.time()
      # The genetic map and the gene annotation are optional. A file that cannot
      # be read, or whose chromosomes match none in the GWAS, is left out with a
      # warning rather than stopping the run or adding a feature of constants.
      read_optional = function(inp, what, reader) {
        if (is.null(inp)) return(NULL)
        d = tryCatch(reader(inp$datapath, inp$name), error = function(e) {
          showNotification(sprintf("The %s could not be read (%s). It was not used.", what, conditionMessage(e)),
                           type = "warning", duration = 15)
          NULL
        })
        if (is.null(d)) return(NULL)
        if (!any(d$chr %in% gwas_chrs)) {
          showNotification(sprintf(paste0("None of the chromosome names in the %s (for example %s) match the GWAS file ",
                                          "(for example %s). It was not used."), what, example_chrs(d$chr), example_chrs(gwas_chrs)),
                           type = "warning", duration = 15)
          return(NULL)
        }
        d
      }
      placement = tryCatch(
        fit_region_placement(
          data.table(chr = gwas_std$CHR, pos = gwas_dt[[column_map["POS"]]],
                     freq = gwas_dt[[column_map["FREQ"]]], beta = gwas_dt[[column_map["BETA"]]],
                     p = gwas_dt[[column_map["P"]]]),
          selscan_dt, tolerance, genetic_map = read_optional(input$genmap_file, "genetic map", read_genetic_map),
          genes = read_optional(input$gene_file, "gene annotation", read_gene_annotation)),
        error = function(e) list(message = paste("The region placement diagnostic could not be fitted:", conditionMessage(e)),
                                 coefficients = data.table(), feature_means = data.table(), auc = NA_real_,
                                 avg_precision = NA_real_, ap_baseline = NA_real_))
      timing["placement"] = as.numeric(difftime(Sys.time(), t_place, units = "secs"))
      state$placement = placement
      state$convergence = conv_res
      state$null_check = null_check
      state$chr_lengths = chr_lengths
      state$selscan = selscan_dt
      state$run_params = data.table(
        parameter = c("distance_threshold_bp", "cojo_p_threshold", "n_perm", "n_input_markers",
                      "n_markers_after_allele_check", "same_source", "independent_signal_method",
                      "n_cojo_signals_selected", "n_cojo_joint_and_marginal", "cojo_median_joint_marginal_ratio",
                      "cojo_max_median_ratio", "n_independent_signals",
                      "null_model_reported", "null_models_agree", "null_baseline_ratio",
                      "placement_auc_by_chromosome", "placement_average_precision", "placement_precision_baseline",
                      "placement_features"),
        value = c(tolerance, cojo_res$threshold_used, n_perm, nrow(gwas_dt),
                  nrow(cojo_res$allele_check$matched), same_source, cojo_res$method,
                  cojo_res$rules$n_cojo, cojo_res$rules$n_retained, round(cojo_res$rules$median_ratio, 2),
                  cojo_res$rules$max_median_ratio, cojo_res$n_signals,
                  "markers", !null_check$conclusions_differ,
                  round(null_check$baseline_ratio, 2),
                  round(placement$auc, 3), round(placement$avg_precision, 3), round(placement$ap_baseline, 3),
                  paste(placement$features, collapse = " "))
      )
      state$run_params = rbind(state$run_params, data.table(
        parameter = paste0("seconds_", names(timing)[!is.na(timing)]),
        value = as.character(round(timing[!is.na(timing)], 1))))
      incProgress(1)
    })
    state$stages = format_stages(timing)
    list(status = "done", stages = state$stages)
  }

  # Runs the pipeline so that an error appears as a notification rather than
  # ending the session, and always returns a status for the timer.
  safe_run = function(...) {
    tryCatch(run_pipeline(...), error = function(e) {
      showNotification(paste("The analysis stopped:", conditionMessage(e)), type = "error", duration = NULL)
      list(status = "error", stages = "")
    })
  }

  observeEvent(input$run_btn, {
    out = list(status = "cancelled", stages = "")
    on.exit(timer_stop(out$status, out$stages))
    req(gwas_raw(), input$selscan_file)
    if (is.na(input$tolerance)) { showNotification("Enter the distance threshold in bp, from LD decay, before running the analysis.", type = "error"); return() }
    t_prep = Sys.time()
    validate_geno = if (input$geno_format == "plink") {
      req(input$bed_file, input$bim_file, input$fam_file)
      prefix = file.path(work_dir, "geno_upload")
      file.copy(input$bed_file$datapath, paste0(prefix, ".bed"), overwrite = TRUE)
      file.copy(input$bim_file$datapath, paste0(prefix, ".bim"), overwrite = TRUE)
      file.copy(input$fam_file$datapath, paste0(prefix, ".fam"), overwrite = TRUE)
      prefix
    } else {
      req(input$vcf_file)
      conv = convert_vcf_to_bed(input$vcf_file$datapath, plink2_path, work_dir, out_prefix = "geno_upload")
      conv$bfile_prefix
    }
    selscan_dt = fread(input$selscan_file$datapath)
    cojo_p = if (input$use_default_cojo_p) NULL else input$cojo_p
    out = safe_run(gwas_raw(), build_column_map(), validate_geno, selscan_dt, input$tolerance, cojo_p, input$n_perm,
                   prep_s = as.numeric(difftime(Sys.time(), t_prep, units = "secs")))
  })

  observeEvent(input$example_btn, {
    out = list(status = "error", stages = "")
    on.exit(timer_stop(out$status, out$stages))
    t_prep = Sys.time()
    # run_cojo() needs the genotype files inside work_dir, as for uploads.
    ex_prefix = file.path(work_dir, "example_geno")
    for (ext in c(".bed", ".bim", ".fam")) {
      file.copy(paste0(example_paths$bfile_prefix, ext), paste0(ex_prefix, ext), overwrite = TRUE)
    }

    # GAPIT's output has no allele columns, only the minor-allele frequency,
    # so the effect allele and its frequency are reconstructed here.
    # The effect allele is the one plink2 --export A counted when the GWAS
    # dosages were built, which the .raw headers show was the .bim A2 for all
    # 142,779 markers, against plink2's documented default, because REF is
    # provisional in these data. Taking the .bim A1 instead would reverse
    # every effect.
    # The frequency must be that allele's, not the minor-allele frequency,
    # which matches it for only about 79% of markers here, so it is computed
    # from the genotypes with plink2 --freq.
    ex_gwas = fread(example_paths$gwas)
    bim = fread(paste0(example_paths$bfile_prefix, ".bim"),
                col.names = c("bim_chr", "SNP", "cm", "bim_pos", "bim_A1", "bim_A2"))

    freq_log = run_plink2_relative(plink2_path,
      c("--bfile", basename(ex_prefix), "--freq", "--out", "example_freq", "--allow-extra-chr"), work_dir)
    afreq_path = file.path(work_dir, "example_freq.afreq")
    if (!file.exists(afreq_path)) {
      showNotification(paste("plink2 --freq failed for the example dataset:",
                             paste(utils::tail(freq_log, 3), collapse = " ")), type = "error", duration = NULL)
      return(invisible())
    }
    afreq = fread(afreq_path)
    setnames(afreq, c("#CHROM", "ID", "REF", "ALT", "PROVISIONAL_REF?", "ALT_FREQS", "OBS_CT"),
             c("chr_f", "SNP", "REF", "ALT", "prov_ref", "alt_freq", "obs_ct"), skip_absent = TRUE)

    ex_gwas = merge(ex_gwas, bim[, .(SNP, bim_A1, bim_A2)], by = "SNP")
    ex_gwas = merge(ex_gwas, afreq[, .(SNP, ALT, alt_freq)], by = "SNP")
    ex_gwas[, A1 := bim_A2]  # the counted/effect allele (verified above)
    ex_gwas[, A2 := bim_A1]
    ex_gwas[, EAF := fifelse(A1 == ALT, alt_freq, 1 - alt_freq)]
    # merge() sorts by marker ID, and the permutations draw markers by row,
    # so genome order is restored to reproduce the paper's result exactly
    setorder(ex_gwas, Chr, Pos)

    gwas_raw(ex_gwas)
    updateNumericInput(session, "tolerance", value = example_paths$tolerance)
    # The GWAS was run on this genotype panel, so no strand flip is possible
    # and A/T and C/G markers need not be resolved by frequency.
    updateCheckboxInput(session, "same_source", value = example_paths$same_source)
    selscan_dt = fread(example_paths$selscan)
    col_map = c(SNP = "SNP", CHR = "Chr", POS = "Pos", A1 = "A1", A2 = "A2",
                FREQ = "EAF", BETA = "effect", P = "P.value", N = "nobs")

    showNotification("The maize example is loaded and running. This takes about two to four minutes.", type = "message")
    # passed explicitly, since input$same_source still holds its old value
    out = safe_run(ex_gwas, col_map, ex_prefix, selscan_dt,
                   example_paths$tolerance, NULL, input$n_perm %||% 1000,
                   same_source = example_paths$same_source,
                   prep_s = as.numeric(difftime(Sys.time(), t_prep, units = "secs")))
  })

  output$allele_check_summary_ui = renderUI({
    req(state$allele_check)
    s = state$allele_check$summary
    n_fmt = function(x) format(x, big.mark = ",")
    # Every drop reason is listed, so the total always equals input minus passing markers.
    breakdown = sprintf("Dropped by reason. Not SNPs (insertions, deletions or missing-allele codes): %s. Not found in the genotype data: %s. Alleles that do not match the genotypes: %s. Allele frequency that disagrees with the genotypes: %s. A/T or C/G markers too close to 50%% frequency to resolve: %s. Each dropped marker is listed in the table below.",
                        n_fmt(s$n_dropped_non_snp), n_fmt(s$n_dropped_not_in_genotypes), n_fmt(s$n_dropped_no_match),
                        n_fmt(s$n_dropped_freq_inconsistent), n_fmt(s$n_dropped_ambiguous_unresolvable))
    tagList(
      layout_column_wrap(width = "10rem", class = "mt-3 mb-2",
        div(card(div(class = "metric-value", n_fmt(s$n_input)), div(class = "metric-label", "Input markers"))),
        div(card(div(class = "metric-value", n_fmt(s$n_matched)), div(class = "metric-label", "Markers passing the check"))),
        div(card(div(class = "metric-value", fmt_match_rate(s$n_matched, s$n_input)), div(class = "metric-label", "Match rate"))),
        div(card(div(class = "metric-value", n_fmt(s$n_dropped_total)), div(class = "metric-label", "Markers dropped (all reasons)"))),
        div(card(div(class = "metric-value", n_fmt(s$n_ambiguous_resolved_by_freq)), div(class = "metric-label", "A/T and C/G markers resolved"))),
        div(card(div(class = "metric-value", n_fmt(s$n_dropped_ambiguous_unresolvable)), div(class = "metric-label", "A/T and C/G markers dropped as unresolvable"))),
        div(card(div(class = "metric-value", n_fmt(s$n_freq_swap_corrected)), div(class = "metric-label", "Frequencies corrected")))
      ),
      p(class = "plot-note mb-3", breakdown)
    )
  })
  output$allele_check_dropped_table = renderDT({ req(state$allele_check); state$allele_check$dropped }, options = list(pageLength = 10))
  output$dl_allele_dropped = downloadHandler("dropped_markers.csv", function(f) fwrite(state$allele_check$dropped, f))

  output$cojo_scatter_plot = renderPlotly({
    req(state$cojo, state$cojo$cojo_output)
    # GCTA's full output is plotted, so the removed joint-only signals stay
    # visible, with tested signals coloured by the convergence result
    plot_cojo_scatter(state$cojo$cojo_output, state$convergence$regions, method = state$cojo$method)
  })
  output$cojo_plot_note = renderUI({
    req(state$cojo, state$cojo$cojo_output)
    txt = if (identical(state$cojo$method, "clumping")) {
      paste("Each point is a signal selected by COJO, placed by its significance on its own (horizontal axis) and in the",
            "joint model (vertical axis). Because the joint model was unstable, none of these signals was used in the test,",
            "which was run on the LD-clumped loci listed below instead.")
    } else {
      paste("Each point is a signal selected by COJO, placed by its significance on its own (horizontal axis) and in the",
            "joint model (vertical axis). Points near the dotted diagonal kept their significance after conditioning on",
            "the other signals. Grey open circles were significant only jointly and were removed. Filled points are the",
            "signals used in the test, coloured by whether they lie within the distance threshold of a selection-scan region.")
    }
    p(class = "plot-note", txt)
  })
  output$dl_cojo_output = downloadHandler("cojo_full_output.csv", function(f) fwrite(state$cojo$cojo_output, f))
  # which rule produced the independent signals, with its numbers
  output$cojo_rule_ui = renderUI({
    req(state$cojo, state$cojo$rules)
    cj = state$cojo; r = cj$rules
    ratio_txt = if (is.na(r$median_ratio)) "not defined" else sprintf("%.2f", r$median_ratio)
    if (cj$method == "cojo") {
      head = "Independent signals from COJO"
      msg = sprintf(paste0("COJO selected %d signals, of which %d were significant both jointly and marginally and ",
                           "were kept. The other %d were significant only jointly and were removed. Such signals can arise ",
                           "when correlated markers are given large offsetting effects. The joint model looks stable. ",
                           "Among the kept signals, the median joint effect is %s times the marginal effect. The ",
                           "model is treated as unstable only when this ratio exceeds %s. Independent loci are then ",
                           "taken by LD clumping instead."),
                    r$n_cojo, r$n_retained, r$n_joint_only, ratio_txt, r$max_median_ratio)
    } else {
      head = "COJO joint model unstable: independent loci taken by LD clumping"
      why = if (r$n_retained == 0) "No signal was significant both jointly and marginally." else
        sprintf(paste0("Among the %d signals significant both jointly and marginally, the median joint effect is %s times ",
                       "the marginal effect. This is above the stability limit of %s."),
                r$n_retained, ratio_txt, r$max_median_ratio)
      msg = sprintf(paste0("COJO selected %d signals. %s This can indicate collinear markers rather than ",
                           "separate associations. Collinear markers are common with long-range LD, strong population ",
                           "structure or selfing. Independent loci were therefore taken by LD clumping, which fits no joint model ",
                           "(index p < %s, r\u00b2 %s within %s kb). Clumping gave %d loci."),
                    r$n_cojo, why, sub("e([-+])0*", "e\\1", format(cj$threshold_used, digits = 3, scientific = TRUE)),
                    cj$clump_settings[["r2"]], cj$clump_settings[["kb"]], cj$n_signals)
    }
    div(class = paste("null-check", if (cj$method == "clumping") "flag"),
      div(class = "null-check-head", head),
      p(class = "null-check-msg", msg))
  })
  output$cojo_signals_table = renderDT({ req(state$cojo); state$cojo$signals }, options = list(pageLength = 10))
  output$dl_cojo_signals = downloadHandler("cojo_signals.csv", function(f) fwrite(state$cojo$signals, f))

  output$convergence_summary_ui = renderUI({
    req(state$convergence)
    cv = state$convergence
    v = convergence_verdict(cv, state$null_check)
    tagList(
      div(class = "result-verdict", v$text,
          if (nzchar(v$caveat)) span(class = "caveat", v$caveat)),
      layout_column_wrap(width = "10rem", class = "mt-3 mb-3",
      div(card(div(class = "metric-value", paste0(cv$observed, " / ", cv$n_regions)), div(class = "metric-label", "Convergent signals"))),
      div(card(div(class = "metric-value", paste0(cv$observed_pct, "%")), div(class = "metric-label", "Observed share"))),
      div(card(div(class = "metric-value", paste0(cv$null_mean_pct, "%")), div(class = "metric-label", "Expected by chance"))),
      div(card(div(class = "metric-value", format(cv$empirical_p, digits = 3)), div(class = "metric-label", "p-value (marker-based null)")))
    ))
  })
  # Null-model check: shows the uniform null beside the reported
  # marker-based result and flags when their conclusions differ.
  output$null_check_ui = renderUI({
    req(state$null_check)
    nc = state$null_check
    head = switch(nc$status,
      uniform_only = "Null models disagree: the uniform null would overstate convergence",
      markers_only = "Null models disagree: the uniform null would understate convergence",
      agree        = "Null models agree")
    rows = lapply(seq_len(nrow(nc$table)), function(i) {
      r = nc$table[i]
      tags$tr(class = if (i == 1) "reported" else NULL,
        tags$td(if (i == 1) "Drawn from tested markers (reported)" else "Drawn uniformly along each chromosome (check)"),
        tags$td(sprintf("%.1f%%", r$null_mean_pct)),
        tags$td(format(signif(r$empirical_p, 3), scientific = FALSE)))
    })
    div(class = paste("null-check", if (nc$conclusions_differ) "flag"),
      div(class = "null-check-head", head),
      p(class = "null-check-msg", nc$message),
      tags$table(class = "null-check-table",
        tags$thead(tags$tr(tags$th("Null model"), tags$th("Expected by chance"), tags$th("p-value"))),
        tags$tbody(rows)))
  })
  # Region placement diagnostic: where the selection-scan regions sit relative
  # to the markers, which explains a disagreement between the two nulls.
  output$placement_ui = renderUI({
    req(state$placement)
    pl = state$placement
    has_fit = nrow(pl$coefficients) > 0
    rows = if (has_fit) lapply(seq_len(nrow(pl$feature_means)), function(i) {
      r = pl$feature_means[i]
      fmt = function(x) formatC(x, digits = 3, format = "fg", big.mark = ",")
      tags$tr(tags$td(r$label), tags$td(fmt(r$mean_near)), tags$td(fmt(r$mean_not_near)),
              tags$td(sprintf("%.2f", r$ratio)))
    })
    stats = if (has_fit && !is.na(pl$auc)) sprintf(paste0(
      "Average precision %.3f, which is %.1f times its baseline of %.3f. The baseline is the share of markers near a region."),
      pl$avg_precision, pl$avg_precision / pl$ap_baseline, pl$ap_baseline)
    card(card_header("Region placement diagnostic"),
      p(class = "null-check-msg", pl$message),
      if (!is.null(stats)) p(class = "plot-note", style = "margin-top: 0;", tags$b(stats)),
      if (has_fit) plotlyOutput("placement_plot", height = paste0(90 + 42 * nrow(pl$coefficients), "px")),
      if (has_fit) tags$table(class = "placement-table",
        tags$thead(tags$tr(tags$th("Feature"), tags$th("Near a region"), tags$th("Elsewhere"), tags$th("Ratio"))),
        tags$tbody(rows)),
      if (has_fit) p(class = "plot-note", paste(
        "The model is a logistic regression of whether each tested marker lies within the distance threshold of a",
        "selection-scan region.",
        "Bars show standardised coefficients. Blue bars mean the feature is higher near regions, and grey bars mean lower.",
        "Longer bars matter more. The lines are 95% intervals from a jackknife over chromosomes, which leaves",
        "out one chromosome at a time. A line that crosses zero means the direction of that feature is not",
        "consistent across chromosomes. The table gives the average of each feature for markers near a region",
        "and elsewhere. The AUC and average precision are measured on each chromosome while it is withheld.",
        "Neighbouring markers share their window features, so a random split would overstate them.",
        sprintf("The model was refitted %d times, once without each chromosome.", pl$n_chr))),
      if (has_fit) p(class = "plot-note", "See section 6 of the Guide tab for how to read this diagnostic."))
  })
  output$placement_plot = renderPlotly({
    req(state$placement, nrow(state$placement$coefficients) > 0)
    plot_placement_coefficients(state$placement$coefficients)
  })
  output$genome_chr_ui = renderUI({
    req(state$convergence, state$markers)
    chrs = chr_order(state$markers$chr)
    selectInput("genome_chr", NULL, choices = c("All chromosomes" = "all", setNames(chrs, paste("Chromosome", chrs))),
                selected = "all", width = "220px")
  })
  output$convergence_genome_plot = renderPlotly({
    req(state$convergence, state$selscan, state$markers)
    this_chr = if (is.null(input$genome_chr) || input$genome_chr == "all") NULL else input$genome_chr
    plot_convergence_manhattan(state$markers, state$convergence$regions, state$selscan,
                               threshold = state$cojo$threshold_used, chr = this_chr)
  })
  output$permutation_null_plot = renderPlotly({
    req(state$convergence)
    cv = state$convergence
    nc = state$null_check
    plot_permutation_null(cv$null_distribution, cv$observed, cv$empirical_p,
                           uniform_null = nc$sequence$null_distribution, uniform_p = nc$sequence$empirical_p,
                           x_label = paste0("# convergent signals (of ", cv$n_regions, ")"))
  })
  output$convergence_regions_table = renderDT({ req(state$convergence); state$convergence$regions }, options = list(pageLength = 10))
  output$dl_convergence_regions = downloadHandler("convergence_regions.csv", function(f) fwrite(state$convergence$regions, f))

  output$dl_guide_pdf = downloadHandler("SweepGWAS_guide.pdf", function(f) file.copy("www/SweepGWAS_guide.pdf", f))

  report_ready = function() req(state$cojo, state$convergence, state$markers)
  output$dl_report_html = downloadHandler("sweepgwas_report.html", function(f) {
    report_ready()
    withProgress(message = "Building the HTML report", value = 0.5, write_report_html(f, state, plink2_path, gcta_path))
  })
  output$dl_report_pdf = downloadHandler("sweepgwas_report.pdf", function(f) {
    report_ready()
    withProgress(message = "Building the PDF report", value = 0.5, write_report_pdf(f, state, plink2_path, gcta_path))
  })

  output$dl_all_zip = downloadHandler(
    filename = "convergence_results.zip",
    content = function(zipfile) {
      req(state$cojo, state$convergence)
      tmp = tempfile(); dir.create(tmp)
      fwrite(state$allele_check$dropped, file.path(tmp, "dropped_markers.csv"))
      fwrite(state$allele_check$matched, file.path(tmp, "matched_markers.csv"))
      fwrite(state$cojo$signals, file.path(tmp, "independent_signals.csv"))
      fwrite(state$cojo$cojo_output, file.path(tmp, "cojo_full_output.csv"))
      fwrite(state$convergence$regions, file.path(tmp, "convergence_regions.csv"))
      fwrite(data.table(null_distribution = state$convergence$null_distribution), file.path(tmp, "permutation_null.csv"))
      if (!is.null(state$null_check)) fwrite(state$null_check$table, file.path(tmp, "null_model_comparison.csv"))
      if (!is.null(state$placement) && nrow(state$placement$coefficients) > 0) {
        fwrite(state$placement$coefficients, file.path(tmp, "region_placement_coefficients.csv"))
        fwrite(state$placement$feature_means, file.path(tmp, "region_placement_feature_means.csv"))
        fwrite(state$placement$auc_by_chr, file.path(tmp, "region_placement_auc_by_chromosome.csv"))
      }
      fwrite(state$run_params, file.path(tmp, "run_parameters.csv"))
      withProgress(message = "Building the reports", value = 0.5, {
        write_report_html(file.path(tmp, "sweepgwas_report.html"), state, plink2_path, gcta_path)
        write_report_pdf(file.path(tmp, "sweepgwas_report.pdf"), state, plink2_path, gcta_path)
      })
      old = setwd(tmp); on.exit(setwd(old))
      zip(zipfile, files = list.files(tmp))
    }
  )
}

shinyApp(ui, server)
