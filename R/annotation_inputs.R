# Readers for the two optional inputs of the region placement diagnostic: a
# genetic map, from which the local recombination rate is taken, and a gene
# annotation, from which gene density is taken. Both return chromosome names
# passed through std_chr(), so that they can be matched to the GWAS.
library(data.table)

#' Read a genetic map: a CSV or tab-separated file with chromosome, physical
#' position (bp) and genetic position (cM), or a PLINK .map file (chromosome,
#' marker, cM, bp, without a header).
#'
#' Columns are found by name when the file has a header, and otherwise taken
#' as the first three in the order chromosome, bp, cM.
#' @param path the file.
#' @param name the original file name, whose extension identifies a .map file.
#' @return a data.table with chr, bp and cm.
read_genetic_map = function(path, name = path) {
  if (grepl("\\.map(\\.gz)?$", name, ignore.case = TRUE)) {
    d = fread(path, header = FALSE)
    if (ncol(d) < 4) stop("a PLINK .map file needs four columns: chromosome, marker, cM and bp")
    d = d[, .(chr = V1, bp = V4, cm = V3)]
  } else {
    d = fread(path)
    nm = tolower(names(d))
    pick = function(pat) { i = which(grepl(pat, nm)); if (length(i)) i[1] else NA_integer_ }
    ic = pick("^(chr|chrom|chromosome|seqid)$")
    ib = pick("^(pos|position|bp|physical|phys_pos|physical_position|base_pair|bp_pos)$")
    im = pick("(^cm$|centimorgan|genetic|^gen_pos$|^map_pos$|cm_pos|^pos_cm$)")
    if (any(is.na(c(ic, ib, im)))) {
      if (ncol(d) < 3) stop("the genetic map needs three columns: chromosome, position in bp and genetic position in cM")
      ic = 1L; ib = 2L; im = 3L
    }
    d = data.table(chr = d[[ic]], bp = d[[ib]], cm = d[[im]])
  }
  d = d[, .(chr = std_chr(chr), bp = suppressWarnings(as.numeric(bp)), cm = suppressWarnings(as.numeric(cm)))]
  d = d[!is.na(bp) & !is.na(cm)]
  if (nrow(d) < 2) stop("the genetic map has fewer than two usable positions")
  d
}

#' Read gene start positions from a GFF3 or GTF annotation, gzipped or not, or
#' from a CSV whose first two columns are chromosome and gene start.
#'
#' Only rows whose feature type is "gene" are used, one per gene. Some GTF
#' files have no gene rows, and then transcript or mRNA rows are used, with
#' repeated starts counted once. The file is read in chunks, and a FASTA
#' section at the end of a GFF3 is skipped.
#' @return a data.table with chr and pos (the gene start).
read_gene_annotation = function(path, name = path) {
  if (!grepl("\\.(gff3?|gtf)(\\.gz)?$", name, ignore.case = TRUE)) {
    d = fread(path)
    if (ncol(d) < 2) stop("the gene file needs two columns: chromosome and gene start")
    d = data.table(chr = std_chr(d[[1]]), pos = suppressWarnings(as.numeric(d[[2]])))
    return(d[!is.na(pos)])
  }
  con = gzfile(path, "r")   # reads plain and gzipped files alike
  on.exit(close(con))
  out = list()
  repeat {
    x = readLines(con, n = 200000, warn = FALSE)
    if (!length(x)) break
    fasta = which(startsWith(x, "##FASTA"))
    if (length(fasta)) x = x[seq_len(fasta[1] - 1)]
    x = x[nzchar(x) & !startsWith(x, "#")]
    f = strsplit(x, "\t", fixed = TRUE)
    f = f[lengths(f) >= 5]
    if (length(f)) {
      type = vapply(f, `[`, "", 3)
      keep = type %in% c("gene", "transcript", "mRNA")
      out[[length(out) + 1]] = data.table(chr = vapply(f[keep], `[`, "", 1), type = type[keep],
                                          pos = as.numeric(vapply(f[keep], `[`, "", 4)))
    }
    if (length(fasta)) break
  }
  a = rbindlist(out)
  if (!nrow(a)) stop("no gene, transcript or mRNA rows were found in the annotation")
  a = if (any(a$type == "gene")) a[type == "gene"] else unique(a[, .(chr, pos)])
  a[, .(chr = std_chr(chr), pos)][!is.na(pos)]
}
