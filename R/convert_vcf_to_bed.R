# Converts a user's VCF to PLINK bed/bim/fam, since run_cojo() accepts only
# that format. Uses run_plink2_relative() from cojo_wrapper.R, which must be
# sourced first.
#
# Multiallelic sites are dropped with --max-alleles 2, because plink2 refuses
# to write them to a .bim file, and they are common enough in real VCFs that
# stopping on them would block most users.
library(data.table)

#' Convert a VCF to PLINK bed/bim/fam, ready for run_cojo()'s bfile_prefix.
#'
#' @param vcf_path path to the input VCF (.vcf or .vcf.gz). It can be
#'   anywhere, since it is copied into work_dir under a fixed name before
#'   plink2 is called.
#' @param plink2_path path to the plink2 binary.
#' @param work_dir a writable directory where plink2 runs and the bed/bim/fam
#'   files are written; normally the same work_dir later passed to run_cojo().
#' @param out_prefix output file name without extension or path.
#'
#' @return a list: $bfile_prefix (full path, ready for run_cojo()),
#'   $n_samples, $n_variants, and $log (plink2's output, for diagnosing a
#'   failed conversion).
convert_vcf_to_bed = function(vcf_path, plink2_path, work_dir, out_prefix = "genotypes") {
  if (!file.exists(vcf_path)) stop("VCF not found: ", vcf_path)
  if (!dir.exists(work_dir)) dir.create(work_dir, recursive = TRUE)

  vcf_local_name = paste0("input_vcf", if (grepl("\\.gz$", vcf_path)) ".vcf.gz" else ".vcf")
  vcf_local_path = file.path(work_dir, vcf_local_name)
  if (normalizePath(vcf_path, mustWork = FALSE) != normalizePath(vcf_local_path, mustWork = FALSE)) {
    file.copy(vcf_path, vcf_local_path, overwrite = TRUE)
  }

  cat("Converting VCF to PLINK bed format (biallelic sites only)...\n")
  log = run_plink2_relative(plink2_path,
    c("--vcf", vcf_local_name, "--max-alleles", "2", "--make-bed",
      "--out", out_prefix, "--allow-extra-chr"),
    work_dir)

  bfile_prefix = file.path(work_dir, out_prefix)
  bim_path = paste0(bfile_prefix, ".bim")
  fam_path = paste0(bfile_prefix, ".fam")
  if (!file.exists(bim_path) || !file.exists(fam_path)) {
    stop("VCF-to-bed conversion failed, no .bim/.fam produced. plink2 output:\n",
         paste(log, collapse = "\n"))
  }

  n_samples = nrow(fread(fam_path, header = FALSE))
  n_variants = nrow(fread(bim_path, header = FALSE))
  log_line = grep("variants? (loaded|remaining)|out of", log, ignore.case = TRUE, value = TRUE)
  cat("Conversion done:", n_samples, "samples,", n_variants, "biallelic variants kept.\n")
  if (length(log_line) > 0) cat(paste(log_line, collapse = "\n"), "\n")

  list(bfile_prefix = bfile_prefix, n_samples = n_samples, n_variants = n_variants,
       log = log)
}
