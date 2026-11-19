# Allele-orientation check, run on every dataset before COJO. GCTA-COJO does
# not check that the summary-statistic effect allele matches the genotype
# coding, and a mismatch gives a confidently wrong result with nothing in the
# output to show it. The maize data in this project had exactly this problem:
# plink2 had counted the REF allele rather than A1 when the GWAS dosages were
# built.
#
# Allele identity alone is not enough, because a supplied frequency can
# belong to the wrong allele even when the alleles match. A MAF column, for
# example, gives the minor allele's frequency, which is not the effect
# allele's frequency whenever the effect allele is the major one. The
# supplied frequency is therefore compared with one computed from the
# genotypes.
#
# The check handles SNPs only. complement_allele() swaps letters without
# reversing the string, which is correct for single-letter alleles but wrong
# for insertions and deletions, so those are dropped rather than
# miscomputed.
library(data.table)

complement_allele = function(a) {
  chartr("ACGT", "TGCA", a)
}

#' Verify and orient GWAS summary-statistic alleles against genotype data.
#'
#' @param gwas_dt data.table with columns SNP, A1 (effect allele), A2 (other
#'   allele), FREQ (supplied frequency of A1), plus the other GWAS columns
#'   (BETA, SE, P, N), already renamed to these standard names.
#' @param bim_path path to the genotype .bim file (chr, SNP, cm, pos, A1, A2).
#' @param freq_path path to a plink2 --freq output (.afreq) computed from the
#'   same genotypes, used as the reference for the supplied FREQ.
#' @param freq_tolerance largest allowed absolute difference between the
#'   supplied and the genotype-derived frequency. The default of 0.05 allows
#'   for sampling differences when the GWAS sample and the genotype panel are
#'   not the same individuals, while still catching a frequency reported for
#'   the other allele.
#' @param ambiguous_freq_margin how close to 0.5 the frequency of an A/T or
#'   C/G SNP must be before its strand is treated as unresolvable and the
#'   marker dropped. Default 0.05.
#' @param same_source TRUE only when the GWAS was run on these same genotypes.
#'   A strand flip is then impossible, so A/T and C/G SNPs are matched by
#'   allele and none are dropped for a frequency near 0.5. Only the user can
#'   know this, so it defaults to FALSE; setting it wrongly brings back the
#'   miscoding risk this check exists to prevent.
#'
#' @return a list: $matched (data.table ready for the .ma file, with the
#'   genotype-derived FREQ), $summary (counts by outcome), and $dropped
#'   (every excluded marker with the reason).
verify_and_orient_alleles = function(gwas_dt, bim_path, freq_path, freq_tolerance = 0.05,
                                      ambiguous_freq_margin = 0.05, same_source = FALSE) {
  bim = fread(bim_path, col.names = c("chr", "SNP", "cm", "pos", "A1_bim", "A2_bim"))
  freq = fread(freq_path)
  setnames(freq, c("#CHROM", "ID", "REF", "ALT", "PROVISIONAL_REF?", "ALT_FREQS", "OBS_CT"),
           c("chr_f", "SNP", "REF", "ALT", "prov_ref", "alt_freq", "obs_ct"), skip_absent = TRUE)

  n_input = nrow(gwas_dt)

  # Non-SNP markers (insertions, deletions, or missing-allele codes such as
  # N) are set aside and reported individually rather than stopping the run,
  # since a few malformed rows should not block an otherwise good dataset.
  is_non_snp = nchar(gwas_dt$A1) != 1 | nchar(gwas_dt$A2) != 1 |
    !grepl("^[ACGTacgt]$", gwas_dt$A1) | !grepl("^[ACGTacgt]$", gwas_dt$A2)
  non_snp_dropped = gwas_dt[is_non_snp][, .(SNP, A1, A2, FREQ, orientation = NA_character_,
    drop_reason = "not a SNP (an insertion or deletion, or a missing-allele code such as N); only SNPs are used")]
  gwas_dt = gwas_dt[!is_non_snp]

  # A duplicated marker ID would turn the merges below into a many-to-many
  # join and inflate every count without warning.
  for (nm in c("gwas_dt", "bim", "freq")) {
    dt = get(nm)
    dup = dt$SNP[duplicated(dt$SNP)]
    if (length(dup) > 0) {
      stop("Duplicate marker ID(s) in ", nm, " (e.g. \"", dup[1], "\", ", length(unique(dup)),
           " duplicated ID(s) in total). Marker IDs must be unique. A common cause is a multiallelic ",
           "site split into several biallelic records that were not given distinct IDs.")
    }
  }

  m = merge(gwas_dt, bim[, .(SNP, A1_bim, A2_bim)], by = "SNP")
  m = merge(m, freq[, .(SNP, REF, ALT, alt_freq)], by = "SNP")

  n_id_matched = nrow(m)
  # The merges above silently lose markers whose ID is absent from the
  # genotypes, so they are listed here; otherwise the drop counts would not
  # add up to the difference between input and passing markers.
  not_found_dropped = gwas_dt[!SNP %in% m$SNP][, .(SNP, A1, A2, FREQ, orientation = NA_character_,
    drop_reason = "marker ID not found in the genotype data")]

  # genotype-derived frequency of the .bim A1 allele, which is not
  # necessarily the effect allele
  m[, freq_A1_bim := ifelse(A1_bim == ALT, alt_freq, 1 - alt_freq)]
  # The .bim alleles must be REF and ALT in some order, since both files come
  # from the same genotypes; if not, they were computed from different data.
  bad_ref_alt = m[!((A1_bim == REF & A2_bim == ALT) | (A1_bim == ALT & A2_bim == REF))]
  if (nrow(bad_ref_alt) > 0) {
    stop(nrow(bad_ref_alt), " marker(s) have .bim alleles that do not match REF/ALT from the --freq file ",
         "(e.g. SNP ", bad_ref_alt$SNP[1], "). This means the .bim and the allele-frequency file were not ",
         "computed from the same genotype data. Check that both came from the same uploaded dataset.")
  }

  m[, A1_comp := complement_allele(A1)]
  m[, A2_comp := complement_allele(A2)]
  # For A/T and C/G SNPs, complementing maps the allele pair onto itself, so
  # allele identity cannot tell the strands apart. The two possibilities
  # differ in which allele the frequency belongs to, so the frequency
  # resolves them unless it is too close to 0.5 (the approach of
  # TwoSampleMR::harmonise_data). On the maize data, dropping all 33,008 such
  # markers instead, of which only 2,973 were unresolvable, reduced the COJO
  # signals from 46 to 25.
  m[, is_ambiguous := (A1 == A2_comp) & (A2 == A1_comp)]

  m[, orientation := fcase(
    !is_ambiguous & A1 == A1_bim & A2 == A2_bim, "direct",
    !is_ambiguous & A1 == A2_bim & A2 == A1_bim, "swapped",
    !is_ambiguous & A1_comp == A1_bim & A2_comp == A2_bim, "flipped",
    !is_ambiguous & A1_comp == A2_bim & A2_comp == A1_bim, "flipped_swapped",
    # For A/T and C/G SNPs the frequency decides: "direct" means the
    # declared A1 lines up with the .bim A1, "swapped" with the .bim A2.
    is_ambiguous & (A1 == A1_bim | A1 == A2_comp & A2_comp == A1_bim) &
      abs(FREQ - freq_A1_bim) <= abs(FREQ - (1 - freq_A1_bim)), "ambiguous_direct",
    is_ambiguous & (A1 == A1_bim | A1 == A2_comp & A2_comp == A1_bim), "ambiguous_swapped",
    is_ambiguous & abs(FREQ - (1 - freq_A1_bim)) <= abs(FREQ - freq_A1_bim), "ambiguous_swapped",
    is_ambiguous, "ambiguous_direct",
    default = "no_match"
  )]

  # Near 0.5 both orientations fit the frequency equally well, so the call
  # would be a guess, unless same_source rules out a strand flip and allele
  # identity settles it.
  m[, ambiguous_unresolvable := if (same_source) FALSE else
      is_ambiguous & abs(FREQ - 0.5) < ambiguous_freq_margin]
  if (same_source) {
    m[is_ambiguous == TRUE, orientation := fifelse(A1 == A1_bim, "ambiguous_direct", "ambiguous_swapped")]
  }

  # expected frequency of the declared effect allele, given its orientation
  m[, freq_expected := fcase(
    orientation %in% c("direct", "flipped", "ambiguous_direct"), freq_A1_bim,
    orientation %in% c("swapped", "flipped_swapped", "ambiguous_swapped"), 1 - freq_A1_bim,
    default = NA_real_
  )]

  m[, freq_diff := abs(FREQ - freq_expected)]
  # A supplied frequency that matches 1 - freq_expected but not freq_expected
  # has been reported for the other allele, typically a MAF column, while the
  # alleles themselves are correct. It is corrected rather than dropped,
  # since the correct value is already known from the genotypes. On the
  # maize data all of the roughly 27,000 such markers matched the complement
  # clearly, with no borderline cases. Corrected markers are flagged in
  # freq_corrected.
  m[, freq_looks_swapped := orientation %in% c("direct", "flipped") &
      freq_diff > freq_tolerance & abs(FREQ - (1 - freq_expected)) <= freq_tolerance]
  m[, freq_corrected := freq_looks_swapped]

  m[, drop_reason := fcase(
    orientation == "no_match", "alleles in the GWAS file do not match the genotype data, even allowing for a strand flip",
    ambiguous_unresolvable, sprintf("A/T or C/G marker with frequency %.3f, too close to 0.5 to tell which strand was used", FREQ),
    freq_looks_swapped, NA_character_,  # corrected, not dropped
    !is.na(freq_diff) & freq_diff > freq_tolerance, sprintf("allele frequency in the GWAS file (%.3f) does not agree with the genotype data (%.3f)", FREQ, freq_expected),
    default = NA_character_
  )]

  dropped = rbind(m[!is.na(drop_reason), .(SNP, A1, A2, FREQ, orientation, drop_reason)], non_snp_dropped,
                  not_found_dropped)
  matched = m[is.na(drop_reason)]
  matched[, FREQ := freq_expected]  # the genotype-derived frequency replaces the supplied one

  summary_dt = data.table(
    n_input = n_input,
    n_id_matched = n_id_matched,
    n_dropped_non_snp = nrow(non_snp_dropped),
    n_dropped_not_in_genotypes = nrow(not_found_dropped),
    n_dropped_total = nrow(dropped),
    n_dropped_no_match = sum(m$orientation == "no_match", na.rm = TRUE),
    n_ambiguous_total = sum(m$is_ambiguous, na.rm = TRUE),
    n_ambiguous_resolved_by_freq = sum(m$is_ambiguous & is.na(m$drop_reason), na.rm = TRUE),
    n_dropped_ambiguous_unresolvable = sum(m$ambiguous_unresolvable & !is.na(m$drop_reason), na.rm = TRUE),
    n_dropped_freq_inconsistent = sum(!is.na(m$drop_reason) & !m$ambiguous_unresolvable & m$orientation != "no_match"),
    n_freq_swap_corrected = sum(matched$freq_corrected, na.rm = TRUE),
    n_matched = nrow(matched),
    # Not rounded here: rounding 0.99998 to 1.0000 would report a clean 100%
    # for a run that dropped markers. The app formats it.
    match_rate = nrow(matched) / n_input
  )

  list(matched = matched[, .(SNP, A1, A2, FREQ, freq_corrected, BETA, SE, P, N)],
       summary = summary_dt,
       dropped = dropped)
}
