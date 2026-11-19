# Example dataset

The maize case study included with SweepGWAS as a worked example. The
**Load maize case-study example** button in the app loads it, and the file
paths are set in `example_paths` in `app.R`.

The `.bed`, `.bim` and `.fam` files and the GWAS summary statistics are
stored with Git LFS (see `.gitattributes` at the repository root), about
87 MB in total. A plain `git clone` without Git LFS installed fetches small
pointer files instead of the data, and the example button then fails.

## Files

| File | What it is |
|---|---|
| `gwas_farmcpu_gdd_dts_kansas.csv` | FarmCPU GWAS results for flowering time (`GDD_DTS`, Kansas), 142,779 markers, 2,012 accessions. Unmodified GAPIT output, kept in its original format so that the example also tests the app's column mapping on real output from another program. |
| `genotypes.{bed,bim,fam}` | The genotype panel the GWAS was run on, filtered at a minor-allele frequency of 0.05: 142,779 markers and 2,012 accessions. |
| `selection_scan_regions.csv` | 6,862 selection-scan regions from pcadapt, XP-nSL and RAiSD, with columns `chr,region_id,start,end,n_markers,method`. |

## Expected result

At the 10 kb distance threshold the example should give the following. If your run
differs, something in the setup is wrong, which is the reason for providing a
known answer.

- **46** signals selected by COJO at the genome-wide Bonferroni threshold
  (p < 3.5e-7), of which **24** are significant both jointly and marginally
  and are kept. Their median ratio of joint to marginal effects is 1.00, so
  the joint model is judged stable and LD clumping is not used.
- **9 of 24 (37.5%)** signals lie within 10 kb of a selection-scan region,
  against **36.4%** expected by chance under the null drawn from tested
  markers, **p = 0.524** (1,000 permutations, seed 20260918). This is
  consistent with chance.
- The null-model check reports that the two nulls disagree: a null drawn
  uniformly along each chromosome gives 10.0% by chance and p = 0.000999.
- The run takes about two to three minutes, most of it in COJO.

## What the example loader reconstructs

GAPIT reports no effect-allele column, only the minor-allele frequency, so
`app.R` derives both the effect allele and its frequency rather than assuming
them. Both are easy to get wrong in ways that still produce a plausible
result.

1. **The effect allele** is the allele that `plink2 --export A` counted when
   the dosage matrix used by GAPIT was built. Its column headers show that
   all 142,779 columns counted the `.bim` file's **A2**, which contradicts
   plink2's documented default of counting the ALT allele, because the REF
   allele is provisional in these genotyping-by-sequencing data (the `.afreq`
   file carries a `PROVISIONAL_REF?` column). The file, not the manual, is
   the authority here.
2. **The effect-allele frequency** is computed from the genotypes with
   `plink2 --freq` rather than taken from GAPIT's minor-allele frequency,
   which matches the effect allele's frequency for only about 79% of markers
   in this panel.

Using the `.bim` file's A1 and GAPIT's minor-allele frequency instead changes
the signals COJO selects (45 rather than 46), with no error anywhere.

## Why the same-genotypes option is ticked

The GWAS was run on this exact genotype panel, so a strand flip is impossible.
The example therefore ticks **The GWAS was run on these same genotypes**,
which matches A/T and C/G markers by allele rather than by frequency and drops
none of them. Without it, the 2,973 such markers whose frequency is too close
to 0.5 are dropped, and COJO selects 45 of the 46 signals.

## Provenance

Derived from the ZeaGBSv2.7 maize genotyping-by-sequencing release (Romay et
al. 2013, *Genome Biology* 14(6), R55). Romay et al. (2013) report
linkage-disequilibrium decay of about 1 kb across the whole collection and
about 10 kb within breeding groups, and the example uses 10 kb.
