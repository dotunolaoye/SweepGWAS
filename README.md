# SweepGWAS

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24380048.svg)](https://doi.org/10.5281/zenodo.24380048)

An R/Shiny application that tests whether independent GWAS signals lie closer
to selection-scan regions than chance would predict. It checks the GWAS
alleles against the genotypes, resolves independent signals, and runs a
permutation test whose chance baseline is drawn from the markers that were
actually tested, with a distance threshold that the user supplies from the
reported linkage-disequilibrium decay of their population.

SweepGWAS runs in the browser, with nothing to install, at
https://dotunolaoye.shinyapps.io/SweepGWAS/. The hosted app suits trying the
tool and moderate datasets. For large datasets, run it on your own computer,
as described below.

## What it does

1. **Allele check.** Checks each marker's effect allele against the genotype
   data, compares the supplied allele frequency with one computed from the
   genotypes, corrects what it can, and reports every dropped marker with the
   reason.
2. **Independent signals.** Runs GCTA-COJO and keeps only signals that are
   significant both jointly and marginally ([Rowan et al. 2024](https://doi.org/10.1111/eva.13666)). If the joint
   model is unstable, meaning that the median ratio of joint to marginal
   effects among the kept signals exceeds 2, independent loci are taken by
   LD clumping in PLINK2 instead.
3. **Convergence test.** Counts the signals that lie within the distance threshold of a
   selection-scan region and compares the count with 1,000 permutations in
   which each signal is redrawn from the tested markers on its own chromosome.
   A null that draws positions uniformly along each chromosome is run beside
   it as a check, and the app states when the two would lead to different
   conclusions.
4. **Region placement diagnostic.** A logistic model of whether each tested
   marker lies within the distance threshold of a selection-scan region, from marker
   density, the gap to the nearest marker, allele frequency, effect size and
   association strength, with the local recombination rate when you upload
   a genetic map (chromosome, bp and cM, as CSV or PLINK .map) and gene
   density when you upload a gene annotation (GFF3 or GTF, or a CSV of gene
   starts). It shows where the regions lie relative to the markers,
   and so in which direction a uniform null would err, with its accuracy
   measured on whole chromosomes withheld in turn.

Each result is explained in a sentence generated from the run, and the app
shows a marginal against joint significance plot for the COJO signals, a
Manhattan plot of every tested marker with the independent signals marked and
the selection-scan regions beneath it, one track per method, for all
chromosomes or any single chromosome, and the permutation distribution. Every
table can be downloaded, the whole analysis can be downloaded as an HTML or PDF
report with a Manhattan plot of each chromosome, and all results, both reports
and a record of the run parameters and stage times come as one zip archive.

## Running it on your own computer

1. **Download SweepGWAS.** On this page, click **Code**, then **Download ZIP**,
   and unzip it. With Git, `git clone` also works, but Git LFS must be
   installed, or the example data arrive as small placeholder files.
2. **Install R** from [CRAN](https://cran.r-project.org), then the packages
   the app uses:

   ```r
   install.packages(c("shiny", "bslib", "DT", "shinyjs", "plotly",
                      "data.table", "ggplot2", "patchwork"))
   ```

3. **Download PLINK 2 and GCTA**, the two programs the app calls, from
   [PLINK 2](https://www.cog-genomics.org/plink/2.0/) and
   [GCTA](https://yanglab.westlake.edu.cn/software/gcta/), and unzip each.
4. **Start the app** in R, with the paths changed to where you put the files:

   ```r
   setwd("path/to/SweepGWAS")
   Sys.setenv(PLINK2_PATH = "path/to/plink2",
              GCTA_PATH = "path/to/gcta64")
   shiny::runApp()
   ```

   If both programs are already on your system path, the `Sys.setenv` line
   is not needed.
5. **Check the setup.** Click **Load maize case-study example** and compare
   the result with the one given in
   [`example_data/README.md`](example_data/README.md). The run takes about two
   to three minutes, most of it in COJO.

## Input requirements

| Input | Notes |
|---|---|
| GWAS summary statistics | Any column names, mapped in the app. Needs marker ID, chromosome, position, effect allele, other allele, effect-allele frequency, effect size, p-value and sample size. The standard error is optional and is calculated from the effect size and p-value when absent. |
| Genotype data | PLINK `bed/bim/fam` or VCF, used as the LD reference and to check alleles. SNPs only; insertions, deletions and multiallelic codes are dropped. |
| Selection-scan regions | CSV with `chr`, `start`, `end` and an optional `method` label. A signal counts as convergent if it lies within the distance threshold of any one supplied region, from any method. |
| Distance threshold | **Required, with no default.** The linkage-disequilibrium decay distance reported for your population, or one estimated with [PopLDdecay](https://github.com/BGI-shenzhen/PopLDdecay). |

Two inputs are easy to get wrong without any visible error:

- **The effect allele** must be the allele the effect sizes refer to, which is
  not necessarily the minor allele or the A1 column of your genotype file.
- **The frequency column** must give that same allele's frequency, not the
  minor-allele frequency. The two differ whenever the effect allele is the
  major allele.

The allele check catches and corrects many such mismatches, but it can only
work from what you supply.

## Design decisions

**The chance baseline is drawn from tested markers.** A GWAS signal can only
arise at a marker that was tested, and a selection-scan region can only be
called where markers exist. A null that places signals uniformly along each
chromosome ignores this and can misstate chance in either direction, for
example where both signals and regions sit where genotyping-by-sequencing
markers are dense.

**The distance threshold has no default.** An arbitrary window can place most of the
genome within range of a region, so that chance convergence is close to
complete and the test can detect nothing. The app stops and asks for a distance threshold rather than guessing one.

**Strand-ambiguous markers are resolved, not discarded.** A/T and C/G SNPs
cannot be matched by allele alone, but their frequency resolves the strand
unless it is too close to 0.5. Dropping all such markers, the usual practice,
can remove a large share of the data and change the signals COJO selects. If
the GWAS was run on the same genotypes you upload, tick
**The GWAS was run on these same genotypes**: a strand flip is then
impossible and none are dropped.

## Not in the app

`R/phenotype_permutation.R` reruns the GWAS on shuffled phenotypes to test the
association signal itself. Each replicate is a full GWAS, so it is run directly from an R session on a local
machine rather than through the app.

## Layout

```
SweepGWAS/
├── app.R                        interface and server
├── R/
│   ├── allele_check.R           allele orientation and frequency check
│   ├── annotation_inputs.R      genetic map and gene annotation inputs
│   ├── cojo_wrapper.R           GCTA-COJO, the signal rules and LD clumping
│   ├── convergence_test.R       permutation test under both nulls
│   ├── convert_vcf_to_bed.R     VCF to PLINK conversion
│   ├── manhattan.R              the Manhattan plots
│   ├── phenotype_permutation.R  phenotype permutation, run locally only
│   ├── plots.R                  the interactive plots
│   ├── region_placement.R       the region placement diagnostic
│   └── report.R                 the HTML and PDF reports
├── www/                         the Guide tab, its figures and the PDF guide
├── example_data/                maize example dataset (Git LFS)
└── setup_binaries.sh            places PLINK2 and GCTA in bin/ for a server
```

The functions in `R/` are general and can be used on any dataset.

## How to cite

Olaoye, D. (2026). *SweepGWAS*. Zenodo.
https://doi.org/10.5281/zenodo.24380048

## Licence

SweepGWAS is released under the MIT License (see [`LICENSE`](LICENSE)).
PLINK2 and GCTA, which the app calls, are separate programs under their own
licences and are not included in this repository.

## Contact

Contact: Dotun Olaoye, University of Missouri-Columbia. Questions and bug reports: dsolaoye@missouri.edu.
