#' Simulate a small breeding population
#'
#' Generates a pedigreed population with markers and a heritable trait, for
#' examples, tests and teaching. Founders are unrelated; later generations are
#' produced by mating sampled parents, and each offspring inherits one allele
#' from each parent at every marker, so the simulated genotypes and the pedigree
#' agree with each other.
#'
#' Set `genotyped` below `n_ind` to leave some individuals ungenotyped, which is
#' the situation [Hmatrix()] is designed for. The genotyped individuals are the
#' most recent ones, so the ungenotyped are the earlier generations, as in a
#' programme that began genotyping part way through.
#'
#' Note that `n_founder` defaults to 10 and is capped at `n_ind`, so a call
#' such as `simulateBreeding(n_ind = 8)` returns eight unrelated founders and
#' no pedigree structure at all. Ask for more individuals than founders to get
#' relatives.
#'
#' The trait is scaled so that the genetic variance is `h2` and the residual
#' variance is `1 - h2`, and `trueBV` in the column data is the genetic
#' component on that scale.
#'
#' The function draws random numbers but does not set the seed itself. Call
#' [set.seed()] beforehand if you need the same population twice.
#'
#' @param n_ind Number of individuals in total.
#' @param n_marker Number of markers.
#' @param n_founder Number of unrelated founders. Capped at `n_ind`; at least
#'   two are needed if there are to be any offspring.
#' @param genotyped Number of individuals with genotypes, taken from the most
#'   recent ones. Defaults to all of them.
#' @param n_qtl Number of markers with an effect on the trait.
#' @param h2 Narrow-sense heritability of the simulated trait, between 0 and 1.
#' @param missing Proportion of genotype calls set to `NA`, between 0 and 1.
#' @return A [BreedingExperiment-class] object with a `phenotype` column and the
#'   true breeding value in `colData`.
#' @seealso [BreedingExperiment()], [Hmatrix()]
#' @examples
#' set.seed(1)
#' be <- simulateBreeding(n_ind = 25, n_marker = 100)
#' be
#' head(as.data.frame(colData(be)))
#'
#' # a population where only the recent individuals were genotyped
#' part <- simulateBreeding(n_ind = 30, n_marker = 100, n_founder = 6,
#'                          genotyped = 20)
#' c(pedigree = nrow(pedigree(part)), genotyped = ncol(part))
#' @export
simulateBreeding <- function(n_ind = 50, n_marker = 200, n_founder = 10,
                             genotyped = NULL, n_qtl = 20, h2 = 0.4,
                             missing = 0) {
    n_ind <- .one_count(n_ind, "n_ind")
    n_marker <- .one_count(n_marker, "n_marker")
    n_founder <- .one_count(n_founder, "n_founder")
    n_qtl <- .one_count(n_qtl, "n_qtl")
    h2 <- .one_proportion(h2, "h2")
    missing <- .one_proportion(missing, "missing")

    n_founder <- min(n_founder, n_ind)
    if (n_ind > n_founder && n_founder < 2L) {
        stop("at least two founders are needed to produce offspring, but ",
             "`n_founder` is ", n_founder, " and `n_ind` is ", n_ind, ".",
             call. = FALSE)
    }
    n_geno <- if (is.null(genotyped)) n_ind else
        .one_count(genotyped, "genotyped")
    if (n_geno > n_ind) {
        stop("`genotyped` cannot exceed `n_ind` (", n_ind, ").",
             call. = FALSE)
    }

    pop <- .simulate_pedigree_genotypes(n_ind, n_marker, n_founder)
    trait <- .simulate_trait(pop$geno, n_qtl, h2)

    keep <- seq.int(n_ind - n_geno + 1L, n_ind)
    g <- pop$geno[, keep, drop = FALSE]
    if (missing > 0) {
        n_na <- round(length(g) * missing)
        if (n_na > 0) g[sample.int(length(g), n_na)] <- NA_real_
    }

    BreedingExperiment(
        genotypes = g,
        phenotypes = data.frame(phenotype = trait$pheno[keep],
                                trueBV = trait$bv[keep]),
        pedigree = pop$pedigree,
        rowRanges = .simple_map(n_marker))
}

#' @keywords internal
#' @noRd
.one_count <- function(value, name) {
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
            value < 1 || value != round(value)) {
        stop("`", name, "` must be a single positive integer.", call. = FALSE)
    }
    as.integer(value)
}

#' @keywords internal
#' @noRd
.one_proportion <- function(value, name) {
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
            value < 0 || value > 1) {
        stop("`", name, "` must be a single number between 0 and 1.",
             call. = FALSE)
    }
    as.numeric(value)
}

## sample() reads a length-one vector as 1:n, so sampling from a pool that has
## narrowed to one candidate would draw an individual that is not in the pool.
#' @keywords internal
#' @noRd
.pick <- function(v) v[sample.int(length(v), 1L)]

## Drop genotypes through a randomly mated pedigree, so markers and pedigree
## tell the same story.
#' @keywords internal
#' @noRd
.simulate_pedigree_genotypes <- function(n_ind, n_marker, n_founder) {
    ids <- sprintf("ind%03d", seq_len(n_ind))
    p <- stats::runif(n_marker, 0.1, 0.9)
    geno <- matrix(0, n_marker, n_ind,
                   dimnames = list(sprintf("snp%04d", seq_len(n_marker)), ids))

    ## founders are unrelated, in Hardy-Weinberg proportions
    for (j in seq_len(n_founder)) {
        geno[, j] <- stats::rbinom(n_marker, 2, p)
    }
    sire <- rep(NA_character_, n_ind)
    dam <- rep(NA_character_, n_ind)

    ## everyone later inherits one allele from each parent
    if (n_ind > n_founder) {
        for (j in seq.int(n_founder + 1L, n_ind)) {
            pool <- seq_len(j - 1L)
            s <- .pick(pool)
            d <- .pick(setdiff(pool, s))
            sire[j] <- ids[s]
            dam[j] <- ids[d]
            geno[, j] <- .gamete(geno[, s]) + .gamete(geno[, d])
        }
    }
    list(geno = geno,
         pedigree = data.frame(id = ids, sire = sire, dam = dam,
                               stringsAsFactors = FALSE))
}

## A trait controlled by a handful of markers, scaled so that the genetic and
## residual variances are h2 and 1 - h2. Scaling both parts, rather than
## inflating the noise, keeps h2 = 0 and h2 = 1 meaningful.
#' @keywords internal
#' @noRd
.simulate_trait <- function(geno, n_qtl, h2) {
    n_marker <- nrow(geno)
    qtl <- sort(sample.int(n_marker, min(n_qtl, n_marker)))
    eff <- stats::rnorm(length(qtl))
    bv <- as.numeric(crossprod(geno[qtl, , drop = FALSE], eff))
    s <- stats::sd(bv)
    bv <- if (is.na(s) || s == 0) rep(0, length(bv)) else
        sqrt(h2) * (bv - mean(bv)) / s
    list(bv = bv,
         pheno = bv + stats::rnorm(ncol(geno), sd = sqrt(1 - h2)))
}

## Placeholder marker coordinates spread over five sequences.
#' @keywords internal
#' @noRd
.simple_map <- function(n_marker) {
    GenomicRanges::GRanges(
        seqnames = rep(paste0("chr", 1 + (seq_len(n_marker) - 1) %% 5),
                       length.out = n_marker),
        ranges = IRanges::IRanges(start = seq_len(n_marker) * 1000L,
                                  width = 1L))
}

#' @keywords internal
#' @noRd
.gamete <- function(dosage) {
    ## a parent of dosage 0 or 2 always transmits the same allele; a
    ## heterozygote transmits either with probability one half
    ifelse(dosage == 1, stats::rbinom(length(dosage), 1L, 0.5), dosage / 2)
}