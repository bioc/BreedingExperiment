#' Pedigree of a breeding experiment
#'
#' Gets or sets the pedigree, a table with columns `id`, `sire` and `dam`. The
#' pedigree may include ancestors that were never genotyped; those individuals
#' are what allow [Hmatrix()] to combine pedigree and genomic information.
#'
#' @param x A [BreedingExperiment-class] object.
#' @param value A data frame with columns `id`, `sire` and `dam`.
#' @param ... Unused.
#' @return `pedigree()` returns a [S4Vectors::DataFrame]; the replacement form
#'   returns the updated object.
#' @examples
#' set.seed(1)
#' be <- simulateBreeding(n_ind = 8, n_marker = 20)
#' head(pedigree(be))
#' @name pedigree
#' @aliases pedigree pedigree<-
NULL

#' @rdname pedigree
#' @export
setMethod("pedigree", "BreedingExperiment", function(x, ...) x@pedigree)

#' @rdname pedigree
#' @export
setReplaceMethod("pedigree", "BreedingExperiment", function(x, value) {
    x@pedigree <- .clean_pedigree(value, colnames(x))
    methods::validObject(x)
    x
})

#' Genotype dosages
#'
#' Gets or sets the matrix of allele dosages, a shorthand for the `genotype`
#' assay. Markers are rows and individuals are columns.
#'
#' @param x A [BreedingExperiment-class] object.
#' @param value A numeric matrix of dosages coded 0, 1 or 2.
#' @param ... Unused.
#' @return `genotypes()` returns a numeric matrix; the replacement form returns
#'   the updated object.
#' @examples
#' set.seed(1)
#' be <- simulateBreeding(n_ind = 6, n_marker = 10)
#' genotypes(be)[1:5, 1:4]
#' @name genotypes
#' @aliases genotypes genotypes<-
NULL

#' @rdname genotypes
#' @export
setMethod("genotypes", "BreedingExperiment", function(x, ...) {
    SummarizedExperiment::assay(x, "genotype")
})

#' @rdname genotypes
#' @export
setReplaceMethod("genotypes", "BreedingExperiment", function(x, value) {
    value <- .as_geno_matrix(value)
    SummarizedExperiment::assay(x, "genotype") <- value
    x
})

#' Which individuals are genotyped
#'
#' Returns the identifiers that appear in the genotype matrix. In a single-step
#' analysis the pedigree usually names more individuals than this, and the
#' difference is what [Hmatrix()] reconciles.
#'
#' @param x A [BreedingExperiment-class] object.
#' @param ... Unused.
#' @return A character vector of individual identifiers.
#' @examples
#' set.seed(1)
#' be <- simulateBreeding(n_ind = 8, n_marker = 20)
#' genotyped(be)
#' @name genotyped
#' @aliases genotyped genotyped,BreedingExperiment-method
#' @export
setMethod("genotyped", "BreedingExperiment", function(x, ...) colnames(x))

#' Allele frequencies and marker summaries
#'
#' Computes, for every marker, the frequency of the counted allele, the minor
#' allele frequency, the call rate and the observed heterozygosity.
#'
#' `freq` and `maf` use every dosage that is present, so an expected dosage
#' from imputation contributes to them on the same footing as a call.
#' `callRate` is the proportion of dosages that are not `NA`; after
#' [imputeMarkers()] that is 1 for every marker, so take call rates before
#' imputing if you need them.
#'
#' `heterozygosity` is different, because a heterozygote is a genotype rather
#' than a quantity. A cell counts towards it only if it is a genotype call:
#' its dosage must be a whole number, and it must not be marked in the
#' `imputed` assay that [imputeMarkers()] writes. Everything else is left out
#' of both the numerator and the denominator. Without those two checks an
#' imputed cell would be counted as a non-heterozygote and heterozygosity
#' would appear to fall every time a data set was imputed. A marker with no
#' genotype call at all gets `NA`.
#'
#' Dosages imputed by other software arrive without the `imputed` assay, so
#' only the whole-number check applies to them; pass such dosages in as
#' fractional values, or carry your own mask as an extra assay, if observed
#' heterozygosity matters to you.
#'
#' @param x A [BreedingExperiment-class] object.
#' @param tol How far a dosage may sit from a whole number and still count as
#'   a genotype call when `heterozygosity` is computed.
#' @param ... Unused.
#' @return A [S4Vectors::DataFrame] with one row per marker and columns `freq`,
#'   `maf`, `callRate` and `heterozygosity`.
#' @seealso [filterMarkers()], [imputeMarkers()]
#' @examples
#' set.seed(1)
#' be <- simulateBreeding(n_ind = 30, n_marker = 50)
#' head(alleleFrequency(be))
#'
#' # imputation does not change the observed heterozygosity of a marker,
#' # because imputed cells are not counted as genotype calls
#' bem <- simulateBreeding(n_ind = 30, n_marker = 50, missing = 0.2)
#' all.equal(alleleFrequency(bem)$heterozygosity,
#'           alleleFrequency(imputeMarkers(bem))$heterozygosity)
#' @name alleleFrequency
#' @aliases alleleFrequency alleleFrequency,BreedingExperiment-method
#' @export
setMethod("alleleFrequency", "BreedingExperiment",
          function(x, tol = 1e-8, ...) {
    g <- genotypes(x)
    p <- rowMeans(g, na.rm = TRUE) / 2
    called <- !is.na(g) & abs(g - round(g)) <= tol & !.imputed_mask(x)
    n_called <- rowSums(called)
    het <- rowSums(called & abs(g - 1) <= tol)
    S4Vectors::DataFrame(
        freq = p,
        maf = pmin(p, 1 - p),
        callRate = rowMeans(!is.na(g)),
        heterozygosity = ifelse(n_called > 0L, het / n_called, NA_real_),
        row.names = rownames(g))
})