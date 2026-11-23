#' Filter markers on quality
#'
#' Removes markers that are uninformative or poorly genotyped: those below a
#' minor allele frequency threshold, those called in too few individuals, and
#' optionally those that are monomorphic.
#'
#' @section Choosing the thresholds:
#' `min_maf` and `min_call_rate` are judgement calls, not constants of nature.
#' Common starting points are a minor allele frequency of 0.01 to 0.05 and a
#' call rate of 0.90 to 0.95. Raise `min_maf` when the panel is large and the
#' sample small, because a rare allele seen in one or two individuals carries
#' almost no information about relationships but still adds noise; lower it
#' when rare variants are the point of the study. Run this before
#' [imputeMarkers()], so that imputation is not asked to rescue markers that
#' should simply be discarded.
#'
#' @param x A [BreedingExperiment-class] object.
#' @param min_maf Minimum minor allele frequency to keep.
#' @param min_call_rate Minimum proportion of individuals with a call.
#' @param drop_monomorphic Drop markers with no variation.
#' @param ... Unused.
#' @return A [BreedingExperiment-class] object with fewer rows.
#' @seealso [imputeMarkers()], [alleleFrequency()]
#' @examples
#' set.seed(1)
#' be <- simulateBreeding(n_ind = 30, n_marker = 100, missing = 0.05)
#' nrow(be)
#' nrow(filterMarkers(be, min_maf = 0.05, min_call_rate = 0.9))
#' @name filterMarkers
#' @aliases filterMarkers filterMarkers,BreedingExperiment-method
#' @export
setMethod("filterMarkers", "BreedingExperiment",
          function(x, min_maf = 0.01, min_call_rate = 0.9,
                   drop_monomorphic = TRUE, ...) {
    st <- alleleFrequency(x)
    keep <- st$maf >= min_maf & st$callRate >= min_call_rate
    if (drop_monomorphic) {
        g <- genotypes(x)
        varies <- apply(g, 1L, function(v) {
            v <- v[!is.na(v)]
            length(v) > 0L && stats::sd(v) > 0
        })
        keep <- keep & varies
    }
    keep[is.na(keep)] <- FALSE
    if (!any(keep)) stop("no markers passed the filter.", call. = FALSE)
    x[keep, ]
})

#' Impute missing genotypes
#'
#' Fills missing dosages with twice the allele frequency of the marker, the
#' simple mean imputation that relationship matrices assume. It is adequate for
#' the low rates of missingness left after quality control; for anything more,
#' impute with dedicated software before building the object.
#'
#' Mean imputation leaves an expected dosage rather than a whole number of
#' copies, which the container accepts. Note what that costs: an imputed cell
#' is no longer a genotype call, and its value alone no longer says so, since
#' a marker whose mean dosage is exactly 1 imputes to exactly 1. The cells
#' that were filled in are therefore recorded in a logical assay named
#' `imputed`, reachable with `assay(x, "imputed")` and carried through
#' subsetting with everything else. [alleleFrequency()] consults it, so
#' observed heterozygosity is unchanged by imputing. `callRate`, on the other
#' hand, counts dosages that are present and so becomes 1 for every marker;
#' record call rates before imputing if you need them.
#'
#' A marker with no calls at all has no frequency to impute from. Such markers
#' are dropped with a warning rather than being filled in, because filling
#' them would invent a genotype for every individual at that marker.
#' [filterMarkers()] removes them too, and running it first is the better
#' habit.
#'
#' @param x A [BreedingExperiment-class] object.
#' @param drop_uncalled Drop markers with no genotype calls at all, with a
#'   warning. If `FALSE`, their presence is an error.
#' @param ... Unused.
#' @return A [BreedingExperiment-class] object with no missing genotypes and
#'   an `imputed` assay marking the cells that were filled in.
#' @seealso [filterMarkers()], [alleleFrequency()]
#' @examples
#' set.seed(1)
#' be <- simulateBreeding(n_ind = 20, n_marker = 50, missing = 0.1)
#' anyNA(genotypes(be))
#' anyNA(genotypes(imputeMarkers(be)))
#' @name imputeMarkers
#' @aliases imputeMarkers imputeMarkers,BreedingExperiment-method
#' @export
setMethod("imputeMarkers", "BreedingExperiment",
          function(x, drop_uncalled = TRUE, ...) {
    g <- genotypes(x)
    uncalled <- rowSums(!is.na(g)) == 0L
    if (any(uncalled)) {
        ids <- .some(rownames(g)[uncalled])
        if (!isTRUE(drop_uncalled)) {
            stop(sum(uncalled), " marker(s) have no genotype calls at all, ",
                 "so there is nothing to impute from: ", ids,
                 ". Remove them, or use drop_uncalled = TRUE.",
                 call. = FALSE)
        }
        warning(sum(uncalled), " marker(s) with no genotype calls were ",
                "dropped rather than imputed: ", ids, call. = FALSE)
        x <- x[!uncalled, ]
        g <- genotypes(x)
    }
    if (!anyNA(g)) return(x)
    mask <- is.na(g)
    mu <- rowMeans(g, na.rm = TRUE)
    idx <- which(mask, arr.ind = TRUE)
    g[idx] <- mu[idx[, "row"]]
    filled <- .merge_imputed(x, mask)
    genotypes(x) <- g
    SummarizedExperiment::assay(x, .IMPUTED) <- filled
    x
})