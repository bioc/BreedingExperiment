#' Create a breeding experiment
#'
#' Builds a [BreedingExperiment-class] from a marker genotype matrix, with
#' optional phenotypes, pedigree and marker coordinates.
#'
#' Genotypes are allele dosages: usually the integers 0, 1 and 2, with `NA` for
#' missing calls. Fractional values in between are accepted, so expected dosages
#' from imputation can be stored directly. Markers are rows and individuals are
#' columns, matching the Bioconductor convention that features are rows.
#'
#' If `rowRanges` is not supplied, markers are given placeholder coordinates on
#' a single sequence so that the object is still a valid ranged experiment; pass
#' real coordinates whenever you have them, because that is what allows a data
#' set to be subset by genomic region.
#'
#' The pedigree may name ancestors that were never genotyped. Those extra
#' individuals are exactly what [Hmatrix()] needs in order to combine pedigree
#' and genomic information.
#'
#' @section How phenotypes are matched to individuals:
#' If `phenotypes` has row names, they are treated as individual identifiers:
#' they must be the same set as the column names of `genotypes`, and the rows
#' are reordered to match. If it has no row names, or only the automatic
#' `1, 2, 3, ...` of a freshly built data frame, the rows are taken in the
#' order given. Row names that are neither absent nor matching are an error
#' rather than something to be overwritten, because the one thing worse than
#' refusing a phenotype table is silently attaching each phenotype to the
#' wrong individual.
#'
#' @param genotypes A numeric matrix of allele dosages, markers in rows and
#'   individuals in columns, with column names identifying the individuals.
#' @param phenotypes Optional data frame or [S4Vectors::DataFrame] of
#'   individual-level variables, with one row per column of `genotypes`.
#' @param pedigree Optional data frame with columns `id`, `sire` and `dam`.
#'   Unknown parents may be `NA` or `0`.
#' @param rowRanges Optional [GenomicRanges::GRanges-class] giving the position
#'   of each marker, of the same length as `nrow(genotypes)`.
#' @param ... Further assays, passed to
#'   [SummarizedExperiment::SummarizedExperiment()].
#' @return A [BreedingExperiment-class] object.
#' @seealso [Gmatrix()], [Amatrix()], [Hmatrix()], [simulateBreeding()]
#' @examples
#' set.seed(1)
#' geno <- matrix(rbinom(40 * 12, 2, 0.3), nrow = 40, ncol = 12,
#'                dimnames = list(paste0("snp", 1:40), paste0("ind", 1:12)))
#' pheno <- data.frame(yield = rnorm(12))
#' be <- BreedingExperiment(geno, phenotypes = pheno)
#' be
#' dim(be)
#'
#' # phenotypes given in a different order are matched by row name
#' shuffled <- data.frame(yield = rnorm(12),
#'                        row.names = rev(colnames(geno)))
#' identical(colnames(BreedingExperiment(geno, phenotypes = shuffled)),
#'           colnames(geno))
#' @export
BreedingExperiment <- function(genotypes, phenotypes = NULL, pedigree = NULL,
                               rowRanges = NULL, ...) {
    genotypes <- .as_geno_matrix(genotypes)
    n_marker <- nrow(genotypes)
    n_ind <- ncol(genotypes)

    if (is.null(colnames(genotypes))) {
        colnames(genotypes) <- paste0("ind", seq_len(n_ind))
    }
    if (is.null(rownames(genotypes))) {
        rownames(genotypes) <- paste0("marker", seq_len(n_marker))
    }
    if (anyDuplicated(colnames(genotypes)) > 0L) {
        stop("`genotypes` has duplicated column names, so individuals ",
             "cannot be told apart.", call. = FALSE)
    }

    cd <- if (is.null(phenotypes)) {
        S4Vectors::DataFrame(row.names = colnames(genotypes))
    } else {
        pd <- methods::as(as.data.frame(phenotypes), "DataFrame")
        if (nrow(pd) != n_ind) {
            stop("`phenotypes` must have one row per individual (",
                 n_ind, ").", call. = FALSE)
        }
        .align_phenotypes(phenotypes, pd, colnames(genotypes))
    }

    rr <- if (is.null(rowRanges)) {
        GenomicRanges::GRanges(
            seqnames = rep("unknown", n_marker),
            ranges = IRanges::IRanges(start = seq_len(n_marker), width = 1L))
    } else {
        if (length(rowRanges) != n_marker) {
            stop("`rowRanges` must have one range per marker (",
                 n_marker, ").", call. = FALSE)
        }
        rowRanges
    }
    names(rr) <- rownames(genotypes)

    se <- SummarizedExperiment::SummarizedExperiment(
        assays = c(list(genotype = genotypes), list(...)),
        rowRanges = rr, colData = cd)

    ped <- .clean_pedigree(pedigree, colnames(genotypes))
    methods::new("BreedingExperiment", se, pedigree = ped)
}

## Match phenotype rows to individuals by name where there are names to match
## on, and by position only where there are not. See the constructor's
## documentation for why a mismatch is an error.
#' @keywords internal
#' @noRd
.align_phenotypes <- function(raw, pd, ids) {
    rn <- rownames(raw)
    ## a data frame built without row.names reports "1", "2", ... which say
    ## nothing about which individual a row describes; .row_names_info()
    ## returns a negative count for exactly that case
    automatic <- is.data.frame(raw) && .row_names_info(raw) < 0L
    if (automatic || is.null(rn)) {
        rownames(pd) <- ids
        return(pd)
    }
    if (setequal(rn, ids)) {
        pd <- pd[match(ids, rn), , drop = FALSE]
        rownames(pd) <- ids
        return(pd)
    }
    stop("the row names of `phenotypes` do not identify the individuals in ",
         "`genotypes`. Either give `phenotypes` row names matching the ",
         "genotype column names, or remove them to match by position. ",
         "Not matched: ",
         .some(c(setdiff(ids, rn), setdiff(rn, ids))), ".", call. = FALSE)
}

#' @keywords internal
#' @noRd
.as_geno_matrix <- function(g) {
    if (is.data.frame(g)) g <- as.matrix(g)
    if (!is.matrix(g)) stop("`genotypes` must be a matrix.", call. = FALSE)
    storage.mode(g) <- "double"
    ## Dosages are usually the integers 0, 1 and 2, but imputation - whether by
    ## imputeMarkers() or by external software reporting expected dosages -
    ## produces values in between, so anything within [0, 2] is accepted.
    bad <- !is.na(g) & (g < 0 | g > 2)
    if (any(bad)) {
        stop("`genotypes` must be allele dosages between 0 and 2 (or NA).",
             call. = FALSE)
    }
    g
}

#' @keywords internal
#' @noRd
.clean_pedigree <- function(pedigree, ids) {
    if (is.null(pedigree)) {
        return(S4Vectors::DataFrame(id = character(), sire = character(),
                                    dam = character()))
    }
    p <- as.data.frame(pedigree, stringsAsFactors = FALSE)
    if (!all(.PED_COLS %in% names(p))) {
        stop("`pedigree` must have columns 'id', 'sire' and 'dam'.",
             call. = FALSE)
    }
    p <- p[, .PED_COLS, drop = FALSE]
    for (nm in .PED_COLS) p[[nm]] <- as.character(p[[nm]])
    ## treat "0" and "" as unknown, the convention in pedigree files
    p$sire[p$sire %in% c("0", "")] <- NA_character_
    p$dam[p$dam %in% c("0", "")] <- NA_character_
    if (anyNA(p$id)) {
        stop("`pedigree` has rows with a missing id.", call. = FALSE)
    }
    p <- .resolve_duplicate_ids(p)

    ## parents that never appear as an id are added as founders, so the
    ## pedigree is self-contained and can be sorted
    parents <- unique(c(p$sire, p$dam))
    parents <- parents[!is.na(parents) & !(parents %in% p$id)]
    if (length(parents)) {
        p <- rbind(data.frame(id = parents, sire = NA_character_,
                              dam = NA_character_, stringsAsFactors = FALSE), p)
    }
    S4Vectors::DataFrame(p)
}

## Repeating a pedigree row changes nothing and is worth no more than a
## warning. Two rows that give the same individual different parents are a
## data error, and keeping whichever happened to come first would hide it.
#' @keywords internal
#' @noRd
.resolve_duplicate_ids <- function(p) {
    dup <- unique(p$id[duplicated(p$id)])
    if (!length(dup)) return(p)
    conflicting <- vapply(dup, function(i) {
        nrow(unique(p[p$id == i, , drop = FALSE])) > 1L
    }, logical(1))
    if (any(conflicting)) {
        stop("`pedigree` gives different parents for the same id: ",
             .some(dup[conflicting]), ".", call. = FALSE)
    }
    warning(length(dup), " pedigree id(s) appeared more than once with the ",
            "same parents; the repeats were removed: ", .some(dup),
            call. = FALSE)
    p[!duplicated(p$id), , drop = FALSE]
}