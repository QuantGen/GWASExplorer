# Cannot assign to global because of staged installation, see
# https://developer.r-project.org/Blog/public/2019/02/14/staged-install/
getDataPath <- function() {
    system.file("extdata", package = "GWASExplorer")
}

loadChromosomeLengths <- function() {
    chromosomeLengths <- arrow::read_parquet(paste0(getDataPath(), "/chromosome_lengths.parquet"))
    chromosomeLengths <- chromosomeLengths[!chromosomeLengths$chromosome %in% c("X", "Y"), ]
    chromosomeLengths$chromosome <- factor(chromosomeLengths$chromosome, levels = 1:22)
    return(chromosomeLengths)
}

loadSegments <- function(pValueRange = c(-log10(5e-8), 60)) {
    segments <- arrow::read_parquet(paste0(getDataPath(), "/segments.parquet"))
    segments$chromosome <- factor(
        segments$chromosome,
        levels = 1:26,
        labels = c(1:22, "X", "Y", "XY", "MT")
    )
    segments <- segments[!segments$chromosome %in% c("X", "Y", "XY", "MT"), ]
    segments$pValue <- -log10(segments$pValue)
    # Remove p-values outside of the lower range
    segments <- segments[segments$pValue >= pValueRange[1], ]
    # Truncate upper p-value range
    segments$pValue[is.infinite(segments$pValue) | segments$pValue > pValueRange[2]] <- pValueRange[2]
    return(segments)
}

extractSegments <- function(segments, trait) {
    segments[segments$trait == trait, ]
}

loadVariants <- function(pValueRange = c(0, 323)) {
    variants <- arrow::read_parquet(paste0(getDataPath(), "/variants.parquet"))
    variants$chromosome <- factor(variants$chromosome, levels = 1:26, labels = c(1:22, "X", "Y", "XY", "MT"))
    variants <- variants[!variants$chromosome %in% c("X", "Y", "XY", "MT"), ]
    variants$LBR_BFDR <- BGLR::BFDR(variants$LBR_d)
    variants$SMR_y <- -log10(variants$SMR_y)
    # Truncate upper limit of p-values
    variants$SMR_y[is.infinite(variants$SMR_y)] <- pValueRange[2]
    return(variants)
}

extractVariants <- function(variants, trait, cutoff = 1) {
    variants[variants$trait == trait & variants$LBR_BFDR < cutoff, ]
}

range01 <- function(x) {
    (x - min(x)) / (max(x) - min(x))
}

shiftceil <- function(lengths) {
    maxLength <- max(lengths)
    nDigits <- floor(log10(maxLength))
    round(maxLength, digits = -(nDigits - 1))
}
