#!/usr/bin/env Rscript
# One chunk of source IDAT pairs -> ssNoob beta with detection-P masking.
# Supports 450K, EPIC and EPICv2 (minfi auto-detects the array from the IDATs).
# Writes beta for all probes, plus per-sample QC. No clocks, no sample selection.
options(warn=1, digits=17)
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 1L) stop("Expected one JSON configuration path")
suppressPackageStartupMessages({
    library(jsonlite)
    library(minfi)
})
config <- fromJSON(args[[1]], simplifyVector=TRUE)
basenames <- config$basenames
gsms <- config$sample_ids
if (!length(basenames) || length(basenames) != length(gsms) || anyDuplicated(gsms))
    stop("Invalid chunk roster")
threshold <- as.numeric(config$detection_p_threshold)

rg <- read.metharray(basenames, extended=FALSE, force=TRUE, verbose=FALSE)
array <- annotation(rg)[["array"]]
colnames(rg) <- gsms
mset <- preprocessNoob(rg, dyeMethod="single")
beta <- getBeta(mset)
detp <- detectionP(rg)
common <- intersect(rownames(beta), rownames(detp))
if (length(common) < 0.9 * nrow(beta)) stop("Detection-P probes do not cover the beta probes")
beta <- beta[common, gsms, drop=FALSE]
detp <- detp[common, gsms, drop=FALSE]
failed <- detp > threshold | !is.finite(beta)
beta[failed] <- NA

qc <- data.frame(
    sample_id=gsms,
    array=array,
    probes=nrow(beta),
    detection_failed=colSums(detp > threshold),
    nonfinite_beta=colSums(!is.finite(getBeta(mset)[common, gsms, drop=FALSE])),
    masked=colSums(failed),
    median_meth=apply(getMeth(mset)[common, gsms, drop=FALSE], 2, median),
    median_unmeth=apply(getUnmeth(mset)[common, gsms, drop=FALSE], 2, median),
    check.names=FALSE
)
con <- gzfile(config$beta_output, "wt")
write.table(data.frame(ID_REF=rownames(beta), beta, check.names=FALSE), con,
            sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE, na="NA")
close(con)
write.table(qc, config$qc_output, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
write_json(list(status="PASS", array=array, samples=length(gsms), probes=nrow(beta),
                minfi=as.character(packageVersion("minfi")),
                annotation=paste(annotation(rg), collapse="/"),
                method="read.metharray(force=TRUE) -> preprocessNoob(dyeMethod='single') -> getBeta; detectionP > threshold masked to NA",
                detection_p_threshold=threshold),
           config$receipt_output, auto_unbox=TRUE, pretty=TRUE)
