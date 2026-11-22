#!/usr/bin/env Rscript
# Explicit sample-level EPIC reprocessing. This does not reconstruct author averages.
options(warn=1, digits=17)
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=1L) stop("Expected one JSON worker configuration")
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(minfi))
suppressPackageStartupMessages(library(IlluminaHumanMethylationEPICmanifest))
config <- fromJSON(args[[1]], simplifyVector=FALSE)
acq <- fromJSON(config$acquisition, simplifyVector=FALSE)
out <- normalizePath(config$output, mustWork=TRUE)
profile <- config$profile
subset_profile <- identical(profile$name,"minfi_ssnoob_beta100_detection01_supported_v2")
if (!subset_profile && profile$name!="minfi_ssnoob_beta100_detection01_v1") stop("Unknown preprocessing profile")
function_snapshot <- list(readIDAT=illuminaio::readIDAT,readIDAT_nonenc=get("readIDAT_nonenc",asNamespace("illuminaio")),
    read_metharray=read.metharray,preprocessRaw=preprocessRaw,preprocessNoob=preprocessNoob,
    detectionP=detectionP,getBeta=getMethod("getBeta","MethylSet"),getAnnotation=getAnnotation,
    getProbeInfo=get("getProbeInfo",asNamespace("minfi")),getOOB=get("getOOB",asNamespace("minfi")))
ids <- vapply(acq$pairs, function(x) x$sample_id, character(1))
if (anyDuplicated(ids) || !length(ids)) stop("Invalid sample roster")
outputs <- list()
register <- function(role,name) {
    if (role %in% names(outputs)) stop("Duplicate output role")
    outputs[[role]] <<- name
    file.path(out,name)
}
write_json_new <- function(x,role,name) {
    p <- register(role,name)
    if (file.exists(p)) stop("Output already exists")
    write_json(x,p,auto_unbox=TRUE,pretty=TRUE,digits=NA,na="null")
}
write_matrix <- function(x,role,name,index_name="probe_id") {
    p <- register(role,name)
    if (file.exists(p)) stop("Output already exists")
    con <- gzfile(p,"wt")
    on.exit(close(con))
    frame <- data.frame(rownames(x),as.data.frame(x),check.names=FALSE)
    colnames(frame)[1] <- index_name
    write.table(frame,con,sep=",",row.names=FALSE,col.names=TRUE,quote=FALSE,na="NA")
}
write_frame <- function(x,role,name) {
    p <- register(role,name)
    if (file.exists(p)) stop("Output already exists")
    write.table(x,p,sep=",",row.names=FALSE,col.names=TRUE,quote=TRUE,na="NA")
}
decoded <- list(); headers <- list(); all_addresses <- NULL; basenames <- character()
for (pair in acq$pairs) {
    gsm <- pair$sample_id
    decoded[[gsm]] <- list(); headers[[gsm]] <- list()
    for (channel in c("Grn","Red")) {
        path <- pair$channels[[channel]]$path
        z <- illuminaio::readIDAT(path)
        q <- z$Quants
        if (z$versionNumber!=3L || nrow(q)!=z$nSNPsRead || anyDuplicated(rownames(q)) ||
            anyNA(rownames(q)) || !all(c("Mean","SD","NBeads") %in% colnames(q)) ||
            any(!is.finite(q)) || any(q<0)) stop("Invalid raw IDAT fields for ",gsm," ",channel)
        # These public EPIC IDATs store RedGreen=0 in both declared channels.
        # Retain that field without interpreting it; channel ownership is the GEO URL/basename.
        if (length(z$RedGreen)!=1L || !is.finite(z$RedGreen)) stop("Invalid uninterpreted RedGreen header field")
        if (is.null(all_addresses)) all_addresses <- rownames(q)
        if (!setequal(rownames(q),all_addresses)) stop("Raw channel/address sets differ; no silent intersection")
        pos <- z$Unknowns$MostlyA
        if (is.null(pos) || length(pos)!=1L || !nzchar(pos)) stop("Array position unavailable")
        h <- list(sample_id=gsm,channel=channel,version=z$versionNumber,nSNPsRead=z$nSNPsRead,
                  RedGreen=z$RedGreen,Barcode=z$Barcode,ChipType=z$ChipType,position=pos)
        headers[[gsm]][[channel]] <- h
        decoded[[gsm]][[channel]] <- q[all_addresses,,drop=FALSE]
        write_json_new(h,paste0("channel_header:",gsm,":",channel),paste0(gsm,"-",channel,"-header.json"))
        write_matrix(q,paste0("channel_decoded:",gsm,":",channel),paste0(gsm,"-",channel,"-decoded.csv.gz"),"address")
    }
    for (field in c("Barcode","ChipType","position","nSNPsRead")) {
        if (!identical(headers[[gsm]]$Grn[[field]],headers[[gsm]]$Red[[field]])) stop("Paired channel header mismatch")
    }
    basenames <- c(basenames,sub("_Grn\\.idat$","",pair$channels$Grn$path))
}
array_ids <- vapply(headers,function(h) paste(h$Grn$Barcode,h$Grn$position,sep=":"),character(1))
if (anyDuplicated(array_ids)) stop("Distinct GEO samples reference the same physical array position")
manifest <- get("IlluminaHumanMethylationEPICmanifest")
pI <- as.data.frame(getProbeInfo(manifest,type="I"))
pII <- as.data.frame(getProbeInfo(manifest,type="II"))
controls <- as.data.frame(getProbeInfo(manifest,type="Control"))
Iout <- data.frame(probe_id=pI$Name,AddressA=pI$AddressA,AddressB=pI$AddressB,Color=pI$Color)
IIout <- data.frame(probe_id=pII$Name,AddressA=pII$AddressA)
write_frame(Iout,"manifest_I","manifest-I.csv")
write_frame(IIout,"manifest_II","manifest-II.csv")
write_frame(controls,"manifest_controls","manifest-controls.csv")
required_addresses <- unique(c(pI$AddressA,pI$AddressB,pII$AddressA,controls$Address))
missing_addresses <- setdiff(required_addresses,all_addresses)
if (!subset_profile && length(missing_addresses)) stop("Source IDAT lacks ",length(missing_addresses)," required manifest/control addresses")
if (any(!controls$Address %in% all_addresses)) stop("Source IDAT lacks a required control address")
if (anyDuplicated(c(pI$Name,pII$Name))) stop("Manifest probe IDs overlap or duplicate")
if (!all(pI$Color %in% c("Grn","Red"))) stop("Unknown Type I manifest color")
if (!all(c("NEGATIVE","NORM_A","NORM_C","NORM_G","NORM_T") %in% controls$Type)) stop("Required control families absent")
fullI <- pI; fullII <- pII
if (subset_profile) {
    pI <- pI[pI$AddressA %in% all_addresses & pI$AddressB %in% all_addresses,,drop=FALSE]
    pII <- pII[pII$AddressA %in% all_addresses,,drop=FALSE]
    if (!all(c("Grn","Red") %in% pI$Color)) stop("An empirical out-of-band pool is unavailable")
}
rg <- read.metharray(basenames,extended=TRUE,force=FALSE,verbose=TRUE)
if (annotation(rg)[["array"]]!=profile$array) stop("Importer array differs from source EPIC declaration")
sampleNames(rg) <- ids
if (!identical(rownames(rg),all_addresses)) stop("Importer lost or reordered raw addresses")
for (gsm in ids) {
    for (channel in c("Grn","Red")) {
        imported <- if (channel=="Grn") getGreen(rg)[,gsm] else getRed(rg)[,gsm]
        if (!isTRUE(all.equal(as.numeric(imported),as.numeric(decoded[[gsm]][[channel]][,"Mean"]),tolerance=0))) {
            stop("Decoded versus imported channel intensity differs")
        }
    }
}
raw <- preprocessRaw(rg)
rawM <- getMeth(raw); rawU <- getUnmeth(raw)
probes <- rownames(rawM)
if (!identical(probes,c(as.character(pI$Name),as.character(pII$Name)))) stop("Raw methylation rows differ from supported manifest order")
for (gsm in ids) {
    directM <- directU <- setNames(rep(NA_real_,length(probes)),probes)
    for (color in c("Grn","Red")) {
        pp <- pI[pI$Color==color,,drop=FALSE]
        directM[pp$Name] <- decoded[[gsm]][[color]][pp$AddressB,"Mean"]
        directU[pp$Name] <- decoded[[gsm]][[color]][pp$AddressA,"Mean"]
    }
    directM[pII$Name] <- decoded[[gsm]]$Grn[pII$AddressA,"Mean"]
    directU[pII$Name] <- decoded[[gsm]]$Red[pII$AddressA,"Mean"]
    if (!isTRUE(all.equal(unname(directM),as.numeric(rawM[,gsm]),tolerance=0)) ||
        !isTRUE(all.equal(unname(directU),as.numeric(rawU[,gsm]),tolerance=0))) stop("Manifest-to-M/U mapping differs")
}
selected <- probes
annotation_package <- paste0("IlluminaHumanMethylationEPICanno.",annotation(rg)[["annotation"]])
if (subset_profile) {
    ann <- getAnnotation(annotation_package)
    ann_groups <- data.frame(probe_id=as.character(ann$Name),Type=as.character(ann$Type),
                             Color=as.character(ann$Color),stringsAsFactors=FALSE)
    if (anyNA(ann_groups) || anyDuplicated(ann_groups$probe_id) ||
        !identical(ann_groups$probe_id,rownames(ann)) ||
        any(!paste0(ann_groups$Type,ann_groups$Color) %in% c("IGrn","IRed","II"))) {
        stop("Installed annotation correction groups are invalid or duplicated")
    }
    membership <- data.frame(probe_id=c(as.character(fullI$Name),as.character(fullII$Name)),
        Type=c(rep("I",nrow(fullI)),rep("II",nrow(fullII))),
        Color=c(as.character(fullI$Color),rep("",nrow(fullII))),stringsAsFactors=FALSE)
    membership$address_supported <- membership$probe_id %in% probes
    ai <- match(membership$probe_id,ann_groups$probe_id)
    membership$annotation_supported <- !is.na(ai)
    present <- which(membership$annotation_supported)
    if (any(membership$Type[present]!=ann_groups$Type[ai[present]]) ||
        any(membership$Color[present]!=ann_groups$Color[ai[present]])) {
        stop("Installed annotation type/color conflicts with the manifest")
    }
    membership$disposition <- ifelse(!membership$address_supported,"address_absent",
        ifelse(!membership$annotation_supported,"annotation_absent","retained"))
    selected <- membership$probe_id[membership$disposition=="retained"]
    if (!length(selected)) stop("No observed probes have a supported correction group")
    # Keep the full RG object: otherwise the exclusion also changes empirical OOB pools.
    write_frame(ann_groups,"annotation_groups","annotation-groups.csv")
    write_frame(membership,"probe_membership","probe-membership.csv")
    write_frame(data.frame(address=sort(as.numeric(missing_addresses))),"missing_addresses","missing-addresses.csv")
}
det_raw <- detectionP(rg,type="m+u")
if (!setequal(rownames(det_raw),probes) || !identical(colnames(det_raw),ids)) stop("Detection identities differ")
det_raw <- det_raw[probes,,drop=FALSE]
det <- det_raw[selected,,drop=FALSE]
corrected <- preprocessNoob(rg,offset=15,dyeCorr=TRUE,dyeMethod="single",verbose=TRUE)
allM <- getMeth(corrected); allU <- getUnmeth(corrected)
allBeta <- getBeta(corrected,type="Illumina",betaThreshold=0)
if (!identical(dimnames(allM),dimnames(rawM)) || !identical(dimnames(allM),dimnames(allU)) ||
    !identical(dimnames(allM),dimnames(allBeta))) stop("Corrected matrix identities differ")
M <- allM[selected,,drop=FALSE]; U <- allU[selected,,drop=FALSE]; beta <- allBeta[selected,,drop=FALSE]
formula <- pmax(M,0)/(pmax(M,0)+pmax(U,0)+100)
if (!isTRUE(all.equal(beta,formula,tolerance=1e-12))) stop("Explicit beta arithmetic differs")
single_delta <- list()
single_beta <- beta
for (i in seq_along(ids)) {
    one <- preprocessNoob(rg[,i,drop=FALSE],offset=15,dyeCorr=TRUE,dyeMethod="single",verbose=FALSE)
    one_beta <- getBeta(one,type="Illumina",betaThreshold=0)
    if (!identical(rownames(one_beta),probes)) stop("Single-array probe order differs")
    one_beta <- one_beta[selected,,drop=FALSE]
    single_beta[,i] <- one_beta[,1]
    if (!identical(is.finite(one_beta[,1]),is.finite(beta[,i]))) stop("Single-array nonfinite pattern differs")
    finite_pair <- is.finite(beta[,i])
    if (!identical(as.numeric(one_beta[!finite_pair,1]),as.numeric(beta[!finite_pair,i]))) {
        stop("Single-array nonfinite values differ")
    }
    if (!any(finite_pair)) stop("Single-array comparison has no finite values")
    delta <- max(abs(one_beta[finite_pair,1]-beta[finite_pair,i]))
    if (!is.finite(delta) || delta>1e-12) stop("Single-array versus joint ssNoob differs")
    single_delta[[ids[[i]]]] <- delta
}
failed_det <- is.finite(det) & det>0.01
nonfinite_det <- !is.finite(det); nonfinite_beta <- !is.finite(beta)
mask <- failed_det | nonfinite_det | nonfinite_beta
candidate <- beta; candidate[mask] <- NA_real_
if (any(colSums(!is.na(candidate))==0L)) stop("A sample has no finite QC-retained methylation values")
write_matrix(rawM,"raw_meth","raw-meth.csv.gz")
write_matrix(rawU,"raw_unmeth","raw-unmeth.csv.gz")
write_matrix(M,"meth_corrected","meth-corrected.csv.gz")
write_matrix(U,"unmeth_corrected","unmeth-corrected.csv.gz")
write_matrix(det,"detection_p","detection-p.csv.gz")
write_matrix(beta,"beta_unmasked","beta-unmasked.csv.gz")
write_matrix(candidate,"matrix","beta-qc.csv.gz")
if (subset_profile) {
    write_matrix(det_raw,"raw_detection_p","raw-detection-p.csv.gz")
    write_matrix(single_beta,"beta_single_array","beta-single-array.csv.gz")
}
qc <- data.frame(sample_id=ids,probe_count=nrow(beta),detection_failed_count=colSums(failed_det),
    detection_nonfinite_count=colSums(nonfinite_det),beta_nonfinite_count=colSums(nonfinite_beta),
    masked_count=colSums(mask),observed_count=colSums(!mask),
    raw_median_M=apply(rawM,2,median),raw_median_U=apply(rawU,2,median),check.names=FALSE)
write_frame(qc,"sample_qc","sample-qc.csv")
si <- register("session_info","session-info.txt")
capture.output(sessionInfo(),file=si)
function_file <- register("function_source","function-source.R.txt")
capture.output(function_snapshot,file=function_file)
compatibility <- list(channel_headers="PASS",address_sets="PASS",address_count=length(all_addresses),
    importer_values="PASS",manifest_mapping="PASS",required_manifest_addresses=length(required_addresses),
    missing_manifest_addresses=length(missing_addresses),probe_count=length(selected),control_count=nrow(controls),
    typeI_count=nrow(fullI),typeII_count=nrow(fullII),raw_array_ids=as.list(unname(array_ids)))
if (subset_profile) {
    compatibility$raw_probe_count <- length(probes)
    compatibility$annotation_excluded_probe_count <- sum(membership$disposition=="annotation_absent")
    compatibility$address_excluded_probe_count <- sum(membership$disposition=="address_absent")
    compatibility$corrected_probe_count <- length(selected)
}
result <- list(schema="clockbase-dnam-idat-processing-v1",status="PASS",profile=profile,sample_ids=as.list(ids),
    array=annotation(rg)[["array"]],annotation=as.list(annotation(rg)),outputs=outputs,
    compatibility=compatibility,
    sample_vs_joint_max_abs_delta=single_delta,
    packages=as.list(vapply(c("minfi","illuminaio","IlluminaHumanMethylationEPICmanifest",paste0("IlluminaHumanMethylationEPICanno.",annotation(rg)[["annotation"]])),
        function(p) if (requireNamespace(p,quietly=TRUE)) as.character(packageVersion(p)) else "annotation label",character(1))),
    scientific_status="CANNOT VERIFY",author_method_reproduction=FALSE,accepted_results=0,model_calls=0,new_clock_runs=0)
write_json(result,file.path(out,"processing.json"),auto_unbox=TRUE,pretty=TRUE,digits=NA,na="null")
message("IDAT numerical candidate complete: ",length(ids)," samples, ",length(selected)," corrected probes; ",length(probes)," raw probes retained")
