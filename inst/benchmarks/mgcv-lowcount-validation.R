#!/usr/bin/env Rscript

# Reproducible bam/fREML entry for the two prespecified low-count validations.
case <- Sys.getenv("MGCVST_BAM_VALIDATION_CASE","power")
if(case=="power") {
  Sys.setenv(MGCVST_BAM_REPS="200",MGCVST_BAM_MEAN="0.3",
    MGCVST_BAM_RHO="0.7",MGCVST_BAM_SEED="361000",
    MGCVST_BAM_CACHE="artifacts/lowcount-investigation/power-200/cache",
    MGCVST_BAM_OUTPUT="artifacts/lowcount-investigation/bam-power-200")
} else if(case=="stress") {
  Sys.setenv(MGCVST_BAM_REPS="200",MGCVST_BAM_MEAN="0.1",
    MGCVST_BAM_RHO="0",MGCVST_BAM_SEED="261000",
    MGCVST_BAM_CACHE="artifacts/lowcount-investigation/stress-200/cache",
    MGCVST_BAM_OUTPUT="artifacts/lowcount-investigation/bam-stress-200")
} else stop("MGCVST_BAM_VALIDATION_CASE must be power or stress")
Sys.setenv(MGCVST_BAM_WORKERS="1",MGCVST_BAM_ORACLE_REPS="2")
sys.source("inst/benchmarks/mgcv-low-count-null.R",envir=new.env(parent=globalenv()))
