suppressPackageStartupMessages(library(mgcvST))
repo <- Sys.getenv("MGCVST_BASELINE")
if (!dir.exists(file.path(repo,"R"))) stop("Set MGCVST_BASELINE to the unmodified b5195dd source.")
out <- Sys.getenv("MGCVST_BENCH_OUT",tempdir())
dir.create(out,showWarnings=FALSE,recursive=TRUE)
old <- new.env(parent=asNamespace("mgcvST"))
for (f in list.files(file.path(repo,"R"),pattern="\\.R$",full.names=TRUE)) sys.source(f,old)
data("MISO_E13",package="mgcvST")
basis <- spde_basis(MISO_E13$meshes$spde,as.matrix(MISO_E13$covariates[,c("x","y")]),kappa=.1)
model <- model.set(response ~ offset(offset0),MISO_E13$covariates,basis,family=mgcv::nb())
# Twelve feature indices, repeating the three bundled genes four times.
Y <- t(as.matrix(MISO_E13$expression))[rep(1:3,4),]
rownames(Y) <- paste0("feature",seq_len(nrow(Y)))
fit <- mgcvST.estimate(Y,model)
pairs <- t(combn(seq_len(nrow(Y)),2))
strip <- function(z) {z$call <- z$timing <- NULL; z}
a <- old$mgcvST.test(fit,pairs=pairs)
b <- mgcvST.test(fit,pairs=pairs)
stopifnot(identical(strip(a),strip(b)))
ta <- tb <- numeric(5)
for (r in 1:5) {
  gc(FALSE); ta[r] <- system.time(old$mgcvST.test(fit,pairs=pairs))[["elapsed"]]
  gc(FALSE); tb[r] <- system.time(mgcvST.test(fit,pairs=pairs))[["elapsed"]]
}
result <- data.frame(engine="single_model",calibration="liu",observations=ncol(Y),
  features=nrow(Y),pairs=nrow(pairs),basis_dimension=ncol(basis$B),
  old_seconds=median(ta),new_seconds=median(tb),speedup=median(ta)/median(tb),
  statistical_output_identical=identical(strip(a),strip(b)))
print(result)
write.csv(result,file.path(out,"hot-path-benchmark.csv"),row.names=FALSE)
