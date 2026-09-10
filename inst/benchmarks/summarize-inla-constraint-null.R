# Summarize the finished, frozen type-I experiments without refitting.
out <- "artifacts/constraint-type1"
estimated <- read.csv(file.path(out,"estimated","type1-summary.csv"))
known <- read.csv(file.path(out,"gaussian-known","rejection-rates.csv"))
case_order <- read.csv(file.path(out,"estimated","design.csv"))$name
stopifnot(setequal(unique(estimated$case), case_order))
primary <- subset(estimated,alpha==.05 & variant %in% c("projected","raw_kernel_only"))
primary$estimate_with_ci <- sprintf("%.1f%% [%.1f, %.1f]",100*primary$rejection_rate,
                                   100*primary$ci_lower,100*primary$ci_upper)
write.csv(primary,file.path(out,"primary-comparison.csv"),row.names=FALSE)
print(primary[,c("case","variant","calibration","valid","failed","estimate_with_ci")],row.names=FALSE)

files <- list.files(file.path(out,"estimated"),pattern="-replicates[.]csv$",full.names=TRUE)
replicates <- do.call(rbind,lapply(files,read.csv))
diag <- data.frame(
  attempted_fitted_replicates=nrow(unique(replicates[c("case","replicate")])),
  maximum_fitted_mean_error=max(replicates$fitted_mean_error,na.rm=TRUE),
  maximum_equivalence_error=max(replicates$equivalence_error,na.rm=TRUE),
  minimum_working_weight=min(replicates$minimum_working_weight,na.rm=TRUE),
  failed_test_cells=sum(!is.finite(replicates$p_value)),
  fallback_cells=sum(replicates$fallback,na.rm=TRUE))
write.csv(diag,file.path(out,"estimated-diagnostics.csv"),row.names=FALSE)
print(diag)

fig <- subset(primary,calibration=="davies")
png(file.path(out,"estimated-type1.png"),width=1450,height=900,res=150)
par(mar=c(6,11,3,1))
limits <- range(c(0,.05,fig$ci_upper),finite=TRUE)
plot(NA,xlim=limits,ylim=c(.5,length(case_order)+.5),yaxt="n",ylab="",
     xlab="Rejection rate at nominal 5% (95% binomial interval)",
     main="INLA fitted null: conditioned SPDE kernel vs raw G")
axis(2,at=seq_along(case_order),labels=case_order,las=1,cex.axis=.8)
abline(v=.05,lty=2,col="gray35")
colors <- c(projected="#236AB9",raw_kernel_only="#C95626")
for (variant in names(colors)) {
  z <- fig[fig$variant==variant,]
  y <- match(z$case,case_order)+if(variant=="projected") -.10 else .10
  segments(z$ci_lower,y,z$ci_upper,y,col=colors[variant],lwd=2)
  points(z$rejection_rate,y,pch=19,col=colors[variant])
}
legend("topright",c("Current conditioned SPDE kernel","Raw G; P removes the intercept"),
       col=colors,pch=19,bty="n",cex=.8)
mtext("All fits mean-zero; log-hyperparameter N(0, 3^2). Davies fallbacks, if any, are tabulated separately.",
      side=1,line=5,cex=.65)
dev.off()

md_rows <- vapply(seq_len(nrow(primary)),function(i) {
  z <- primary[i,]
  sprintf("| %s | %s | %s | %s | %d/%d | %d |",z$case,z$variant,z$calibration,
          z$estimate_with_ci,z$valid,z$attempted,z$fallback)
},character(1L))
writeLines(c("| 场景 | 测试核 | 校准 | 拒绝率与95%区间 | 成功/尝试 | 回退 |",
             "| --- | --- | --- | --- | --- | --- |",md_rows),
           file.path(out,"estimated-table.md"))
