library(csaw) 
library(edgeR)
library(statmod)
library(locfit)
library(knitr)
library(ATACseeker)

# divergent from app 
base_dir <- "~/baseSpaceBAMs/CohesinATACSeq/"
setwd(base_dir)

control.label <- "wildtype"
compare.label <- "cohesinMutant"

# function is custom for this particular ATAC analysis.
getName <- function(bam, s="-") strsplit(bam, s)[[1]][2]

control.files <- list.files(pattern="WT.*.bam$")
names(control.files) <- sapply(control.files, getName)
controls <- names(control.files)

compare.files <- list.files(pattern="MUT.*.bam$")
names(compare.files) <- sapply(compare.files, getName)
compares <- names(compare.files)

bam.files <- append(control.files, compare.files)
samples <- c(controls, compares)

# function is custom for this particular ATAC analysis.
getStatus <- function(name, s="_") strsplit(name, s)[[1]][3]
getSubject <- function(name, s="_") strsplit(name, s)[[1]][1]

treatment <- as.factor(sapply(samples, getStatus))
treatment <- relevel(treatment, which(levels(treatment) == "WT"))
subject <- as.factor(sapply(samples, getSubject))
design <- model.matrix(~ treatment + subject)
rownames(design) <- samples

# for a two-factor design matrix...
kable(design, caption = "Design Matrix")

# estimate from data??
frag.len <- 150

# Do we use the minq=50 parameter ? (I don't)
pe.first <- readParam(dedup=TRUE, pe="first")
pe.param <- readParam(max.frag = frag.len, pe = "both", dedup=TRUE)

# parallelize:
library(parallel)
library(ATACseeker) # for getQC

# use mcmapply(), but only if we haven't already... lots and lots of bitching
if (!file.exists("cohesinControlQC.rds")) {
  controlQC <- mcmapply(getQC, name=controls, pe.bam=control.files, SIMPLIFY=F)
  saveRDS(controlQC, "cohesinControlQC.rds")
} else { 
  controlQC <- readRDS("cohesinControlQC.rds")
}
if (!file.exists("cohesinCompareQC.rds")) {
  compareQC <- mcmapply(getQC, name=compares, pe.bam=compare.files, SIMPLIFY=F)
  saveRDS(compareQC, file="cohesinCompareQC.rds")  
} else { 
  compareQC <- readRDS("cohesinCompareQC.rds")
}

# use ATACseeker::fragPlot
library(ATACseeker)

# side by side 
par(mfrow=c(1,2))

# plot fragment sizes for control samples
frag.dists.control <- lapply(controlQC, `[[`, "frag.dist")
fragPlot(frag.dists.control, control.label)

# plot fragment sizes for comparison samples
frag.dists.compare <- lapply(compareQC, `[[`, "frag.dist")
fragPlot(frag.dists.compare, compare.label)

dev.copy2pdf(file="cohesinFragPlots.pdf")

# now dump a diagnostic table 
combinedQC <- append(controlQC, compareQC)
readqc <- as.data.frame(do.call(rbind, lapply(combinedQC, `[[`, "qc")))
colnames(readqc) <- c("Sample", "Total Reads", "Mapped Reads", "Single Read", "Mate Unmapped", "Unoriented", "Inter Chromosomal", "Frag. Size > 400 bp")
kable(readqc, caption="Quality metrics for PE reads")

# cross correlate 
max.delay <- 400 # why?
x <- correlateReads(bam.files, max.delay, param=pe.first)
plot(0:max.delay, x, type="l", ylab="CCF", xlab="Delay (bp)")
title("Cross Correlation Plot")
saveRDS(x, file="correlateReadsOutput.rds")

# coverage (parallelize again?)
pe.param <- readParam(max.frag = frag.len, pe = "both", dedup=TRUE)

# added to ATACseeker as ATACseeker::collectStats
collectStats <- function(curbam, pe.param) { # {{{
  wcnt <- windowCounts(curbam, spacing=50, width=50, param=pe.param, filter=20)
  rsums <- rowSums(assays(wcnt)$counts)
  maxed <- findMaxima(rowRanges(wcnt), range=1000, metric=rsums)
  rmax <- rowRanges(wcnt)[maxed]
  weight <- 1 / rsums[maxed]
  profileSites(curbam, rmax, param=pe.param, weight=weight)
} # }}}
collected <- mclapply(bam.files, collectStats, pe.param=pe.param)
saveRDS(collected, file="collectedStats.rds")

# coverage plot at various window sizes
xranged <- as.integer(names(collected[[1]]))
for (i in seq(length(samples))){
  if (i == 1) {
    plot(xranged, collected[[i]], type="l", col=i, lwd=1,
     xlim=c(-1000, 1000), ylim=c(0, 1.0), 
     xlab="Distance (bp)", ylab="Relative coverage per base")
  } else {
    lines(xranged, collected[[i]], col=i, lwd=1)
  }
}
legend("topright", col=seq(length(samples)), samples, pch=16)
title("Coverage Profile")
abline(v=c(-150,150), col="gray", lty=3, lwd=4)
dev.copy2pdf(file="cohesinCoveragePlots.pdf")

# window 
window.width <- 75
data <- windowCounts(bam.files, width = window.width, param = pe.param)
kable(as.data.frame(data$totals, row.name=samples), 
      caption="Total number of reads in each library")
binned <- windowCounts(bam.files, bin=TRUE, width=1000, param=pe.param)

# filter 
filter.stat <- filterWindows(data, background=binned, type="global")
keep <- filter.stat$filter > log2(3) # what's this? compute from data?

# FIXME: what about library complexity based on 5' cut sites? 
#        it would be nice to have that back, perhaps also useful in QC!

# FIXME: also, what about TSS +/- 2000bp plots for ATACseq QC?
#        this is different from the usual ChIPseq QC but important

# one plot-per-window from here on out
par(mfrow=c(1,1))
hist(filter.stat$back.abundances, xlab="Adjusted bin log-CPM", breaks=100, 
     main="", xlim=c(min(filter.stat$back.abundances), 0))
global.bg <- filter.stat$abundances - filter.stat$filter
abline(v=global.bg[1], col="red", lwd=2)
abline(v=global.bg[1] + log2(3), col="blue", lwd=2)
legend("topright", lwd=2, col=c("red", "blue"), 
       legend=c("Background", "Threshold"))
dev.copy2pdf(file="cohesinATACthresholdPlot.pdf")
filtered.data <- data[keep,]
saveRDS(filtered.data, file="filtered.data.rds")
kable(as.data.frame(filtered.data$totals, row.name=samples), 
      caption="Total number of reads post filtering.")  # senescent2 is fucked

# Normalization factors
Norm.factors <- normOffsets(binned)
saveRDS(Norm.factors, file="Norm.factors.rds") 
kable(as.data.frame(Norm.factors, row.names=samples), 
      caption="Normalizing factors")
par(mfrow=c(1,2))
plot(Norm.factors, ylab="Normalization factor", xlab="Sample", 
     main="Count normalization factor, by sample", pch=16,
     col=ifelse(Norm.factors < 2, "darkgreen", "red"))
barplot(filtered.data$totals, names.arg=seq_along(samples), xlab="sample",
        col="gray80", ylab="Reads left after filtering", 
        main="Filtered read counts, by sample")
dev.copy2pdf(file="NormFactors.pdf")
par(mfrow=c(1,1))

# MDS plot -- any particular reason to leave the width at 1000?
binned.2 <- windowCounts(bam.files, bin=TRUE, width=1000, param=pe.param)
bin.adjc <- cpm(asDGEList(binned), log=TRUE)
saveRDS(bin.adjc, file="bin.adjc.rds")
plotMDS(bin.adjc, labels=samples, main="Multidimensional scaling plot (for QC)",
        col=c(rep("darkblue",length(controls)),rep("darkred",length(compares))))
dev.copy2pdf(file="Cohesin_ATAC_MDSplot.pdf")

# differential accessibility testing
y <- asDGEList(filtered.data, norm.factors=Norm.factors)
y <- estimateDisp(y, design) # trended/common/robust? need to check 
saveRDS(y, file="y.rds")

par(mfrow=c(1,2))
o <- order(y$AveLogCPM)
plot(y$AveLogCPM[o], sqrt(y$trended.dispersion[o]), type="l", lwd=2,
     ylim=c(0,1), xlab="Average log2 CPM", 
     ylab="Biological coefficient of variation")
fit <- glmQLFit(y, design, robust=TRUE, coef=2)
plotQLDisp(fit)
dev.copy2pdf(file="QLdisp.pdf")

# uncorrected for FDR: 
results <- glmQLFTest(fit, coef=2)
kable(topTags(results), caption="Top Tags")

# corrected: 
merged <- mergeWindows(rowRanges(filtered.data), tol=1000L)
tabcom <- combineTests(merged$id, results$table)
kable(head(tabcom), caption="Combined Tests")

# annotated:
library(Homo.sapiens)
#is.sig <- tabcom$FDR <= 0.05
anno <- detailRanges(merged$region, txdb=Homo.sapiens, orgdb=Homo.sapiens)
combined <- data.frame(as.data.frame(merged$region)[,1:3], tabcom, anno)
diffAccRanges <- makeGRangesFromDataFrame(combined, keep=TRUE)

# export to a BED file
library(rtracklayer)

# -log10(p) as bed file "score", could also export to bigWig
score(diffAccRanges) <- -1 * log10(diffAccRanges$PValue)

# Output sites as a BED file (post-FDR-correction, chop at 0.1 or 0.01?) 
export(diffAccRanges, "diffAccRanges.CohesinMutants.hg19.bed")

# print out the significant (FDR <= 0.1) results 
kable(as(diffAccRanges, "data.frame"))

# run LOLA? Yes!
library(LOLA) 
dbPath <- system.file("extdata", "hg19", package="LOLA")
regionDB <- loadRegionDB(dbPath)
data("sample_input", package="LOLA") # load userSets
data("sample_universe", package="LOLA") # load userUniverse
if (FALSE) show(userUniverse)
locResults <- runLOLA(diffAccRanges, userUniverse, regionDB, cores=1)
kable(locResults[order(qValue),]) # this is just a demo, fix for app

if(FALSE) {

  # for later
  library(rtracklayer)
  CTCFsites <- import("ubiquitousCTCFsites.hg19.bed")
  repeatRegions <- import("repeatmasker.hg19.bed")

  # FIXME: use erma for this 
  hspcChromHMM <- import("CD34.chromImpute.hg19.bed")
}

# Motif enrichment (PWM-based, will re-use this for differential ChIP-seq app)
library(PWMEnrich)
library(PWMEnrich.Hsapiens.background)
library(BSgenome.Hsapiens.UCSC.hg19)
data(PWMLogn.hg19.MotifDb.Hsap)
seqs <- getSeq(Hsapiens, diffAccRanges)
res <- motifEnrichment(seqs, PWMLogn.hg19.MotifDb.Hsap)
saveRDS(res, file="motifEnrichment.cohesinMutants.hg19.rds")

# Top 10 motifs by enrichment
top10p <- head(motifRankingForGroup(res), 10)
top10 <- data.frame(motifName=names(top10p), pvalue=top10p)
kable(top10, digits=20, 
      caption="Top 10 motifs within differentially accessible regions")

# Dump session info for versioning, debugging, CITATIONS (!), etc.
sessionInfo()

# Explicitly remind user: cite us and/or the various components of the pipeline!
message("Please do not forget to appropriately cite this work!")
citation("ATACseekeR")
citation("csaw")
citation("PWMEnrich")
citation("LOLA")
citation()
