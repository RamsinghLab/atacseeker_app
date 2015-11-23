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

# now dump a diagnostic table 
combinedQC <- append(controlQC, compareQC)
readqc <- as.data.frame(do.call(rbind, lapply(combinedQC, `[[`, "qc")))
colnames(readqc) <- c("Sample", "Total Reads", "Mapped Reads", "Single Read", "Mate Unmapped", "Unoriented", "Inter Chromosomal", "Frag. Size > 400 bp")
kable(readqc, caption="Quality metrics for PE reads")
