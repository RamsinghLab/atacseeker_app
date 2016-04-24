# fix this 
library(ATACseeker)
BAMs <- c(N1="~/BAMs/ATACseq/041814-1-N_TAAGGC_L005.hg19.bam",
          S1="~/BAMs/ATACseq/041814-1-S_CGTACT_L005.hg19.bam",
          N2="~/BAMs/ATACseq/041814-2-N_AGGCAG_L005.hg19.bam",
          S2="~/BAMs/ATACseq/041814-2-S_TCCTGA_L005.hg19.bam")

# someChrMReads <- extractReads(bam.file=BAMs[1],
#                               region=GRanges("chrM", IRanges(1, 16517)))
firstBp <- 1
lastBp <- 249250621
chr1_5primeCuts <- lapply(BAMs, function(BAM) 
                          resize(extractReads(bam.file=BAM,
                                              region=GRanges("chr1", 
                                                             IRanges(firstBp,
                                                                     lastBp))),
                           1, fix="start"))
ests <- lapply(chr1_5primeCuts, getEsts)

# normal HSPCs from a donor 
plotComplexity(chr1_5primeCuts$N1, ests=ests$N1)
# add senescent HSPCs from the same subject
plotComplexity(chr1_5primeCuts$S1, ests=ests$S1, add=TRUE) 
# normal HSPCs from another donor 
plotComplexity(chr1_5primeCuts$N2, ests=ests$N2, add=TRUE) 
# add senescent HSPCs from the other donor 
plotComplexity(chr1_5primeCuts$S2, ests=ests$S2, add=TRUE) 
