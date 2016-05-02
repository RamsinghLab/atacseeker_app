library(Homo.sapiens)
proms <- promoters(Homo.sapiens, upstream=2000, downstream=2000) # for hg19

library(rtracklayer)
blacklist <- import("wgEncodeDacMapabilityConsensusExcludable.hg19.bed.gz")

# on BAMs
library(csaw)
bams <- list.files(pattern="bam$")
names(bams) <- bams # somewhat shady examples in my case... 

atac.param <- readParam(pe = "none", dedup=TRUE, minq=10, discard=blacklist,
                        restrict=paste0("chr", c(1:22, "X", "Y")))
wincounts <- windowCounts(bams, ext=0, shift=4, bin=TRUE, param=atac.param)
tssWindows <- queryHits(findOverlaps(wincounts, proms))
nonTssWindows <- setdiff(seq_len(nrow(wincounts)), tssWindows)

TSS.bamCounts <- colSums(assays(wincounts[tssWindows,])$counts)
nonTSS.bamCounts <- colSums(assays(wincounts[nonTssWindows,])$counts)
TSSvsNonTSS <- TSS.bamCounts / nonTSS.bamCounts # ideally > 10...
