library(csaw)
library(Gviz)
library(rtracklayer)
diffAccessible <- import("diffAccRanges.NS.fdrLt01.hg19.bed", genome="hg19")
bigDiff <- diffAccessible[max(diffAccessible$score)]
start(bigDiff) <- start(bigDiff) - 100
end(bigDiff) <- end(bigDiff) + 100

setwd("~/baseSpaceBAMs/ATACseq/")
BAMs <- lf(patt="bam$")
names(BAMs) <- BAMs

frag.len <- 150
param <- readParam(max.frag = frag.len, pe = "both", dedup=TRUE)
data <- windowCounts(BAMs, width=75, param=param)
saveRDS(data, file="data.rds")

plotRegion <- function(region) { 
  olReads <- lapply(BAMs, extractReads, region=region, param=param)
  for (x in names(olReads)) 
    metadata(olReads[[x]])$totals <- colData(data[,x])$totals
  covs <- GRangesList(lapply(olReads, function(x) 
                             as(coverage(x) / metadata(x)$totals, "GRanges")))
  maxDepth <- max(sapply(covs, function(x) max(score(x))))
  cols <- c("darkred","darkgreen")[seq_along(BAMs) %% 2 + 1]
  names(cols) <- BAMs
  collected <- list()
  for(i in BAMs) {
    covr <- covs[[i]]
    ctrack <- DataTrack(covr, type="histogram", lwd=0, fill=cols[i], name=i, 
                        ylim=c(0,maxDepth), col.axis="black", col.title="black")
    collected[[i]] <- OverlayTrack(trackList=list(coverage=ctrack))
  }
  gax <- GenomeAxisTrack(col="black")
  plotTracks(c(gax, collected), from=start(region), to=end(region))
  invisible(covs)
}

# need to add this to discard.param in beginning of workflow
blacklist <- import("~/Dropbox/Asif/blacklist/wgEncodeDacMapabilityConsensusExcludable.hg19.bed.gz", genome="hg19")

# also restrict.param <- readParam(restrict=c(paste0("chr", 1:22), "chrX"))

accessByScore <- diffAccessible[rev(order(diffAccessible$score))]
plotRegion(accessByScore[1]) # blacklisted!
plotRegion(accessByScore[2]) # blacklisted! 
plotRegion(accessByScore[3]) # blacklisted! 
plotRegion(accessByScore[4])
plotRegion(accessByScore[5]) # blacklisted! 
plotRegion(accessByScore[6]) # blacklisted! 
plotRegion(accessByScore[7]) # blacklisted! 
plotRegion(accessByScore[8]) 
plotRegion(accessByScore[9]) # blacklisted! 
plotRegion(accessByScore[13]) 
plotRegion(accessByScore[15])

# running IGV from R to check on this:
library(SRAdb)
IGVsock <- IGVsocket()
IGVgoto(IGVsock, as.character(accessByScore[15]))

plotRegion(accessByScore[24]) 
IGVgoto(IGVsock, as.character(accessByScore[24]))

plotRegion(accessByScore[25])
IGVgoto(IGVsock, as.character(accessByScore[25]))
