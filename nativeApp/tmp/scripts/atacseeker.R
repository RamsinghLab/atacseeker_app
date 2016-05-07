###############################################################
#
# Convenience functions for the ATACseeker Native App
#
###############################################################


## Wrapper for getPESizes function. 
getQC <- function(name, pe.bam) { 
  out <- getPESizes(pe.bam)
  too.large <- sum(out$sizes > 400)
  qc <- c(sample=name, out$diagnostics, too.large=too.large)
  frag.dist <- density(out$sizes[out$sizes <= 800])
  results <- list(qc=qc, frag.dist=frag.dist)
  return(results)
} 


## QC plot function for distirbution of fragment sizes. 
qcPlot <- function(frag.dists, label) { 
  dens.max <- max(do.call(c, lapply(frag.dists, `[[`, "y")))
  plot(frag.dists[[1]], lwd=3, ylim=c(0, dens.max * 1.1),
       xlab="Fragment sizes (base pairs)", ylab="Density", 
       main=paste(label, "sample fragment sizes"), col=1)
  for (i in seq_along(frag.dists)[-1]) {
    lines(frag.dists[[i]], col=i, lwd=3)
  }
  getMode <- function(dens) dens$x[which.max(dens$y)]
  frag.mode <- median(round(do.call(c, lapply(frag.dists, getMode))))
  abline(v=frag.mode, col="gray80", lwd=3, lty=3)
  legend("topright", legend=names(frag.dists), col=seq(length(frag.dists)), 
         lwd=1, bty="n")
}



## plot the top differentially accessible regions using Gviz
## plot the top differentially accessible regions using Gviz
plotRegion <- function(region) { 
    olReads <- lapply(bam.files, extractReads, region=region, param=discard.se.param)
    for (x in seq(num_samples)) 
        metadata(olReads[[x]])$totals <- colData(data[,x])$totals
    covs <- GRangesList(lapply(olReads, function(x) 
        as(coverage(x) / metadata(x)$totals, "GRanges")))
    maxDepth <- max(sapply(covs, function(x) max(score(x))))
    cols <- c("darkred","darkgreen")[seq_along(bam.files) %% 2 + 1]
    names(cols) <- samples
    collected <- list()
    for(i in bam.files) {
        covr <- covs[[i]]
        ctrack <- DataTrack(covr, type="histogram", lwd=0, fill=cols[i], name=i, 
                            ylim=c(0,maxDepth), col.axis="black", col.title="black")
        collected[[i]] <- OverlayTrack(trackList=list(coverage=ctrack))
    }
    gax <- GenomeAxisTrack(col="black")
    plotTracks(c(gax, collected), from=start(region), to=end(region))
    invisible(covs)
}

