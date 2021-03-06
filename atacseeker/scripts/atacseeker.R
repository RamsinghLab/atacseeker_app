###############################################################
#
# Convenience functions for the ATACseeker Native App
#
###############################################################


## copy folder with files
dir.copy <- function(folder, to){
    dir.create(file.path(to, folder), recursive = T)
    fs = list.files(folder, full.names=T, recursive = T)
    file.copy(from=fs, to=file.path(to, folder), recursive = F, copy.mode = T)
}


## garbage clean up files
cleanMem <- function(n=10) { for (i in 1:n) gc(T,T) }
collectGarbage <- function() {while (gc()[2,4] != gc()[2,4] | gc()[1,4] != gc()[1,4]){}}


## get bam files specs
get_bams <- function(IDs, bam_names){
    group.in = c()
    group.files = c()
    for (f in seq(length(bam_names))){
        fastq_dir = file.path(base_dir, IDs[f])
        bam.tmp = list.files(fastq_dir, "*bam$", full.names = TRUE)
        bam.dmp = file.path(dump_dir, basename(bam.tmp))
        group.in = c(group.in, bam.tmp)
        group.files = c(group.files, bam.dmp)
    }
    return(list(group.in = group.in, group.files = group.files))
}


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
plotRegion <- function(region) { 
    olReads <- lapply(bam.files, extractReads, region=region, param=discard.se.param)
    for (x in seq(num_samples)) 
        metadata(olReads[[x]])$totals <- colData(data[,x])$totals
    covs <- GRangesList(lapply(olReads, function(x) 
        as(coverage(x) *1e6/ metadata(x)$totals, "GRanges")))
    maxDepth <- max(sapply(covs, function(x) max(score(x))))
    ##cols <- c("darkred","darkgreen")[seq_along(bam.files) %% 2 + 1]
    cols <- c(rep("blue", length(control)), rep("red",length(compare)))
    names(cols) <- samples
    collected <- list()
    for(i in seq(num_samples)) {
        covr <- covs[[i]]
        ctrack <- DataTrack(covr, type="histogram", lwd=0, fill=cols[i], name=samples[i], 
                            ylim=c(0,maxDepth), col.axis="black", col.title="black")
        collected[[i]] <- OverlayTrack(trackList=list(coverage=ctrack))
    }
    gax <- GenomeAxisTrack(col="black")
    plotTracks(c(gax, collected), from=start(region), to=end(region), chromosome=as.character(seqnames(region)))
    invisible(covs)
}

