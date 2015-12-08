###############################################################
#
# Convenience functions for the ATACseeker Native App
#
###############################################################


## Load bed file and generate GRanges file
bed2GRanges <- function(file){

  df <- read.table(file, header=F, stringsAsFactors=F)
  header <- c('chr','start','end','id','score','strand')
  names(df) <- header
  
  if('strand' %in% colnames(df)){
    df$strand <- gsub(pattern="[^+-]+", replacement = '*', x = df$strand)
  }
  
  gr <- makeGRangesFromDataFrame(df, starts.in.df.are.0based=TRUE)
  return(gr)
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



