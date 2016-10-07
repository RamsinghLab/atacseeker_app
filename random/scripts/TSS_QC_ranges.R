library(Homo.sapiens) 

# set up the disjoin() on 5k and 10k, then we will drop out the windows w/in 5k 
standardChroms <- paste0("chr", c(1:22,"X","Y"))
proms5k <- promoters(Homo.sapiens, up=5000, down=5000)
proms10k <- promoters(Homo.sapiens, up=10000, down=10000)

# TSS distal (5kb-10kb)
tss5k10k <- keepSeqlevels(trim(disjoin(c(proms10k, proms5k))), standardChroms)
tss5k10k <- unique(tss5k10k[-queryHits(findOverlaps(tss5k10k, proms5k))])
stopifnot(length(findOverlaps(tss5k10k, proms5k)) == 0) # ensure we got 'em all

# TSS proximal
tss2k <- unique(keepSeqlevels(trim(promoters(Homo.sapiens,up=2000,down=2000)),
                              standardChroms))

# test the rest, I can't right now; should be correct though  
wincounts <- windowCounts(bam.files, bin=TRUE, param=discard.se.param)
tssWindows <- queryHits(findOverlaps(wincounts, tss2k))
nonTssWindows <- queryHits(findOverlaps(wincounts, tss5k10k))
TSS.bamCounts <- colSums(assays(wincounts[tssWindows,])$counts)
nonTSS.bamCounts <- colSums(assays(wincounts[nonTssWindows,])$counts)
TSSvsNonTSS <- TSS.bamCounts / nonTSS.bamCounts

