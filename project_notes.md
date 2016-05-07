# ATACseeker App #

## mt-DNA Variant Calling ##

### mt-DNA Reference ###

This [page](http://haplogrep.uibk.ac.at/blog/rcrs-vs-rsrs-vs-hg19) has some information on the different references available for mtDNA

### ATACseq mt-DNA ###

Because in ATACseq a lot of the reads are biased for the mitochondrial DNA, it seems like an obvious thing to use ATACseq reads to call variants on mt-DNA. There are a few caveats though. The mt-DNA is haploid, inherited completely from the mother. However, most callers like `gatk` and `samtools` assume a diploid genome. 

This [post](http://gatkforums.broadinstitute.org/gatk/discussion/1214/can-i-use-gatk-on-non-diploid-organisms) says that post-version 3.3 `gatk` should be able to handle non-diploid genome. But a later [post](http://gatkforums.broadinstitute.org/gatk/discussion/3345/question-about-ploidy-when-mtdna-variants-calling) categorically denied any work by gatk relating to mt-DNA. Also, this later [post](http://gatkforums.broadinstitute.org/gatk/discussion/3651/gatk-for-mitochondrial-dna-mutations) claims that the user might have to play around with the ploidy settings of `gatk` but does not offer any conclusive answers. 

`samtools mpileup` explicitly says that it is for diploid genomes - as seen here:

```shell
Asifs-iMac:~ asifzubair$ samtools mpileup

Usage: samtools mpileup [options] in1.bam [in2.bam [...]]

  ...
  ...
  ...

Notes: Assuming diploid individuals.

```

However, a question that begs to be asked is that - if indeed we go out looking for heteroplasmy, shouldn't a heterozygous position by `samtools` give us evidence of the same. I fear my reasoning might be naive. 

In any case, I tried using samtools to call variants following the [CUREFFI](http://www.cureffi.org/2012/09/07/an-alternative-exome-sequencing-pipeline-using-bowtie2-and-samtools) folks. However, as I am interested in only reads that align with the mt-DNA, I filtered for those first. 

```shell
## Get chrM reads
for file in *bam; do  
    samtools view -b $file chrM > $(basename $file .bam)".chrM.bam"; 
done

## call variants
samtools mpileup -d 200 -D -B -f ../../tmp/reference/chrM.fa -b bamlist.txt -u | \
bcftools view -  > variants.vcf
```

However, when I looked at the VCF all calls were reported as missing. Also, in some cases `samtools` complained that the mate was missing. I then got counsel from this biostar [page](https://www.biostars.org/p/74386/) which says to use the `-r` option in place of filtering, which made my command even simpler.

```shell
samtools mpileup -d 200 -D -B -f ../../tmp/reference/hg19.fa -r chrM -b bamlist.txt -u | \
bcftools view -  > variants.vcf
```

Still, no cigar ! 

Based on a recommendation by [GGD](http://www.gettinggeneticsdone.com/2014/03/visualize-coverage-exome-targeted-ngs-bedtools.html), I checked that `bedtools` can be used for checking [coverage](https://github.com/arq5x/bedtools-protocols/blob/master/bedtools.md) and I tried to do the same for one ATACseq sample. I can't remember the exact command, but it was probably something like this.

```shell
# Calculate a histogram of coverage for each chromosome
    # as well as genome-wide.
    # Estimate: 30 minutes
    bedtools genomecov \
             -ibam NA19146.bam \
    > NA19146.coverage.hist.txt
```

I remember that when I ran the above command, per nucleotide coverage in chrM was quite high - above 100 bases. 

__4/13__: I finally got [MToolbox](https://github.com/mitoNGS/MToolBox) to work. The trick was getting the right indexes for `gsnap`. I used `gsnapv7` which is the same as the `MToolbox`'s page. A thing I had difficulty with is that my compiled version of `gsnap` will not work on all nodes on the cluster. 

To get all the mtDNA reads from aligned bams, I did the following:
```shell
samtools view -h atac.bam "chrM" > out.sam
samtools view -Sb out.sam > out.bam
```
[Biostar](https://www.biostars.org/p/56246) has some comments on using the `-f 4` flag with `samtools`. __must__ check with Tim.

## ATACseeker Pipeline ##

### preseqR ###

[preseqR](https://cran.r-project.org/web/packages/preseqR/index.html) is a software package, provided by the smith lab here at usc, that computes library complexity. it seems sensible to estimate complexity of your library before starting any analysis. 

before we get into `preseqR`, maybe it is better to figure out what the complexity beast is. 

#### library complexity ####

Library complexity refers to the number of unique fragments present in a given library.  

Complexity is affected by:
* Amount of starting material
* Amount of DNA lost during cleanups and size selection
* Amount of duplication introduced via PCR

This [page](https://www.kapabiosystems.com/ngs/guide-ngs-coverage-uniformity-bias-library-complexity) by `kapabiosystems` talks about library complexity from the experimental perspective. Mostly about how protocol change can help improve the complexity of the library. i found the section on - good starting amount of DNA but low coverage - most interesting.  

Even with good amount of starting dna, one can get uneven coverage because sequences with particular characteristics pass through each step of the sample prep and sequencing workflow with different efficiencies and some become relatively enriched at the expense of others, which are left behind.  

while a lot of experimental biases can not be overcome, suboptimal pcr during library amplification can introduce a lot of bias. also, base composition bias is another important consideration - GC vs. AT.  

this mit [lecture](http://ocw.mit.edu/courses/biology/7-91j-foundations-of-computational-and-systems-biology-spring-2014/video-lectures/lecture-5-library-complexity-and-short-read-alignment-mapping) also introduces complexity. it introduces the NB model for estimating library complexity. Also, how a simple poisson model maybe wrong, mostly because it doesn't capture the over-dispersion in the library.  

#### back to preseqR ####

with that very biref introduction, and more importantly for me - atleast a basic model understanding - let's see what `preseqR` actually does. paper is [here](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3612374/).  

seems like `preseqR` uses a [empirical Bayes](http://www.r-bloggers.com/understanding-empirical-bayes-estimation-using-baseball-statistics) approach to estimate complexity and is based on rational function approximation, [RFA](https://en.wikipedia.org/wiki/Pad%C3%A9_approximant), to the power series of [Good and Toulmin](http://biomet.oxfordjournals.org/content/43/1-2/45.abstract).

towards the more implementation side, tim D has a good [thread](http://seqanswers.com/forums/showthread.php?t=18439&goto=nextnewest) on `seqanswers` on `preseqR`. importantly, i might have to subsample reads before i do the library complexity.  

Also - it seems that `EstimateLibraryComplexity.jar` of [picard](http://broadinstitute.github.io/picard) also does complexity estimation. however, like tim D pointed out in the thread - it appears that `estimateLibrarySize` assumes a simple Lander-Waterman model, which would correspond to a simple Poisson model. The zero-truncated negative binomial (ZTNB) model is much broader class that includes the simple Poisson model (taking alpha -> 0). Therefore, the estimates from such a model can only be more biased than the ZTNB estimates.  

#### library complexity for atacseq data ####

something that tim T wants to do is look at 5' cut sites for library complexity for atacseq data. notably `preseR` references this [paper](http://www.nature.com/nmeth/journal/v9/n1/full/nmeth.1778.html) for identifying unique molecules. Also - this from tim D's paper - In sequencing applications that identify genomic intervals such as protein-binding sites in chromatin immunoprecipita- tion and sequencing (ChIP-seq) or expressed exons in RNA sequencing (RNA-seq), the number of distinct molecules in the library may be of secondary interest to the `number of distinct genomic intervals identified` after processing mapped reads.  

i feel this is a good direction for the pipeline, and i should figure out how to implement tim's idea. however, a thing that is troubling me that perhaps `preseqR` needs raw reads and the way the pipeline is set up right now is that it takes `bam` files.  

tim T said that the best thing to do is just use our [atacseeker](https://github.com/RamsinghLab/ATACseeker), which is great because it is in-house. i'm quite excited about this direction.  


### CSAW ###

`csaw` actually was written for ChIP-seq data but we appropriate it for atacseq data as the same assumptions apply.  

in it's basic form, `csaw` does read counting in windows along the chromosome and then uses `limma`'s negative-binomial model to find differential binding regions across experiments.  

we've re-formulated the probelm so that `csaw` now gives us differentially accessible regions based on atacseq data.  

#### QC ####

[aaron](https://github.com/LTLA) really did a good job putting together `csaw`. the documentation is great and the options carefully chosen.  

one of the first things about the atacseq data is that it has a lot of reads biased to the mt-DNA. thus, we restrict analysis to the autosomes and sex-chromosomes. also, we take out blacklisted regions as defined [here](https://sites.google.com/site/anshulkundaje/projects/blacklists).  

programmatically, this can be implemented in `csaw` as 
```R
readParam(pe = "none", restrict = chroms, discard = blacklist, minq = 10, dedup = TRUE)
``` 

an important thing here is that we have fixed `pe="none"`, effectively, treating the data as single-end irrespective of what type it is. i should clear this up with tim, but the idea here is that read extension might give us spurious results especially when a large fragment, which might be straddled by nucleosomes, is sequenced. while the ends of the fragment are accessible, it is **wrong** if we do read extension and think that the region with nucleosomes is also accessible.

we also, filter for bad quality reads and **remove duplicates**.  

aaron claims that removing dupes is not as straightforward as we think. per aaron, for the ChIP-seq case, it also caps the number of reads at each position. This can lead to loss of DB detection power in high abundance regions. Spurious differences may also be introduced when the same upper bound is applied to libraries of varying size. Thus, duplicate removal is not recommended for routine DB analyses. Of course, it may be unavoidable in some cases, e.g., involving libraries generated from low quantities of DNA. i **need** to figure out if this applies to atacseq data as well.  

the first QC metrics we compute is the plot of `TSSvsNonTSS` counts. if i remember correctly, tim T said that this should be above 10 for good quality data. 

`csaw`'s `windowCounts` is the main function which does read counting over chromosomes. it does some read extension for forward and reverse reads to account for fragment length. for the atacseq data, i don't want to do any read extension and i think i need to set `ext` option to `NA` (it defaults to a 100). another thing i want to add to the pipeline output is the options used for `windowCounts`, something like `metadata(data)` should do this for me. this is cool stuff - read extension was bothering me for some time now.  

in light of the above, i think we need to change our `windowCounts` call from:

```R
windowCounts(bam.files, ext=0, shift=4, bin=TRUE, param=discard.se.param)
```

to:

```R
windowCounts(bam.files, ext=NA, shift=4, bin=TRUE, param=discard.se.param)
```

the `shift` parameter is for `specifying how much the start of each window should be shifted to the left`.By default, the first window on a chromosome starts at base position 1. This can be shifted to the left by specifying an appropriate value for `shift`.  

however, I am not sure why Tim set it to 4 and i should ask him. 

also, when `bin` is TRUE, then the following flags are set:

```R
	# A convenience flag, which assigns sensible arguments to everything else.
	spacing <- as.integer(width)
	left <- as.integer(shift)
	right <- spacing - 1L - left
	ext <- 1L
	final.ext <- NA
	filter <- min(1, filter)
```

as we can see, `ext` is set to `1L`, `final.ext` is `NA` and also `filter` is set to `0/1`, `width` will default to 100. So I guess, even this should work:

```R
windowCounts(bam.files, shift=4, bin=TRUE, param=discard.se.param)
```

i guess onething, i **must** do is check `metadata` with these settings.  

Also, i got curious about `ext` and dug up the funciton to which `ext` is passed to get the `final.ext`:

```R
.collateExt <- function(nbam, ext)
# Collates the extension parameters into a set of ext and remainder values.
# The idea is to extend each read directionally to 'ext', and then extend in
# both directions by 'remainder' to reach the desired fragment length.

{
	final.ext <- attributes(ext)$final.ext # Do this, before attributes are lost.
	if (is.null(final.ext)) { final.ext <- NA }
	final.ext <- as.integer(final.ext)
	if (length(final.ext)!=1L || (!is.na(final.ext) && final.ext <= 0L)) { 
		stop("final extension length must be a positive integer or NA") 
	}
	
	if (length(ext)==1L) { 
		ext <- rep(ext, nbam)
	} else if (length(ext)!=nbam) {
		stop("length of extension vector is not consistent with number of libraries")
	}
	ext <- as.integer(ext)
	if (any(!is.na(ext) & ext <= 0L)) { stop("extension length must be NA or a positive integer") }

	list(ext=ext, final=final.ext)
}
```

i think i will stick with setting `ext=NA` but less sure about how to interpret it now.  

i used `Rsamtools`'s `testPairedEndBam` function to finding out if the user data is PE and plot fragment size density. While I am not using fragment sizes for the pipeline, the idea here is that there should be within group similarity and outside group dissimilarity for the distribution. if this is not the case, something might be weird - who knows - library swap (?) etc.  

__cross correlation plot__  

`csaw` cites this [paper](http://www.nature.com/nbt/journal/v26/n12/full/nbt.1508.html) for the cross-correlation plot. For ChIP-seq data, this plot provide a measure of the immunoprecipitation (IP) efficiency of a ChIP-seq experiment. Efficient IP should yield a smooth peak at a delay distance corresponding to the average fragment length.

For atacseq data, however, the fragmet sizes are greatly varied and so I don't expect this plot to be smooth. Interestingly, for the truncated data that I used for the initial analysis - i see two well-separated sharp peaks.

A sharp spike may also observed in the plot at a distance corresponding to the read length. This is thought to be an **artifact**, caused by the preference of aligners towards uniquely mapped reads. Duplicate removal is typically required here (`dedup=TRUE` in `readParam`) to reduce the size of this spike. Otherwise, the fragment length peak will not be visible as a separate entity. The size of the smooth peak can also be compared to the height of the spike to [assess the signal-to-noise ratio](http://genome.cshlp.org/content/22/9/1813.long) of the data. Poor IP efficiency will result in a smaller or absent peak as bimodality is less pronounced. However, still has to be assessed what this will mena for atacseq data.  

`csaw` also supports paired-end data, whereby correlations are computed using only those reads in proper pairs. This may be less meaningful as the presence of proper pairs will inevitably result in a strong peak at the fragment length. Instead, IP efficiency can be diagnosed by treating paired-end data as single-end, e.g., with `pe="first"` in `readParam`.

currently, i have `pe="none"`. confused whether `pe="both"` will give me back the fragment-size distribution plot ?

__coverage plot__

The coverage profile around potential binding sites can be obtained with the `profileSites` function. Here, the binding sites are defined by taking high-abundance `50 bp` windows and identifying those that are __locally maximal__ using `findMaxima`. For each selected window, `profileSites` records the coverage across the flanking regions as a function of the distance from the edge of the window. This is divided by the count for the window itself to obtain a relative coverage, based on the specification of `weight`. The values are then averaged across all windows to obtain an aggregated coverage profile for each library.

The version of this plot in the manual for ChIP-seq data seems to be quite smooth. however, my plot is a little jagged and i am not sure if this is due to truncated data. 

also, the `windowCounts` function curently looks like this:
```R
windowCounts(curbam, spacing=150, width=150, param=discard.se.param, filter=20, ext=1)
```

i've set the `spacing` and width set to 150 but I think I should set it to 50. also, `ext=NA` should be done. but, i am not sure if I should `bin` it. also, i am calling this function in a loop, but I think there should be a way to vectorise this. __must__ ask Tim. 

#### Filtering ####

Removing such uninteresting or ineffective tests reduces the severity of the multiple testing correction, increases detection power amongst the remaining tests and reduces computational work. Filtering is valid so long as [it is independent of the test statistic](http://www.pnas.org/content/107/21/9546.long) under the null hypothesis. In the negative binomial (NB) framework, this (probably) corresponds to filtering on the overall NB mean.

three approaches are suggested:
- by count size  
- by proportion  
	- here we assume that only a small proportion of the dna is accessible and then retain only those windows with rank ratios above the unbound proportion of the genome
- by global enrichment  
	- involves choosing a filter threshold based on the fold change over the level of non-specific enrichment. The degree of background enrichment can be estimated by counting reads into large bins across the genome.
	- The effect of filtering can also be visualized with a histogram. This allows users to confirm that the bulk of (assumed) background bins are discarded upon filtering. Note that bins containing genuine binding sites will usually __not__ be visible on such plots.

currently, i am using the global enrichment strategy to filter windows. However, as a general rule, I should filter __less aggressively__ if there is any uncertainty about the features of interest. In particular, the thresholds shown in this chapter for each filtering statistic are fairly mild.

#### Normalisation ####

For a Differential Accessibility analysis I want to do, library-specific biases are of particular interest as they can introduce spurious differences between conditions. This includes 
- composition biases
- efficiency biases
- trended biases  

Thus, normalization between libraries is required to remove these biases prior to any statistical analysis.

in the pipeline right now, i've implemented composition bias normalisation. Highly enriched regions consume more sequencing resources and thereby suppress the representation of other regions. Differences in the magnitude of suppression between libraries can lead to spurious DA calls. This is a typical normalisation that i have seen in rna-seq data as well.  

#### Clustering ####

`csaw` recommends generating the [MDS](https://en.wikipedia.org/wiki/Multidimensional_scaling) plot to investigate how the samples are clustered. Currently, I am using binned counts for making the plots. However, it seems that `csaw` doesn't do this in the manual.  

```R
> par(mfrow=c(2,2), mar=c(5,4,2,2))
> adj.counts <- cpm(y, log=TRUE)
> for (top in c(100, 500, 1000, 5000)) {
+ out <- plotMDS(adj.counts, main=top, col=c("blue", "blue", "red", "red"), 
+ labels=c("es.1", "es.2", "tn.1", "tn.2"), top=top)
+}
```
Also, I should use the `col` argument to label samples.  

It is a little confusing to me wheher I should use bins or not. `csaw` says this:  

> The distance between each pair of libraries is computed as the __square root of the mean squared log-fold change__ across the top set of bins with the highest absolute log-fold changes. A small top set visualizes the most extreme differences whereas a large set visualizes overall differences. Checking a range of top values may be useful when the scope of DB is unknown. Again, __counting with large bins is recommended__ as fold changes will be undefined in the presence of zero counts.

#### Differential Accesibility Testing ####

`csaw` uses the [quasi-likelihood framework](http://scripts.mit.edu/~pazoulay/docs/wooldridge_QML.pdf) of  `edgeR` for testing DB. NB distribution is used to model counts and the usual argument of NB distributions accounting for overdispersion is given.  

quasi-likelihood estimation is one way of allowing for __overdispersion__, that is, greater variability in the data than would be expected from the statistical model used. It is most often used with models for count data or grouped binary data, i.e. data that otherwise be modelled using the Poisson or binomial distribution.  

Instead of specifying a probability distribution for the data, only a __relationship between the mean and the variance__ is specified in the form of a __variance function__ giving the variance as a function of the mean. Generally, this function is allowed to include a multiplicative factor known as the __overdispersion parameter or scale parameter__ that is estimated from the data. Most commonly, the variance function is of a form such that __fixing the overdispersion parameter at unity results in the variance-mean relationship of an actual probability distribution such as the binomial or__ [Poisson](https://en.wikipedia.org/wiki/Generalized_linear_model#Count_data).

__Stabilising estimates with empirical Bayes__

Under the QL framework, both the QL and NB dispersions are used to [model biological variability in the data](http://www.statsci.org/smyth/pubs/QuasiSeqPreprint.pdf). The former ensures that the NB mean-variance relationship is properly specified with appropriate contributions from the Poisson and Gamma components. The latter accounts for variability and uncertainty in the dispersion estimate.  

The effect of EB stabilisation can be visualized by examining the biological coefficient of variation (for the NB dispersion) and the quarter-root deviance (for the QL dispersion). These plots can also be used to __decide whether the fitted trend is appropriate__. Sudden irregulaties may be indicative of an __underlying structure in the data__ which cannot be modelled with the mean-dispersion trend. Discrete patterns in the raw dispersions are indicative of low counts and suggest that __more aggressive filtering__ is required.

__Important:__ A strong trend may also be observed where the NB dispersion drops sharply with increasing average abundance. It is difficult to accurately fit an empirical curve to these strong trends. As a consequence, the dispersions at high abundances may be overestimated. __Filtering__ of low-abundance regions provides some protection by removing the strongest part of the trend.  

Users can compare raw and filtered results to see whether it makes any difference. Example code below:  
```R
> relevant <- rowSums(assay(data)) >= 20 # some filtering; otherwise, it takes too long.
> yo <- asDGEList(data[relevant], norm.factors=normfacs)
> yo <- estimateDisp(yo, design)
> oo <- order(yo$AveLogCPM)
> plot(yo$AveLogCPM[oo], sqrt(yo$trended.dispersion[oo]), type="l", lwd=2,
+ ylim=c(0, max(sqrt(yo$trended))), xlab=expression("Ave."~Log[2]~"CPM"),
+   ylab=("Biological coefficient of variation"))
> lines(y$AveLogCPM[o], sqrt(y$trended[o]), lwd=2, col="grey")
> legend("topright", c("raw", "filtered"), col=c("black", "grey"), lwd=2)
```

__Note:__ Rafa had a good [post](http://simplystatistics.org/2014/10/13/as-an-applied-statistician-i-find-the-frequentists-versus-bayesians-debate-completely-inconsequential) on bayesian vs. frequentist approaches where he briefly mentioned empirical Bayes approach. It might also be instructive to go over this [r-bloggers](http://www.r-bloggers.com/understanding-empirical-bayes-estimation-using-baseball-statistics) example.  

#### Multiple Testing Correction ####



### LOLA ###



## Docker ##


## Native App ##


### testing basespace app ###

`ilmn` has an  interesting way of testing these apps. essentially, i downloaded the native app development [image](https://da1s119xsxmu0.cloudfront.net/sites/developer/native/nativeappsvm/BaseSpace%20Native%20Apps%20VM%20v14.ova) and then loaded it in virtual box. You need to login to the image via the terminal - `ssh basespace@localhost -p2222` . The password is `basespace`. 

This image contains all software to test native apps - docker, spacedock etc. Also the versions should be up-to-date with the ones `ilmn` has. 

Associated with every app, is a `spacedock` command which in my case is - `sudo spacedock -a d31c319955504b5ba162a51a0026d693 -m https://mission.basespace.illumina.com` . The idea is to run this command locally inside the VM and then trigger a test by doing a `Send To Local Agent`.

My first attempt wasn't so successful. I checked the `spacedock` log and saw this relevant part. 

```shell
2016-02-03 22:46:13.719 [THREAD Threadpool worker] [DEBUG] Illumina.SpaceDock.JobExecutorLogic                Error during JobExecutionStage DOWNLOADING
Illumina.BaseSpace.SDK.BaseSpaceException: GET:v1pre3/appresults/28995055 status: 410 (Gone) Message: This resource and its underlying data has been removed by its owner (BASESPACE.COMMON.RESOURCE_REMOVED) ---> ServiceStack.ServiceClient.Web.WebServiceException: Gone
  at ServiceStack.ServiceClient.Web.ServiceClientBase.ThrowWebServiceException[ErrorResponse] (System.Exception ex, System.String requestUri) [0x00000] in <filename unknown>:0 
  at ServiceStack.ServiceClient.Web.ServiceClientBase.ThrowResponseTypeException[GetAppResultResponse] (System.Object request, System.Exception ex, System.String requestUri) [0x00000] in <filename unknown>:0 
  at ServiceStack.ServiceClient.Web.ServiceClientBase.HandleResponseException[GetAppResultResponse] (System.Exception ex, System.Object request, System.String requestUri, System.Func`1 createWebRequest, System.Func`2 getResponse, Illumina.BaseSpace.SDK.ServiceModels.GetAppResultResponse& response) [0x00000] in <filename unknown>:0 
  at ServiceStack.ServiceClient.Web.ServiceClientBase.Send[GetAppResultResponse] (System.String httpMethod, System.String relativeOrAbsoluteUrl, System.Object request) [0x00000] in <filename unknown>:0 
  at Illumina.BaseSpace.SDK.ServiceModels.AbstractRequest`1+<>c__DisplayClass1[Illumina.BaseSpace.SDK.ServiceModels.GetAppResultResponse].<GetSendFunc>b__0 () [0x00000] in <filename unknown>:0 
  at Illumina.BaseSpace.SDK.JsonWebClient+<>c__DisplayClass5`1[Illumina.BaseSpace.SDK.ServiceModels.GetAppResultResponse].<Send>b__2 () [0x00000] in <filename unknown>:0 
  at Illumina.BaseSpace.SDK.RetryLogic.DoWithRetry (UInt32 maxAttempts, System.String description, System.Action op, ILog logger, Double retryIntervalBaseSecs, System.Action error, System.Func`2 retryHandler) [0x00000] in <filename unknown>:0 
  --- End of inner exception stack trace ---
  at Illumina.BaseSpace.SDK.JsonWebClient.Send[GetAppResultResponse] (Illumina.BaseSpace.SDK.ServiceModels.AbstractRequest`1 request, IRequestOptions options) [0x00000] in <filename unknown>:0 
  at Illumina.BaseSpace.SDK.BaseSpaceClient.GetAppResult (Illumina.BaseSpace.SDK.ServiceModels.GetAppResultRequest request, IRequestOptions options) [0x00000] in <filename unknown>:0 
  at Illumina.SpaceDock.Core.BaseSpaceLogic+<>c__DisplayClass8.<CreateAppResultFile>b__7 (IBaseSpaceClient client) [0x00000] in <filename unknown>:0 
  at Illumina.SpaceDock.Core.BaseSpaceLogic.GetPropertyLimits (System.Uri apiUri, System.String accessToken, System.Func`2 request) [0x00000] in <filename unknown>:0 
  at Illumina.SpaceDock.Core.BaseSpaceLogic.CreateAppResultFile (System.Uri apiServerUri, System.String appResultId, System.String filePath, System.String accessToken) [0x00000] in <filename unknown>:0 
  at Illumina.SpaceDock.JobExecutorLogic.CreateAppResultJsons (System.Collections.Generic.List`1 inputAppResults, System.Uri apiUri, System.String accessToken) [0x00000] in <filename unknown>:0 
  at Illumina.SpaceDock.JobExecutorLogic.OnDownloading () [0x00000] in <filename unknown>:0 
  at Illumina.SpaceDock.JobExecutor+<>c__DisplayClass2.<ExecuteJob>b__0 () [0x00000] in <filename unknown>:0 
  at Illumina.SpaceDock.JobExecutorLogic.ExecuteAndCleanUp (System.Action jobWorkflow) [0x00000] in <filename unknown>:0 
2016-02-03 22:46:13.721 [THREAD Threadpool worker] [DEBUG] Illumina.BaseSpace.SDK.JsonWebClient               POST:v1pre3/appsessions/32462305/properties
```

I think this has to be related to the app results I chose. I might not have permissions for them. 

When I used only bams from the BWA Aligner app I had the following errors:

```shell
  ...
  ...
2[bam_index_build2] fail to create the index file.
2[bam_index_build2] fail to create the index file.
2[bam_index_build2] fail to create the index file.
2[bam_index_build2] fail to create the index file.
  ...
  ...
  ordinary text without R code

K
  |                                                                       K
  |....................                                             |  31%
..label:.. ..unnamed-chunk-6 (with options) 
List of 1
 $ echo: logi FALSE

.Quitting from lines 125-136 (atacseeker.Rmd) 
Error in plot.window(xlim, ylim, log = log, ...) : 
  need finite 'ylim' values
Calls: <Anonymous> ... eval -> eval -> barplot -> barplot.default -> plot.window

Execution halted
cp: 1cannot stat '/atacseeker/scripts/atacseeker.html': No such file or directory
```

So, it seems that creating index is failing. This might be because I have to write everything to scratch - which is pretty agonizing. 

I spoke to Anthony and my hunch was right. When `bsfs` is enabled, the `spacedock` directories are read only. This means I'll have to copy everything to `/data/scratch` and build my indexes there.

After I fixed the `/data/scratch` issue, I ran into this error:

```R
Quitting from lines 138-176 (atacseeker.Rmd) 
Error in mcmapply(getQC, name = control, pe.bam = control.files, SIMPLIFY = F) : 
  'names' attribute [2] must be the same length as the vector [1]
Calls: <Anonymous> ... withCallingHandlers -> withVisible -> eval -> eval -> mcmapply
```
 
... which is weird, because I am not seeing this when I run the app locally. My hunch here is that this could be a memory issue with `mcmapply`. 

So, my hunch was correct! Since these jobs are sent to my computer to be run locally, I could increase the memory available to the VM and that should facilitate job completion ... and it did! 

### report builder ###

`ilmn` uses report builder which relies on liquid/javascript etc. but quite frankly fails to impress. I thought of publishing my report as a `html` document and then using iframes for serving this to the user. 

```js
{% comment %}
  BaseSpace Report Builder!
  <td>{{ result.files["atacseeker.html"] | append: "bar" }}</td>
{% endcomment %}

{% for key in result.files %}
  <iframe src={{ result.files[key].href }} frameborder="0" style="overflow:hidden;height:100%;width:100%" height="100%" width="100%">
  </iframe>
{% endfor %}
```

This hack works pretty well, especially after `ilmn` changed their layout. However, I found a problem with the links I put in my report. I want those to open in a new window, otherwise it goes into unending loading phase. [SO](http://stackoverflow.com/questions/3492153/markdown-open-a-new-window-link), like always, has some ideas about that. 

Essentially, since `Rmarkdown` can parse `html` also, I could simply use `html` syntax - `<a href="http://example.com/" target="_blank">example</a>` or this hybrid syntax that I am not too sure about `[link](url){:target="_blank"}` 

A final idea for the analysis report is that it can be useful to researchers to use in grant applications/supplements "as is". I have seen software that include this as a feature. We could include a downloadable PDF report in the output. This should be trivial once we have the `html` document - `pandoc` can do it for us as mentiond in this r-bloggers [page](http://www.r-bloggers.com/converting-a-markdown-file-to-pdf-using-pandoc).   

## Releases ##

- `asifzubair/atacseeker:v5`
- `asifzubair/atacseeker:v6`
	- move output files instead of copy
	- added mtDNA assembly code
	- complexity code: OFF
	- fixed axis labels for barplots: used `cex.names = 0.8`
	- LOLA code: OFF
- `asifzubair/atacseeker:v7`
	- prevent downloading of all files by fixing report_builder
	- added chunk names
	- have a flag to do mtDNA analysis
	- fix `assembleMTgenome.py` to not sort and index
	- fix plot size for fragment size distribution
	- write UP and DOWN to separate HTMLs
- `asifzubair/atacseeker:v8`
	- `assembleMTgenome.sh` will clean up after itself
	- install Gviz
	- fix plots to adjust `cex.names`
	- add max.width for merge
	- visualization for top region
	- bam2bigwig converter
		- installed `bedtools` & [bedGraphToBigWig](http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64)
		- got chrom.sizes from [here](https://genome.ucsc.edu/goldenpath/help/hg19.chrom.sizes)
	- add best window information to diffAccRegions

## ATACseq Data Analysis ##

Notes from `ILMN` [interview](http://www.illumina.com/content/dam/illumina-marketing/documents/icommunity/greenleaf-stanford-interview-miseq-hiseq-cancer-immune-1070-2015-003.pdf?mkt_tok=3RkMMJWWfF9wsRokv6%2FBdu%2FhmjTEU5z16eglWK%2B0hIkz2EFye%2BLIHETpodcMTcdgM7DYDBceEJhqyQJxPr3DLNANwtBlRhjgDw%3D%3D) of Greenleaf:

> You can determine whether the reads that you’re getting are from __nucleosomal or nucleosome-free regions__ based on the predicted length. 

> There’s no reason not to do paired-end sequencing when generating ATAC-seq libraries. Both ends tell you the insertion point of a transposase — so it really doubles your data. It also tells you the length of each of the fragments generated from your library. When you align each of the ends to the genome, you can determine how long that fragment is. That length can be used to identify, for example, nucleosome positions. We recently reported a method for using the length distribution of the ATAC-seq fragments to call nucleosome positions at high resolution.

### data and experiment ###
