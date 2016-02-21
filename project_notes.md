# ATACseeker App #

## Mt-DNA Variant Calling ##

Because in ATACseq a lot of the reads are biased for the mitochondrial DNA, it seems like an obvious thing to use ATACseq reads to call variants on mt-DNA. There are a few caveats though. The mt-DNA is haploid, inherited completely from the mother. However, most callers like `gatk` and `samtools` assume a diploid genome. 

This [post](http://gatkforums.broadinstitute.org/gatk/discussion/1214/can-i-use-gatk-on-non-diploid-organisms) says that post-version 3.3 `gatk` should be able to handle non-diploid genome. But a later [post](http://gatkforums.broadinstitute.org/gatk/discussion/3345/question-about-ploidy-when-mtdna-variants-calling) categorically denied any work by gatk relating to mt-DNA. Also, this later [post](http://gatkforums.broadinstitute.org/gatk/discussion/3651/gatk-for-mitochondrial-dna-mutations) claims that the user might have to play around with the ploidy settings of `gatk` but does not offer any conclusive answers. 

`samtools mpileup` explicitly says that it is for diploid genomes - as seen here:

```
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

```
samtools mpileup -d 200 -D -B -f ../../tmp/reference/hg19.fa -r chrM -b bamlist.txt -u | \
bcftools view -  > variants.vcf
```

Still, no cigar ! 

Based on a recommendation by [GGD](http://www.gettinggeneticsdone.com/2014/03/visualize-coverage-exome-targeted-ngs-bedtools.html), I checked that `bedtools` can be used for checking [coverage](https://github.com/arq5x/bedtools-protocols/blob/master/bedtools.md) and I tried to do the same for one ATACseq sample. I can't remember the exact command, but it was probably something like this.

```
# Calculate a histogram of coverage for each chromosome
    # as well as genome-wide.
    # Estimate: 30 minutes
    bedtools genomecov \
             -ibam NA19146.bam \
    > NA19146.coverage.hist.txt
```

I remember that when I ran the above command, per nucleotide coverage in chrM was quite high - above 100 bases. 

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

#### qc, read-counting ####

[aaron](https://github.com/LTLA) really did a good job putting together `csaw`. the documentation is great and the options carefully chosen.  

one of the first things about the atacseq data is that it has a lot of reads biased to the mt-DNA. thus, we restrict analysis to the autosomes and sex-chromosomes. also, we take out blacklisted regions as defined [here](https://sites.google.com/site/anshulkundaje/projects/blacklists).  

programmatically, this can be implemented in `csaw` as `readParam(pe = "none", restrict = chroms, discard = blacklist, minq = 10, dedup = TRUE)` . we also, filter for bad quality reads and remove duplicates.  

aaron claims that removing dupes is not as straightforward as we think. since, it also caps the number of reads at each position. This can lead to loss of DB detection power in high abundance regions. Spurious differences may also be introduced when the same upper bound is applied to libraries of varying size. Thus, duplicate removal is not recommended for routine DB analyses. Of course, it may be unavoidable in some cases, e.g., involving libraries generated from low quantities of DNA.  

the first QC metrics we compute is the plot of `TSSvsNonTSS` counts. if i remember correctly, tim T said that this should be above 10 for good data. 

`csaw`'s `windowCounts` is the main function which does read counting over chromosomes. it does some read extension for forward and reverse reads to account for fragment length. for the atacseq data, i don't want to do any read extension and i think i need to set `ext` option to `NA` (it defaults to a 100). another thing i want to add to the pipeline output is the options used for `windowCounts`, something like `metadata(data)` should do this for me. this is cool stuff - read extension was bothering me for some time now.  

in light of the above, i think we need to change our `windowCounts` call from:
```
windowCounts(bam.files, ext=0, shift=4, bin=TRUE, param=discard.se.param)
```
to:
```
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

as we can see that `ext` is set to `1L` and also `filter` is set to `0/1`. So I guess, even this should work:

```
windowCounts(bam.files, shift=4, bin=TRUE, param=discard.se.param)
```

i got curious about `ext` and dug up the funciton to which `ext` is passed to get the `final.ext`:
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

i think i will stick with setting `ext=NA` but less sure now about how to interpret it.  


### LOLA ###



## Docker ##




## Native App ##


**testing basespace app**

`ilmn` has an  interesting way of testing these apps. essentially, i downloaded the native app development [image](https://da1s119xsxmu0.cloudfront.net/sites/developer/native/nativeappsvm/BaseSpace%20Native%20Apps%20VM%20v14.ova) and then loaded it in virtual box. You need to login to the image via the terminal - `ssh basespace@localhost -p2222` . The password is `basespace`. 

This image contains all software to test native apps - docker, spacedock etc. Also the versions should be up-to-date with the ones `ilmn` has. 

Associated with every app, is a `spacedock` command which in my case is - `sudo spacedock -a d31c319955504b5ba162a51a0026d693 -m https://mission.basespace.illumina.com` . The idea is to run this command locally inside the VM and then trigger a test by doing a `Send To Local Agent`.

My first attempt wasn't so successful. I checked the `spacedock` log and saw this relevant part. 

```
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

```
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

```
Quitting from lines 138-176 (atacseeker.Rmd) 
Error in mcmapply(getQC, name = control, pe.bam = control.files, SIMPLIFY = F) : 
  'names' attribute [2] must be the same length as the vector [1]
Calls: <Anonymous> ... withCallingHandlers -> withVisible -> eval -> eval -> mcmapply
```
 
... which is weird, because I am not seeing this when I run the app locally. My hunch here is that this could be a memory issue with `mcmapply`. 

So, my hunch was correct! Since these jobs are sent to my computer to be run locally, I could increase the memory available to the VM and that should facilitate job completion ... and it did! 

**report builder**

`ilmn` uses report builder which relies on liquid/javascript etc. but quite frankly fails to impress. I thought of publishing my report as a `html` document and then using iframes for serving this to the user. 

```javascript
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

## ATACseq Data Analysis ##

###data and experiment###