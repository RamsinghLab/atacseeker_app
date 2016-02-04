# ATACseeker App #

## Mt-DNA Variant Calling ##

Because in ATACseq a lot of the reads are biased for the mitochondrial DNA, it seems like an obvious thing to use ATACseq reads to call variants on mt-DNA. There are a few caveats though. The mt-DNA is haploid, inherited completely from the mother. However, most callers like `gatk` and `samtools` assume a diploid genome. 

This [post](http://gatkforums.broadinstitute.org/gatk/discussion/1214/can-i-use-gatk-on-non-diploid-organisms) says that post-version 3.3 `gatk` should be able to handle non-diploid genome. But a later [post](http://gatkforums.broadinstitute.org/gatk/discussion/3345/question-about-ploidy-when-mtdna-variants-calling) categorically denied any work by gatk relating to mt-DNA. 

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

```
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

I checked that `bedtools` can be used for checking [coverage](https://github.com/arq5x/bedtools-protocols/blob/master/bedtools.md) and I tried to do the same for one ATACseq sample. I can't remember the exact command, but it was probably something like this.

```
# Calculate a histogram of coverage for each chromosome
    # as well as genome-wide.
    # Estimate: 30 minutes
    bedtools genomecov \
             -ibam NA19146.bam \
    > NA19146.coverage.hist.txt
```

I remember that when I ran the above command, per nucleotide coverage in chrM was quite high - above 100 bases. 

## ATACseq Analysis Pipeline ##

### preseqR ###

### csaw ###

### lola ###

## Docker ##

## Native App ##


**testing basespace app**

`ilmn` has a weird way of testing these apps. essentially, i downloaded the native app development [image](https://da1s119xsxmu0.cloudfront.net/sites/developer/native/nativeappsvm/BaseSpace%20Native%20Apps%20VM%20v14.ova) and then loaded it in virtual box. You need to login to the image via the terminal - `ssh basespace@localhost -p2222` . The password is `basespace`. 

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
  ordinary text without R code

K
  |                                                                       K
  |....................                                             |  31%
label: unnamed-chunk-6 (with options) 
List of 1
 $ echo: logi FALSE

.Quitting from lines 125-136 (atacseeker.Rmd) 
�Error in plot.window(xlim, ylim, log = log, ...) : 
  need finite 'ylim' values
Calls: <Anonymous> ... eval -> eval -> barplot -> barplot.default -> plot.window

Execution halted
cp: 1cannot stat '/atacseeker/scripts/atacseeker.html': No such file or directory
```

So, it seems that creating index is failing. This might be because I have to write everything to scratch - which is pretty agonizing. 
