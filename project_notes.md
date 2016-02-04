# ATACseeker App #

## Mt-DNA Variant Calling ##

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
