function launchSpec(dataProvider){

   var retval = {
 
            commandLine: ["Rscript","-e 'rmarkdown::render('/atacseeker/scripts/atacseeker.Rmd')'"],
            containerImageId: "asifzubair/atacseeker:v1",
            options:["bsfs.enabled=true"]
            
           };
           
   return retval;  
     
}
