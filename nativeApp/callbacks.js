function launchSpec(dataProvider){

   var retval = {
 
            commandLine: ["/bin/bash","/atacseeker/scripts/atacseeker.sh"],
            containerImageId: "asifzubair/atacseeker:v2",
            options:["bsfs.enabled=true"]
            
           };
           
   return retval;  
     
}
