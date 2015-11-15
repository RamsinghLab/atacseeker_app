function launchSpec(dataProvider){

   var retval = {
 
            commandLine: ["/bin/bash","/atacseeker/scripts/atacseeker.sh"],
            containerImageId: "asifzubair/atacseeker:v5",
            options:["bsfs.enabled=true"]
            
           };
           
   return retval;  
     
}
