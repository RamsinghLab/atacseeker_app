####
## This file parses the AppSession JSON. 
####

require("rjson")

data <- fromJSON(file="../../data/input/AppSession.json")
numberOfPropertyItems = length(data[['Properties']][['Items']])

##
# Let's fetch the property items
##

params = c()

controlID = c()
controlHref = c()
controlName = c()
controlDir = c()

compareID = c()
compareHref = c()
compareName = c()
compareDir = c()

##
# Collect all the AppSession properties
##

compute_library_complexity = FALSE
do_mtdNA_analysis = FALSE
do_motif_analysis = FALSE
generate_BigWigs = FALSE

for (index in seq(numberOfPropertyItems)){
  
    if (data[['Properties']][['Items']][[index]]['Name'] == 'Input.app-session-name'){
        app_session_name = data[['Properties']][['Items']][[index]][['Content']]
        params = c(params, app_session_name)    
    }

    if (data[['Properties']][['Items']][[index]][['Name']] == 'Input.control_label'){
        control = data[['Properties']][['Items']][[index]][['Content']]
        params = c(params, control)    
    }

    if (data[['Properties']][['Items']][[index]][['Name']] == 'Input.compare_label'){
        comparison = data[['Properties']][['Items']][[index]][['Content']]
        params = c(params, comparison)    
    }

    if (data[['Properties']][['Items']][[index]][['Name']] == 'Input.Projects'){
        project_id = data[['Properties']][['Items']][[index]][['Items']][[1]][['Id']]
        params = c(params, project_id)    
    }
 
     if (data[['Properties']][['Items']][[index]][['Name']] == 'Input.library-complexity'){
        compute_library_complexity = TRUE
    }
    
    if (data[['Properties']][['Items']][[index]][['Name']] == 'Input.mtdna-analysis'){
        do_mtdNA_analysis = TRUE
    }

    if (data[['Properties']][['Items']][[index]][['Name']] == 'Input.motif-analysis'){
        do_motif_analysis = TRUE
    }
    
    if (data[['Properties']][['Items']][[index]][['Name']] == 'Input.bigwigs'){
        generate_BigWigs = TRUE
    }
}
params = c(params, do_mtdNA_analysis) 

##
# Collect all sample names and IDs
##

for (index in seq(numberOfPropertyItems)){

    ## Control samples  
    if (data[['Properties']][['Items']][[index]][['Name']] == 'Input.control-app-result-id'){
        for (sample in seq(length(data[['Properties']][['Items']][[index]][['Items']]))){
            controlID = c(controlID, data[['Properties']][['Items']][[index]][['Items']][[sample]][['Id']])
            controlHref = c(controlHref, data[['Properties']][['Items']][[index]][['Items']][[sample]][['Href']])
            controlName = c(controlName, data[['Properties']][['Items']][[index]][['Items']][[sample]][['Name']])
        }
    }

    ## Compare samples
    if (data[['Properties']][['Items']][[index]][['Name']] == 'Input.compare-app-result-id'){
        for (sample in seq(length(data[['Properties']][['Items']][[index]][['Items']]))){
            compareID = c(compareID, data[['Properties']][['Items']][[index]][['Items']][[sample]][['Id']])
            compareHref = c(compareHref, data[['Properties']][['Items']][[index]][['Items']][[sample]][['Href']])
            compareName = c(compareName, data[['Properties']][['Items']][[index]][['Items']][[sample]][['Name']])
        }
    }
}

