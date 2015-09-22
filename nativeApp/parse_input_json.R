# This file parses the input json. 

require("jsonlite")


data <- fromJSON(file="AppSession.json")
numberOfPropertyItems = length(data[['Properties']][['Items']])


# Let's fetch the property items
params = c()

controlID = c()
controlHref = c()
controlName = c()
controlDir = c()

compareID = c()
compareHref = c()
compareName = c()
compareDir = c()

for (index in seq(numberOfPropertyItems)){
  
  if (data[['Properties']][['Items']][[index]]['Name'] == 'Input.app-session-name'){
    app_session_name = data[['Properties']][['Items']][[index]][['Content']]
    params = c(params, app_session_name)    
  }
  
  if (data[['Properties']][['Items']][[index]]['Name'] == 'Input.control_label'){
    control = data[['Properties']][['Items']][[index]][['Content']]
    params = c(params, control)    
  }
  
  if (data[['Properties']][['Items']][[index]]['Name'] == 'Input.compare_label'){
    comparison = data[['Properties']][['Items']][[index]][['Content']]
    params = c(params, comparison)    
  }
  
  if (data[['Properties']][['Items']][[index]]['Name'] == 'Input.Projects'){
    project_id = data[['Properties']][['Items']][[index]][['Items']][[1]][['Id']]
    params = c(params, project_id)    
  }
  
}




