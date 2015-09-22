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
    parameter = data[['Properties']][['Items']][[index]][['Content']]
    params = c(params, parameter)    
  }
  
  if (data[['Properties']][['Items']][[index]]['Name'] == 'Input.control_label'){
    parameter = data[['Properties']][['Items']][[index]][['Content']]
    params = c(params, parameter)    
  }
  
  if (data[['Properties']][['Items']][[index]]['Name'] == 'Input.compare_label'){
    parameter = data[['Properties']][['Items']][[index]][['Content']]
    params = c(params, parameter)    
  }
  
  if (data[['Properties']][['Items']][[index]]['Name'] == 'Input.Projects'){
    parameter = data[['Properties']][['Items']][[index]][['Items']][[1]]['Id']
    params = c(params, parameter)    
  }
  
}

app_session_name = params[1]
control = params[2]
comparison = params[3]
project_id = params[4]


