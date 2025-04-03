library(mongolite)

# Connecting to MongoDB
connect_mongo <- function(collection_name, user = "root", pass = "password", host = "localhost:27017") {
  uri <- sprintf("mongodb://%s:%s@%s/", user, pass, host)
  return (mongo(url = uri, db = "aml-bet", collection = collection_name))
  
}

# Overloaded method for RStudio use
#connect_mongo <- function(collection_name) {
#  return(mongo(collection = collection_name))
#}