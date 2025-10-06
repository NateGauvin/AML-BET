library(mongolite)

# Connecting to MongoDB for GitHub actions
connect_mongo <- function(collection_name, user = "root", pass = "password", host = "aml-bet_db:27017") {
 uri <- sprintf("mongodb://%s:%s@%s/", user, pass, host)
 return (mongo(url = uri, db = "aml-bet", collection = collection_name))

}

# Overloaded method for RStudio use
# connect_mongo <- function(collection_name) {
#   return(mongo(collection = collection_name, db = "aml-bet", url = "mongodb://localhost:27017"))
# }

# Overloaded method for docker use
# connect_mongo <- function(collection_name) {
#   return(mongo(collection = collection_name, db = "aml-bet", url = "mongodb://aml-bet_db:27017"))
# }