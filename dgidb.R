#Systematic interrogation of DGIdb
	#review DGIdb terms of service before using 

library(jsonlite)

#function returns a list with known drug-gene interactions
gene_drug_interactions<-function(gen){
  interrogador<-paste0("http://dgidb.genome.wustl.edu/api/v1/interactions.json?genes=", gen)
  query_dgidb<-fromJSON(interrogador) #USE jsonlite version
  if(length(query_dgidb$matchedTerms)==0){
    return(NULL)
  }
  name_drugs<-character()
  drug_interaction<-character()
  for(i in seq_along(query_dgidb$matchedTerms$interactions[[1]]$drugName)){
    name_drugs<-c(name_drugs, query_dgidb$matchedTerms$interactions[[1]]$drugName[i])
    drug_interaction<-c(drug_interaction, query_dgidb$matchedTerms$interactions[[1]]$interactionType[i])
  }
  names(drug_interaction)<-name_drugs
  return(drug_interaction)
}

#For a set of genes 
#For instance, those identified to be deregulated in BrCan
genes # character vector with genes 
interacciones <- apply(X = genes, FUN = gene_drug_interactions)
