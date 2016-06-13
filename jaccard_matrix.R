#Construction of a pathway crosstalk network
#where nodes are pathways 
#and edges represent crosstalk between pathways
	#weighted by Jaccard index between pathways 


#Functions

##jaccard: calculates jaccard index between two sets
jaccard<-function(a,b)
{
  x<-intersect(a,b)
  y<-union(a,b)
  nx<-length(x)
  ny<-length(y)
  J<-as.numeric(nx/ny)
  print(J)
  return(J)
}

#jaccard_matrix: calculate Jaccard index between two lists of sets
				#returns a Jaccard matrix 
				#with J between each set in list 
jaccard_matrix<-function(x,y){
  matriz<-matrix(nrow = length(x), ncol = length(y))  
  colnames(matriz)<-names(y)
  rownames(matriz)<-names(x)
  for (i in seq_along(x)){
    for (j in seq_along(y)){
      alfa<-x[[i]]
      beta<-y[[j]]
      matriz[i,j]<-jaccard(alfa,beta)
    }
  }
  print(matriz)
  return(matriz)  
}


##
##dir with pathways in .txt files
file_list<-list.files()
Pathways<-list()
for(i in file_list){
  x<-as.matrix(read.table(i, header=TRUE))
  assign(paste(colnames(x)), i) ##nombre de columna = nombre variable
  y<-(read.table(i, header=TRUE))
  Pathways<-c(Pathways, y)
  rm(x)
  rm(i)
}

jaccard_matrix(x = Pathways, y = Pathways)
#returns Jaccard matrix of Pathways v Pathways 
#Equivalent to an adjacency matrix for a pathway crosstalk network
	#If J = 0, no edge 
	#weighted edge otherwise

