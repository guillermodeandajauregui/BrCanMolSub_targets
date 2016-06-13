#GAGE pathway enrichment

library("gage")
matriz_expresion <- read.table(file = "EXPMATRIX.txt") #exp matrix
subtipos<-read.table(file = "subtipos_factores.txt") #table with sample names and subtypes (or control)
#to be used with a list of pathways 
#dir with pathways in .txt files
file_list<-list.files()
Pathways<-list()
for(i in file_list){
  x<-as.matrix(read.table(i, header=TRUE))
  assign(paste(colnames(x)), i) 
  y<-(read.table(i, header=TRUE))
  Pathways<-c(Pathways, y)
  rm(x)
  rm(i)
}

#groups
sanos = which(colnames(matriz_expresion)%in%subtipos[subtipos$subtipos=="SANOS",1]) #control group, "healthy"
lumA = which(colnames(matriz_expresion)%in%subtipos[subtipos$subtipos=="lumA",1])
lumB = which(colnames(matriz_expresion)%in%subtipos[subtipos$subtipos=="lumB",1])
Basal = which(colnames(matriz_expresion)%in%subtipos[subtipos$subtipos=="Basal",1])
HER2 = which(colnames(matriz_expresion)%in%subtipos[subtipos$subtipos=="HER2",1])

#GAGE function does the comparing
luma_v_sanos <-gage(
  data = matriz_expresion, #exp matrix
  gsets = Pathways, # pathways to be analyzed 
  ref = sanos,
  samp = lumA,
  compare = "unpaired" #groups are of unequal size, 
                       #sample-to-sample comparisons 
)

#retrieve only significant pathways 
luma_v_sanos_sig<-sigGeneSet(luma_v_sanos)
#repeat for each comparison

