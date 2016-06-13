# Microarray data processing 
# libraries used

library(affy) 
library(frma)
#library(rma) if rma is prefered 

# Gather CEL files in a dir
Data = ReadAffy()

# Robust Multi-Array summarization/normalization
#Background correction, signal normalization - quantiles 
frmaData <- frma(Data, summarize="robust_weighted_average") 

#or rmaData <- rma(Data, summarize="robust_weighted_average") 
edata<-exprs(frmaData) # expression matrix
                       # each column is a sample (microarray)
                       # each row a transcript (probeset-level)
                       # each cell a normalized signal

################# BATCH EFFECT correction #########################
# libraries used
library(sva)

N=length(Data@phenoData@data$sample)
GSMnames <- colnames(edata)
GSMnames <- data.frame(lapply(GSMnames, function(v) {
  if (is.character(v)) return(toupper(v))
  else return(v)
}))
GSMnames <- as.character(GSMnames)
colnames(edata) <- GSMnames
batch <- c((rep(0,N))) # lotes 
# Cases ## data curated by Dr. Hugo Tovar
GSE1456<-read.table("/lists_of_names/all_GSE1456.txt", colClasses = "character") # This .txt files contain the list of GSM ID for any GSE
GSE1561<-read.table("/lists_of_names/all_GSE1561.txt", colClasses = "character") # This .txt files contain the list of GSM ID for any GSE
GSE2603<-read.table("/lists_of_names/all_GSE2603.txt", colClasses = "character") # This .txt files contain the list of GSM ID for any GSE
GSE2990<-read.table("/lists_of_names/all_GSE2990.txt", colClasses = "character") # This .txt files contain the list of GSM ID for any GSE
GSE3494<-read.table("/lists_of_names/all_GSE3494.txt", colClasses = "character") # This .txt files contain the list of GSM ID for any GSE
GSE4922<-read.table("/lists_of_names/all_GSE4922.txt", colClasses = "character") # This .txt files contain the list of GSM ID for any GSE
GSE7390<-read.table("/lists_of_names/all_GSE7390.txt", colClasses = "character") # This .txt files contain the list of GSM ID for any GSE
GSE1456 <- GSE1456[[1]]
GSE1561 <- GSE1561[[1]]
GSE2603 <- GSE2603[[1]]
GSE2990 <- GSE2990[[1]]
GSE3494 <- GSE3494[[1]]
GSE4922 <- GSE4922[[1]]
GSE7390 <- GSE7390[[1]]
n1456 <- which(GSMnames %in% GSE1456)
n1561 <- which(GSMnames %in% GSE1561)
n2603 <- which(GSMnames %in% GSE2603)
n2990 <- which(GSMnames %in% GSE2990)
n3494 <- which(GSMnames %in% GSE3494)
n4922 <- which(GSMnames %in% GSE4922)
n7390 <- which(GSMnames %in% GSE7390)
batch[n1456] = 1
batch[n1561] = 2
batch[n2603] = 3
batch[n2990] = 4
batch[n3494] = 5
batch[n4922] = 6
batch[n7390] = 7
case <- c(n1456, n1561, n2603, n2990, n3494, n4922, n7390)
# Controls ##data curated by Dr. Hugo Tovar
GSE15852 <-read.table("/lists_of_names/all_GSE15852.txt", colClasses = "character") # This .txt files contain the list of GSM ID for any GSE
GSE6883 <-read.table("/lists_of_names/all_GSE6883.txt", colClasses = "character") # This .txt files contain the list of GSM ID for any GSE
GSE9574 <-read.table("/lists_of_names/all_GSE9574.txt", colClasses = "character") # This .txt files contain the list of GSM ID for any GSE
GSE15852 <- GSE15852[[1]]
GSE6883 <- GSE6883 [[1]]
GSE9574 <- GSE9574 [[1]]
n15852 <- which(GSMnames %in% GSE15852)
n6883  <- which(GSMnames %in% GSE6883)
n9574  <- which(GSMnames %in% GSE9574)
batch[n15852] = 8
batch[n6883] = 9
batch[n9574] = 10
control <- c(n15852, n6883, n9574)
# matrices for control and cases, separated
caseExp <- edata[,-control]
controlExp <- edata[,-case]
# 
caseBatch <- c(batch[-control])
controlBatch <- c(batch[-case])

# Batch effect correction
case_combat = ComBat(dat=caseExp, batch=caseBatch, mod=NULL)
control_combat = ComBat(dat=controlExp, batch=controlBatch, mod=NULL)
combat_2ways <- matrix(rep(0,19609040), ncol=880)
colnames(combat_2ways) <- colnames(edata)
rownames(combat_2ways) <- rownames(edata)
x<- combat_2ways
write.table(x, file = "expmatrix_affyid.txt", sep = "\t")

############gene symbol annotation ######################
library("annotate")
library("hgu133a.db")
x<-as.data.frame(x)
genesymbols<-data.frame()
genesymbols<-(getSYMBOL(as.character(rownames(x)), "hgu133a.db"))
sym<-as.data.frame(genesymbols)
rownames(x)<-rownames(sym)
x$GS<-genesymbols
x$GS<-ifelse(is.na(x$GS), as.vector(rownames(x)), as.vector(x$GS))
x.med <- aggregate(. ~ GS, data = x, median) #Agregar valores de Genesymbol por mediana#
write.csv(x, file = "EXPMATRIX.txt")
#output EXPMATRIX
##############################################################################
