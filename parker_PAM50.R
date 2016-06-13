###
# PAM50 classification of Breast cancer samples. 
#Adapted from material available in https://genome.unc.edu/pubsup/breastGEO/PAM50.zip
	#by Parker et al. 
		#for use on breast cancer samples used by Tovar et al 2015, de Anda-Jáuregui et al 2015
		# to be run as Rscript
###
library(ctc)
library(heatmap.plus)
paramDir<- "/home/bioclassifier_R" # location of functions distributed by 
				   #Parker et al in https://genome.unc.edu/pubsup/breastGEO/PAM50.zip
inputDir<- "/home"  # Exp Matrix localization
inputFile<- "EXPMATRIX.txt" # TAB SEPARATED exp matrix 
short<-"KDB" # output file names 
calibrationParameters<- NA 	# to adjust centroids based on the training set
hasClinical<-FALSE 	# dataset does not include clinical data
collapseMethod<-"median" # expmatrix previously collapsed.  
####
# run the assignment algorithm
####
source(paste(paramDir,"subtypePrediction_functions.R",sep="/"))
source(paste(paramDir,"subtypePrediction_distributed.R",sep="/"))
##Output: table with samples and classification
##	  PCA visualization
