#Classification of samples by the citbcsmt algorithm
	#developed by Marisa et al. 
	#for use on breast cancer samples used by Tovar et al 2015, de Anda-Jáuregui et al 2015
#R code 
library(citbcmst) 
expmat_affy <- read.table(file = "EXPMATRIX_AFFYID.txt")  # expression matrix 
														  # not collapsed,
														  # probeset-level 
anotation <- data.frame(id=rownames(exp.norm.bertheau07), stringsAsFactors=FALSE,
                                   row.names=rownames(exp.norm.bertheau07) ) #annotation provided by Marisa et al 

asignacion_citbcmst <- cit.assignBcmst(  data = expmat_affy, 
                       data.annot = anotation, 
                       data.colId="id", 
                       data.colMap="id", 
                       citbcmst.annot=NULL,
                       dist.method="dlda",
                  plot=TRUE)
#returns dataframe with samples assigned to a subtype
