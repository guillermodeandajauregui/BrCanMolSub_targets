
library(limma)

design<-read.table(file = "subtipos_design.txt") # file with samples and subtypes
###like this###
# sample ctrl lumA  lumB  Basal HER2
# EJE000  1     0     0     0     0
# EJE001  0     1     0     0     0
# EJE003  0     0     0     1     0
# EJE002  0     0     1     0     0
# EJE004  0     0     0     0     1

exp_matrix <- read.table(file = "EXPMATRIX_AFFYID.txt") #expression matrix

#contrast matrix defining comparisons to be made
cont.matrix<- makeContrasts('lumA - sanos', 'lumB - sanos', 'Basal - sanos', 'HER2 - sanos', levels = design)

#adjust data to a lineal model
fit <- lmFit(object = exp_matrix, design = design, method = robust)
#reorient lineal model by comparisons 
fit2<-contrasts.fit(fit = fit, contrasts = cont.matrix)
# Calculate Bayesian statistics
final_fit<-eBayes(fit = fit2)
# write out 
write.table(file = "contrastes.txt", x = topTable(final_fit, coef = 1, adjust = 'fdr', n = length(row.names(matriz_expresion))))
#do for each coefficient to have all contrasts
