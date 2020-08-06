#--------------------------------------------------------------------------------------
# TEAM AND ACKNOWLEDGMENTS
#--------------------------------------------------------------------------------------

# This code is based on code and activities given by professors Edgar and Jesus (ITESM) 
# for the class Análisis de Biología Computacional - 2020. Modified and adapted for academic  
# purposes by the following team of students:

# Armando De la Rosa 	     
# Daniel Zamacona	               
# Ricardo Arriaga                
# Maribel Alvarez Flores          

#--------------------------------------------------------------------------------------
# LOADING DATA
#--------------------------------------------------------------------------------------

load("Data.RData")

# This .RData file has the following Data:

# geoGSE39084_clin  - Patients information Data Frame
# geoGSE39084_exp   - Gene Sample Matrix

#--------------------------------------------------------------------------------------
# FIXING DATA
#--------------------------------------------------------------------------------------

# RENAMING OF THE MATRIX AND DATAFRAME
patients           = geoGSE39084_clin       
samples            = geoGSE39084_exp

# INDEXES
patients[,3]       = as.numeric(patients[,3])
indexY             = which(patients$Age < 50)
indexO             = which(patients$Age > 50)

# CLASSES COLUMNS
Young              = samples[,indexY] 
Old                = samples[,indexO] 

# EMPTY MATRIX
analysis           = matrix(data=NA, nrow(samples), ncol=4) 

# ROWNAMES AND COLNAMES
rownames(analysis) = c(rownames(samples))
colnames(analysis) = c("Young","Old","P - Value","Mean Difference")

#--------------------------------------------------------------------------------------
# Average per sample (Young, Old), P-Value and Mean Difference per gene (rows) MATRIX
#--------------------------------------------------------------------------------------

for(i in c(1:nrow(analysis))) 
{
  test        = t.test(Young[i,], Old[i,])
  analysis[i,] = c(test$estimate[1], test$estimate[2], test$p.value, 
                   (test$estimate[1] -test$estimate[2])) 
}

#--------------------------------------------------------------------------------------
# Analysis P-values                                               |||||||||| INCREASING
#--------------------------------------------------------------------------------------

# NOT COMPLETELY, but still Significant Genes TOP 10
order(analysis[,3])[1:10]
analysis[order(analysis[,3])[1:10],] 

# CLEANING DATA (deleting rows with low expression genes either Young or Old) 
#THIS SHOULD YIELD 100% significant genes.

indexP              = which(analysis[,1:2] < 6)
new_analysis        = analysis[-indexP,]
order_analysis      = new_analysis[order(new_analysis[,3]),] 
View(order_analysis) 

# P-VALUE CONSIDERING GENE EXPRESSION, SIGNIFICANT GENES
index_order_pvalue  = order(order_analysis[,3])
order_pvalue        = order_analysis[index_order_pvalue,]
View(order_pvalue) 

#--------------------------------------------------------------------------------------
# Analysis  MEAN DIFFERENCE of YOUNG - OLD                        |||||||||| DECREASING
#--------------------------------------------------------------------------------------

# MEAN DIFFERENCE CONSIDERING GENE EXPRESSION, SIGNIFICANT GENES
index_order_mean    = order(abs(order_analysis[,4]),decreasing=T)
order_mean          = order_analysis[index_order_mean,]
View(order_mean) 

#--------------------------------------------------------------------------------------
# COMPARISION TOP 20, MEAN vs PVALUES SIGNIFICANT GENES, considering GENE EXPRESSION
#--------------------------------------------------------------------------------------

print(rownames(order_mean)[1:20]) 
print(rownames(order_pvalue)[1:20]) 

#--------------------------------------------------------------------------------------
# Signaling pathways (MEAN DIFFERENCE SIGNIFICANT GENES AREN'T CONSIDERED)
#--------------------------------------------------------------------------------------

# THE FOLLOWING TOP 50 are used to calculate the enrichment for cellular signaling 
# pathways using the MSigDB database.

# www.gsea-msigdb.org/gsea/msigdb/index.jsp

#_________________________________IMPORTANT FOR FINDING SIGNALING PATHS____________________________________________

# TOP 50 MORE EXPRESSED GENES IN YOUNG PEOPLE
write.table(rownames(order_analysis[which(order_analysis[,4] > 0),])[1:50], sep="\t",quote=F, row.names=F, col.names=F)

# TOP 50 LESS EXPRESSED GENES IN YOUNG PEOPLE
write.table(rownames(order_analysis[which(order_analysis[,4] < 0),])[1:50], sep="\t",quote=F, row.names=F, col.names=F)

#--------------------------------------------------------------------------------------
# VISUAL COMAPARISION OF MEAN Vs MEAN, P-VALUE Vs P-Value and Mean Vs P-Value
#--------------------------------------------------------------------------------------

par(mfrow = c(1,3))

# MEAN (eje x) vs MEAN (eje y)
plot(index_order_mean[1:100],index_order_mean[1:100],xlab="Mean",ylab="Mean")
abline(a=1,b=1)

# P-Value (eje x) vs P-Value (eje y)
plot(index_order_pvalue[1:100],index_order_pvalue[1:100],xlab="P-Value",ylab="P-Value")
abline(a=1,b=1)

# MEAN (eje x) VS P-Value (eje y)
plot(index_order_mean[1:100],index_order_pvalue[1:100],xlab="Mean",ylab="P-value") 
abline(a=1,b=2)



