# Created on June 13, 2023
# Author: Jung Hun Oh
# github: https://github.com/MSK-MOI/NWK


# Install packages
install.packages("igraph")
install.packages("Matrix")
install.packages("qgraph")

library("qgraph")
library("igraph")
library("Matrix")

# Read radiomic data
data_raw <- read.csv("lung_data.csv",header=T) 
id <- data_raw[,1]
data<-data_raw[,2:392]
m<-length(data) # number of features
variables<-colnames(data) # feature names

# Remove features that are highly correlated with the tumor volume
vol <- data[,1]
COR_VOl_THRESH=0.8 # Spearman correlation threshold
cor_vol=numeric(m)

for (i in 1:m){
  spear_cor = cor(vol, data[,i], method = "spearman")
  cor_vol[i]<-spear_cor
}
ind <- which(abs(cor_vol)  < COR_VOl_THRESH) 
ind<-c(1,ind) # 1: vol index

# Extract a subset of features 
data_red <- data[,ind]
m_red <- length(data_red)
variables_red <- variables[ind]

# Normalization of features with negative values for unbalanced OMT
data_red_nor <- data_red
ind1<-c()
for (i in 1:m_red){
  if (min(data_red[,i]) < 0){
    ind1<-c(ind1,i)
    x<-data_red[,i]
    nor_x<-(x-min(x))/(max(x)-min(x))
    data_red_nor[,i]<-nor_x
  }
}

# To ensure the matrix is positive definite 
cormatrix <- cor(data_red,method='pearson') 
cor_PD<-nearPD(cormatrix,corr=TRUE)
CorMat<-as.matrix(cor_PD$mat)

# EBIC LASSO (partial correlation-based graphical LASSO)
EBICgraph <- EBICglasso(CorMat, nrow(data_red), 0.5, threshold = TRUE)

# Draw graph
graph <- qgraph(EBICgraph, layout = "spring", details = TRUE)
g = as.igraph(graph, attributes=TRUE)
gg<-as.undirected(g)
com<-components(gg, mode = c("weak", "strong"))

# Save EBIC results and data
write.csv(EBICgraph, file = "EBICgraph.csv", row.names = variables_red)
write.csv(cbind(id, data_red_nor), file = "lung_data_red.csv", row.names = FALSE)
