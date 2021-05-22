# BiocManager::install("GenomicDataCommons")
setwd("C:/Users/52442/Desktop/Proyecto_genomica/")

# BiocManager::install("TCGAbiolinks")
library(WGCNA)
library(preprocessCore)
memory.size(max = T)
# ==============================================================================================================================
# Cargar base de datos y pre-procesamiento
# ==============================================================================================================================

data <- as.matrix(read.csv(file="GDS4299.soft",skip=90,row.names=1,sep="\t",header=T))
      # En la base de datos a partir de la linea 90 empiezan los nombres de las columnasx
      # Pero tambien tiene donde termina
head(data)


tail(data) # !dataset_table_end en el ultimo renglon
dim(data)  # Conociendo las dimensiones podemos eliminarlo
data <- data[-54716,]
tail(data) # Confirmamos



rawData <- matrix(as.numeric(data[,-1]), nrow = 54715) # Vamos a quedarnos con los datos crudos
                                                       # Solamente los valores dumericos
dimnames(rawData) <- dimnames(data[,-1]) # Volvemos a asignar los nombres, nos ahorramos col/rownames con esta funcion





normData<-normalize.quantiles(log2(rawData)) # Normalizamos los datos crudos
dimnames(normData) <- dimnames(rawData) # Asignamos nombres


ubic <- normData[apply(normData, 1, var) != 0, ] # Ubicar cuales son los renglones que tienen varianza de cero
# Comparamos dimensiones
dim(ubic)
dim(normData)
# Ningun renglon tiene varianza de cero

anno <- as.matrix(data[,1])




head(normData)
datC1 <- t(normData[ , c(8:47)])    # Muestras que corresponden a nonETP
datC2 <- t(normData[ , c(1:7, 48:52)]) # Muestras que corresponden a ETP




# ==============================================================================================================================
# Aplicar DiffCoEx
# ==============================================================================================================================
beta1 <- 6 #user defined parameter for soft thresholding

# AdjMatC1 <-sign(cor(datC1,method="spearman"))*(cor(datC1,method="spearman"))^2
# AdjMatC2 <-sign(cor(datC2,method="spearman"))*(cor(datC2,method="spearman"))^2

# diag(AdjMatC1)<-0
# diag(AdjMatC2)<-0

# save(AdjMatC1, file = "AdjMatC1.RData")
# save(AdjMatC2, file = "AdjMatC2.RData")

# load("AdjMatC1.RData")
# load("AdjMatC2.RData")

# collectGarbage()

# dissTOMC1C2 <- AdjMatC1-AdjMatC2

# rm(AdjMatC1)
# rm(AdjMatC2)

# dissTOMC1C2 <- TOMdist((abs(dissTOMC1C2)/2)^(beta1/2))

# dissTOMC1C2 <- TOMdist((abs(AdjMatC1-AdjMatC2)/2)^(beta1/2))

# save(dissTOMC1C2, file = "dissTOMC1C2.RData")

load("dissTOMC1C2.RData")
collectGarbage()

#Hierarchical clustering is performed using the Topological Overlap of the adjacency difference as input distance matrix
geneTreeC1C2 <- flashClust(as.dist(dissTOMC1C2), method = "average");

# Plot the resulting clustering tree (dendrogram)
png(file="hierarchicalTree.png",height=1000,width=1000)
plot(geneTreeC1C2, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",labels = FALSE, hang = 0.04);
dev.off()




### AQUI VAMOS
#We now extract modules from the hierarchical tree. This is done using cutreeDynamic. Please refer to WGCNA package documentation for details
dynamicModsHybridC1C2 <- cutreeDynamic(dendro = geneTreeC1C2, distM = dissTOMC1C2,method="hybrid",cutHeight=.996,deepSplit = T, pamRespectsDendro = FALSE,minClusterSize = 20);

#Every module is assigned a color. Note that GREY is reserved for genes which do not belong to any differentially coexpressed module
dynamicColorsHybridC1C2 <- labels2colors(dynamicModsHybridC1C2)


#the next step merges clusters which are close (see WGCNA package documentation)
mergedColorC1C2<-mergeCloseModules(rbind(datC1,datC2),dynamicColorsHybridC1C2,cutHeight=.2)$color
colorh1C1C2<-mergedColorC1C2


#reassign better colors
colorh1C1C2[which(colorh1C1C2 =="midnightblue")]<-"red"
colorh1C1C2[which(colorh1C1C2 =="lightgreen")]<-"yellow"
colorh1C1C2[which(colorh1C1C2 =="cyan")]<-"orange"
colorh1C1C2[which(colorh1C1C2 =="lightcyan")]<-"green"
# Plot the dendrogram and colors underneath
png(file="module_assignment.png",width=1000,height=1000)
plotDendroAndColors(geneTreeC1C2, colorh1C1C2, "Hybrid Tree Cut",dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05,main = "Gene dendrogram and module colors cells")
dev.off()

#We write each module to an individual file containing affymetrix probeset IDs
modulesC1C2Merged<-extractModules(colorh1C1C2,datC1,anno,dir="modules",file_prefix=paste("Output","Specific_module",sep=''),write=T)
write.table(colorh1C1C2,file="module_assignment.txt",row.names=F,col.names=F,quote=F)

#We plot to a file the comparative heatmap showing correlation changes in the modules
#The code for the function plotC1C2Heatmap and others can be found below under the Supporting Functions section

plotC1C2Heatmap(colorh1C1C2,AdjMatC1,AdjMatC2, datC1, datC2)
png(file="exprChange.png",height=500,width=500)
plotExprChange(datC1,datC2,colorh1C1C2)

























data<-as.matrix(read.csv(file="GDS2901.soft",skip=166,row.names=1,sep="\t",header=T))
data<-data[-15924,]
rawData<-matrix(as.numeric(data[,-1]),nrow=15923)
dimnames(rawData)<-dimnames(data[,-1])
anno<-as.matrix(data[-2475,1]) 
normData<-normalize.quantiles(log2(rawData))
dimnames(normData)<-dimnames(rawData)

uwu2 <- apply(normData, 2, var)
uwu2
