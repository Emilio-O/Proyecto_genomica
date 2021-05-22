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
AdjMatC1<-sign(cor(datC1,method="spearman"))*(cor(datC1,method="spearman"))^2
AdjMatC2<-sign(cor(datC2,method="spearman"))*(cor(datC2,method="spearman"))^2
diag(AdjMatC1)<-0
diag(AdjMatC2)<-0
collectGarbage()

























data<-as.matrix(read.csv(file="GDS2901.soft",skip=166,row.names=1,sep="\t",header=T))
data<-data[-15924,]
rawData<-matrix(as.numeric(data[,-1]),nrow=15923)
dimnames(rawData)<-dimnames(data[,-1])
anno<-as.matrix(data[-2475,1]) 
normData<-normalize.quantiles(log2(rawData))
dimnames(normData)<-dimnames(rawData)

uwu2 <- apply(normData, 2, var)
uwu2
