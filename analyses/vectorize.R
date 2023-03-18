# ----------------------------------------------------------------------------
# [Title]: vectorize function.
# [Date]:  12 Mayo 2008    [Loc]: Sevilla
# Pedro Jordano.
# ----------------------------------------------------------------------------
## First version DATE. Revised DATE
# ----------------------------------------------------------------------------
# STACK: Turn table (C) into (D):
# C
#    a  b  c  d
# A  3  2  .  .
# B  .  .  1  1

# D
# A  a  3
# A  b  2
# A  c  .
# A  d  .
# B  a  .
# B  b  .
# B  c  1
# B  d  1
# The input matrix should be read with row.names=1 option 
# and a blank entry at the upper left corner.
#-------------------------------------------------------------
vectorize <- function(mat)
{
  mat <- t(mat)
  cbind(expand.grid(dimnames(mat))[2:1], as.vector(mat))
}

# Here vectorize for each population (JV, SAUC4, SAUC and ZAP): 
# a) genetic distance matrix
# b) spatial distance matrix
# c) kinship matrix
# d) MEMM pollen network
    # in the matrix, columns: mothers & rows: fathers
    # Previously binarizing the matrix:
quant2bin<-function(matr)
{
  cn=colnames(matr)
  rn=rownames(matr)
  ij<-which(matr>=0.05,arr.ind=T) # Adjust here the threshold value for mating.
  #  ij<-which(matr!=0,arr.ind=T)
  matrb<-matrix(0,dim(matr)[1],dim(matr)[2])
  matrb[ij]<-1
  colnames(matrb)=cn
  rownames(matrb)=rn
  return(matrb)
}
setwd("~/Desktop/PhDRocio/MS/MS_Fa_mating_patterns/Datos/dist_matrix")
setwd("~/Desktop/PhDRocio/MS/MS_Fa_mating_patterns/Datos")

## JV ##
# a)
mat <- read.table("JV_GenDist.txt", header=TRUE, row.names=1, sep="\t", dec=".", na.strings=".")
JV <- vectorize(mat)
write.table(as.matrix(JV),row.names=FALSE,col.names=FALSE,dec=".", 
            file="~/Desktop/PhDRocio/MS/MS_Fa_mating_patterns/Datos/dist_matrix/vectorize/JV_GenDist_v.txt")
# b)
mat <- read.table("JV_SpDist.txt", header=TRUE, row.names=1, sep="\t", dec=".", na.strings=".")
JV <- vectorize(mat)
write.table(as.matrix(JV),row.names=FALSE,col.names=FALSE,dec=".", 
            file="~/Desktop/PhDRocio/MS/MS_Fa_mating_patterns/Datos/dist_matrix/vectorize/JV_SpDist_v.txt")
#c)
mat <- read.table("JV_kinship.txt", header=TRUE, row.names=1, sep="\t", dec=".", na.strings=".")
JV <- vectorize(mat)
write.table(as.matrix(JV),row.names=FALSE,col.names=FALSE,dec=".", 
            file="~/Desktop/PhDRocio/MS/MS_Fa_mating_patterns/Datos/dist_matrix/vectorize/JV_kin_v.txt")
# d.1) row: father / column: mother
matingJVq<-as.matrix(read.table("JV_pollenNetwork.txt",header=T,sep="\t",dec=".",na.strings="NA", row.names=1))
matingJVfm<-quant2bin(matingJVq)
JVfm <- vectorize(matingJVfm)
write.table(as.matrix(JVfm),row.names=FALSE,col.names=FALSE,dec=".", 
            file="~/Desktop/PhDRocio/MS/MS_Fa_mating_patterns/Datos/dist_matrix/vectorize/JV_binpolnet_fm_v.txt")
# d.2) row: mother / column: father
matingJVmf <- t(matingJVfm)
JVmf <- vectorize(matingJVmf)
write.table(as.matrix(JVmf),row.names=FALSE,col.names=FALSE,dec=".", 
            file="~/Desktop/PhDRocio/MS/MS_Fa_mating_patterns/Datos/dist_matrix/vectorize/JV_binpolnet_mf_v.txt")


## SAUC4 ##
# a)
mat <- read.table("SAUC4_GenDist.txt", header=TRUE, row.names=1, sep="\t", dec=".", na.strings=".")
SAUC4 <- vectorize(mat)
write.table(as.matrix(SAUC4),row.names=FALSE,col.names=FALSE,dec=".", 
            file="~/Desktop/PhDRocio/MS/MS_Fa_mating_patterns/Datos/dist_matrix/vectorize/SAUC4_GenDist_v.txt")
# b)
mat <- read.table("SAUC4_SpDist.txt", header=TRUE, row.names=1, sep="\t", dec=".", na.strings=".")
SAUC4 <- vectorize(mat)
write.table(as.matrix(SAUC4),row.names=FALSE,col.names=FALSE,dec=".", 
            file="~/Desktop/PhDRocio/MS/MS_Fa_mating_patterns/Datos/dist_matrix/vectorize/SAUC4_SpDist_v.txt")
# c)
mat <- read.table("SAUC4_kinship.txt", header=TRUE, row.names=1, sep="\t", dec=".", na.strings=".")
SAUC4 <- vectorize(mat)
write.table(as.matrix(SAUC4),row.names=FALSE,col.names=FALSE,dec=".", 
            file="~/Desktop/PhDRocio/MS/MS_Fa_mating_patterns/Datos/dist_matrix/vectorize/SAUC4_kin_v.txt")
# d.1) row: father / column: mother
matingSauc4q<-as.matrix(read.table("Sauc4_pollenNetwork.txt",header=T,sep="\t",dec=".",na.strings="NA",row.names=1))
matingSauc4fm<-quant2bin(matingSauc4qfm) 
SAUC4fm <- vectorize(matingSauc4)
write.table(as.matrix(SAUC4fm),row.names=FALSE,col.names=FALSE,dec=".", 
            file="~/Desktop/PhDRocio/MS/MS_Fa_mating_patterns/Datos/dist_matrix/vectorize/SAUC4_binpolnet_fm_v.txt")
# d.2) row: mother / column: father
matingSauc4mf <- t(matingSauc4fm)
SAUC4mf <- vectorize(matingSauc4mf)
write.table(as.matrix(SAUC4mf),row.names=FALSE,col.names=FALSE,dec=".", 
            file="~/Desktop/PhDRocio/MS/MS_Fa_mating_patterns/Datos/dist_matrix/vectorize/SAUC4_binpolnet_mf_v.txt")


## SAUC ##
# a)
mat <- read.table("SAUC_GenDist.txt", header=TRUE, row.names=1, sep="\t", dec=".", na.strings=".")
SAUC <- vectorize(mat)
write.table(as.matrix(SAUC),row.names=FALSE,col.names=FALSE,dec=".", 
            file="~/Desktop/PhDRocio/MS/MS_Fa_mating_patterns/Datos/dist_matrix/vectorize/SAUC_GenDist_v.txt")
# b)
mat <- read.table("SAUC_SpDist.txt", header=TRUE, row.names=1, sep="\t", dec=".", na.strings=".")
SAUC <- vectorize(mat)
write.table(as.matrix(SAUC),row.names=FALSE,col.names=FALSE,dec=".", 
            file="~/Desktop/PhDRocio/MS/MS_Fa_mating_patterns/Datos/dist_matrix/vectorize/SAUC_SpDist_v.txt")
# c)
mat <- read.table("SAUC_kinship.txt", header=TRUE, row.names=1, sep="\t", dec=".", na.strings=".")
SAUC <- vectorize(mat)
write.table(as.matrix(SAUC),row.names=FALSE,col.names=FALSE,dec=".", 
            file="~/Desktop/PhDRocio/MS/MS_Fa_mating_patterns/Datos/dist_matrix/vectorize/SAUC_kin_v.txt")
# d.1) row: father / column: mother
matingSaucq<-as.matrix(read.table("Sauc_pollenNetwork.txt",header=T,sep="\t",dec=".",na.strings="NA",row.names=1))
matingSaucfm<-quant2bin(matingSaucq) 
SAUCfm <- vectorize(matingSaucfm)
write.table(as.matrix(SAUCfm),row.names=FALSE,col.names=FALSE,dec=".", 
            file="~/Desktop/PhDRocio/MS/MS_Fa_mating_patterns/Datos/dist_matrix/vectorize/SAUC_binpolnet_fm_v.txt")
# d.2) row: mother / column: father
matingSaucmf <- t(matingSaucfm) 
SAUCmf <- vectorize(matingSaucmf)
write.table(as.matrix(SAUCmf),row.names=FALSE,col.names=FALSE,dec=".", 
            file="~/Desktop/PhDRocio/MS/MS_Fa_mating_patterns/Datos/dist_matrix/vectorize/SAUC_binpolnet_mf_v.txt")


## ZP ##
# a)
mat <- read.table("ZAP_GenDist.txt", header=TRUE, row.names=1, sep="\t", dec=".", na.strings=".")
ZAP <- vectorize(mat)
write.table(as.matrix(ZAP),row.names=FALSE,col.names=FALSE,dec=".", 
            file="~/Desktop/PhDRocio/MS/MS_Fa_mating_patterns/Datos/dist_matrix/vectorize/ZAP_GenDist_v.txt")
# b)
mat <- read.table("ZAP_SpDist.txt", header=TRUE, row.names=1, sep="\t", dec=".", na.strings=".")
ZAP <- vectorize(mat)
write.table(as.matrix(ZAP),row.names=FALSE,col.names=FALSE,dec=".", 
            file="~/Desktop/PhDRocio/MS/MS_Fa_mating_patterns/Datos/dist_matrix/vectorize/ZAP_SpDist_v.txt")
# c)
mat <- read.table("ZAP_kinship.txt", header=TRUE, row.names=1, sep="\t", dec=".", na.strings=".")
ZAP <- vectorize(mat)
write.table(as.matrix(ZAP),row.names=FALSE,col.names=FALSE,dec=".", 
            file="~/Desktop/PhDRocio/MS/MS_Fa_mating_patterns/Datos/dist_matrix/vectorize/ZAP_kin_v.txt")
# d.1) row: father / column: mother
matingZAPq<-as.matrix(read.table("ZAP_pollenNetwork.txt",header=T,sep="\t",dec=".",na.strings="NA",row.names=1))
matingZAPfm<-quant2bin(matingZAPq)
ZAPfm <- vectorize(matingZAPfm)
write.table(as.matrix(ZAPfm),row.names=FALSE,col.names=FALSE,dec=".", 
            file="~/Desktop/PhDRocio/MS/MS_Fa_mating_patterns/Datos/dist_matrix/vectorize/ZAP_binpolnet_fm_v.txt")
# d.2) row: mother / column: father
matingZAPmf <- t(matingZAPfm)
ZAPmf <- vectorize(matingZAPmf)
write.table(as.matrix(ZAPmf),row.names=FALSE,col.names=FALSE,dec=".", 
            file="~/Desktop/PhDRocio/MS/MS_Fa_mating_patterns/Datos/dist_matrix/vectorize/ZAP_binpolnet_mf_v.txt")



