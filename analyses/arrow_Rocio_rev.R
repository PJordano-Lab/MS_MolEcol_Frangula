# setwd("~/Desktop/PhDRocio/MS/MS_Fa_mating_patterns/Datos")
library(here)
setwd("~/Documents/Working/MS_Tesis_Rocio/Data")

library(network);library(sna);library(bipartite);library(igraph);library(diagram);
library(adegenet);library(spdep);library(prabclus); library(ade4); library(mpmcorrelogram)

# I can use JV_all_net.txt, it is weigthed matrix (proportion MEMM * harvest fruits) or 
# JV_pollenNetwork.txt (expected proportion by MEMM; values between 0-1)

quant2bin<-function(matr)
{
  cn=colnames(matr)
  rn=rownames(matr)
  ij<-which(matr>=0.05,arr.ind=T) # Adjust here the threshold value for mating.
  # MEMM: green <0.05; 0.05 < yellow< 0.1; 0.1 < brown < 0.5; 0.5 < white < 1
  #  ij<-which(matr!=0,arr.ind=T)
  matrb<-matrix(0,dim(matr)[1],dim(matr)[2])
  matrb[ij]<-1
  colnames(matrb)=cn
  rownames(matrb)=rn
  return(matrb)
}
quant2binPJ<-function(matr)
  # Version modified by PJ to get quant matrix trimming threshold values.
  # This sets aij to 0 whenever the actual P value for a mating event
  # is P<=threshold_value. Then sets the other quantitative P's' to the same
  # value as the oriignal matrix.
{
  cn=colnames(matr)
  rn=rownames(matr)
  ij<-which(matr>=0.50,arr.ind=T) # Adjust here the threshold value for                    # mating.
  #  ij<-which(matr!=0,arr.ind=T)
  matrb<-matrix(0,dim(matr)[1],dim(matr)[2])
  matrb[ij]<-matr[ij]
  #  matrb[ij]<-1
  colnames(matrb)=cn
  rownames(matrb)=rn
  return(matrb)
}

#############################
#####     JUAN VELA     #####
#############################
# Data. Input the population tables (pollination expected proportion-MEMM) as matrix here.
matingJVq<-as.matrix(read.table("./Datos/JV_pollenNetwork.txt",header=T,sep="\t",dec=".",na.strings="NA", row.names=1))
# Set mating links according to threshold values.
# Adjust threshold values in function {quant2bin}.
matingJV<-quant2bin(matingJVq)     # Binary matrix with threshold value
matingJVPJ<-quant2binPJ(matingJVq) # Quant. matrix with threshold value
# Data. Input the population tables (genotypes) as matrix here.
JV<-read.table("./Datos/JV_genot.txt",header=T,row.names=1,sep="\t",dec=".",na.strings="NA")
JV.genot<-df2genind(JV, sep="/", ncode=NULL, ind.names=NULL, loc.names=NULL, pop=NULL,                       # missing=NA, 
                    ploidy=2, type=c("codom","PA"))
# Coordinates
JVcoor<-read.table("./Datos/JV_xy.txt",header=T,row.names=1,sep="\t",dec=".",na.strings="NA")
JV.genot@other<-JVcoor

pp<-as.matrix(JVcoor[,1:2])  # Creates dataset for plotting on the map
n<-length(JVcoor[,1]) # Number of trees in dataset

plot(JV.genot@other$x,JV.genot@other$y, xlim=c(263600,263850),ylim=c(4043950,4044200),
     xlab='x', ylab='y', main='JV 0.3')  # This plots the map 
# pp2<-as.matrix(cbind(JVcoor[,1],0)) # Creates dataset for plotting along line

# Function to iterate
foo <- function (j) {
  zz<- pp[j,]
  for (i in 1:n) {
    if (i!=j)
      if (matingJV[j, i]!= 0)  # GIVE THRESHOLD VALUE HERE!!!
        curvedarrow(from=zz, to=pp[i,],
                    lwd=0.5, endhead=F, arr.width=0.05, arr.length=0.01, 
                    lcol=rgb(0,0,1,0.4), curve=0.35, 
                    arr.pos=1,arr.lwd=0.01, arr.adj=0.5)
  }
}
# Now iterate over all pollen sources. This is to avoid the error with the arrow function 
# aborting in the double loop.
for (i in 1:n) {
  foo(i)
}


###############################
#####     LA SAUCEDA 4    #####
###############################
# Data. Input the population tables (pollination expected proportion-MEMM) as matrix here.
matingSauc4q<-as.matrix(read.table("Sauc4_pollenNetwork.txt",header=T,sep="\t",dec=".",na.strings="NA",row.names=1))
# Set mating links according to threshold values.
# Adjust threshold values in function {quant2binPJ}.
matingSauc4<-quant2bin(matingSauc4q)     # Binary matrix with threshold value
# matingSauc4PJ<-quant2binPJ(matingSauc4q) # Quant. matrix with threshold value
# Data. Input the population tables (genotypes) as matrix here.
SAUC4<-read.table("SAUC4_genot.txt",header=T,row.names=1,sep="\t",dec=".",na.strings="NA")
SAUC4.genot<-df2genind(SAUC4, sep="/", ncode=NULL, ind.names=NULL, loc.names=NULL, pop=NULL, missing=NA,
                       ploidy=2, type=c("codom","PA"))
# Coordinates
SAUC4coor<-read.table("SAUC4_xy.txt",header=T,row.names=1,sep="\t",dec=".",na.strings="NA")
SAUC4.genot@other<-SAUC4coor

plot(SAUC4.genot@other$x,SAUC4.genot@other$y, xlim=c(266700,266950),ylim=c(4046450,4046900),
     xlab='x', ylab='y', main='SAUC4 0.5')  
# pp2<-as.matrix(cbind(JVcoor[,1],0)) # Creates dataset for plotting along line
pp<-as.matrix(SAUC4coor[,1:2])  # Creates dataset for plotting on the map
n<-length(SAUC4coor[,1]) # Number of trees in dataset

# Function to iterate
foo <- function (j) {
  zz<- pp[j,]
  for (i in 1:n) {
    if (i!=j)
      if (matingSauc4[j, i]!= 0)  # GIVE THRESHOLD VALUE HERE!!!
        curvedarrow(from=zz, to=pp[i,],
                    lwd=0.5, endhead=F, arr.width=0.05, arr.length=0.01, 
                    lcol=rgb(0,0,1,0.4), curve=0.25, 
                    arr.pos=1,arr.lwd=0.01, arr.adj=0.5)
  }
}
# Now iterate over all pollen sources. This is to avoid the error with the arrow function 
# aborting in the double loop.
for (i in 1:n) {
  foo(i)
}


##############################
#####     LA SAUCEDA     #####
##############################
# Data. Input the population tables as matrix here.
matingSaucq<-as.matrix(read.table("./Datos/Sauc_pollenNetwork.txt",header=T,sep="\t",dec=".",na.strings="NA",row.names=1))
# Set mating links according to threshold values.
# Adjust threshold values in function {quant2binPJ}.
matingSauc<-quant2bin(matingSaucq)     # Binary matrix with threshold value
matingSaucPJ<-quant2binPJ(matingSaucq) # Quant. matrix with threshold value
# Data. Input the population tables as matrix here.
SAUC<-read.table("./Datos/SAUC_genot.txt",header=T,row.names=1,sep="\t",dec=".",na.strings="NA")
SAUC.genot<-df2genind(SAUC, sep="/", ncode=NULL, ind.names=NULL, loc.names=NULL, 
                      pop=NULL, # missing=NA,
                      ploidy=2, type=c("codom","PA"))
# Coordinates
SAUCcoor<-read.table("./Datos/SAUC_xy.txt",header=T,row.names=1,sep="\t",dec=".",na.strings="NA")
SAUC.genot@other<-SAUCcoor

plot(SAUC.genot@other$x,SAUC.genot@other$y, xlim=c(266750,268100),ylim=c(4044400,4045200),
     xlab='x', ylab='y', main='S2 0.5') 

pp<-as.matrix(SAUCcoor[,1:2])  # Creates dataset for plotting on the map
n<-length(SAUCcoor[,1]) # Number of trees in dataset

# Function to iterate
foo <- function (j) {
  zz<- pp[j,]
  for (i in 1:n) {
    if (i!=j)
      if (matingSauc[j, i]!= 0)  # GIVE THRESHOLD VALUE HERE!!!
        curvedarrow(from=zz, to=pp[i,],
                    lwd=0.5, endhead=F, arr.width=0.05, arr.length=0.01, 
                    lcol=rgb(0,0,1,0.4), curve=0.25, 
                    arr.pos=1,arr.lwd=0.01, arr.adj=0.5)
  }
}
# Now iterate over all pollen sources. This is to avoid the error with the arrow function 
# aborting in the double loop.
for (i in 1:n) {
  foo(i)
}



#############################
#####     EL ZAPATO     #####
#############################
# Data. Input the population tables as matrix here.
matingZAPq<-as.matrix(read.table("./Datos/ZAP_pollenNetwork.txt",header=T,sep="\t",dec=".",na.strings="NA",row.names=1))
# Set mating links according to threshold values.
# Adjust threshold values in function {quant2binPJ}.
matingZAP<-quant2bin(matingZAPq)     # Binary matrix with threshold value
# matingZAPPJ<-quant2binPJ(matingZAPq) # Quant. matrix with threshold value
# Data. Input the population tables as matrix here.
ZAP<-read.table("./Datos/ZAP_genot.txt",header=T,row.names=1,sep="\t",dec=".",na.strings="NA")
ZAP.genot<-df2genind(ZAP, sep="/", ncode=NULL, ind.names=NULL, loc.names=NULL, pop=NULL, missing=NA,
                     ploidy=2, type=c("codom","PA"))
# Coordinates
ZAPcoor<-read.table("./Datos/ZAP_xy.txt",header=T,row.names=1,sep="\t",dec=".",na.strings="NA")
ZAP.genot@other<-ZAPcoor

pp<-as.matrix(ZAPcoor[,1:2])  # Creates dataset for plotting on the map
n<-length(ZAPcoor[,1]) # Number of trees in dataset

plot(ZAP.genot@other$x,ZAP.genot@other$y, xlim=c(265200,266600),ylim=c(4040400,4040850),
     xlab='x', ylab='y', main='ZP 0.5')

# Function to iterate

foo <- function (j) {
  zz<- pp[j,]
  for (i in 1:n) {
    if (i!=j)
      if (matingZAP[j, i]!= 0)  # GIVE THRESHOLD VALUE HERE!!!
        curvedarrow(from=zz, to=pp[i,],
                    lwd=0.5, endhead=F, arr.width=0.05, arr.length=0.01, 
                    lcol=rgb(0,0,1,0.4), curve=0.3, 
                    arr.pos=1,arr.lwd=0.1, arr.adj=0.5)
  }
}
for (i in 1:n) {
  foo(i)
}
