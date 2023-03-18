############################################################################
# Analysis of mating networks. Frangula trees. 
# Example script with JV dataset.
# Pedro. 16 Jan 2012, Sevilla.
#---------------------------------------------------------------------------
setwd("~/Desktop/PhDRocio/MS/MS_Fa_mating_patterns/Datos")
library(network);library(sna);library(bipartite);library(igraph);


#---------------------------------------------------------------------------
# FUNCTIONS
#---------------------------------------------------------------------------
#Transformation of quantitative matrix into binary
#source("~/R/diego/R/quant2bin.R")
quant2bin<-function(matr)
{
  cn=colnames(matr)
  rn=rownames(matr)
  ij<-which(matr>=0.10,arr.ind=T) # Adjust here the threshold value for 
  # mating.
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
  ij<-which(matr>=0.50,arr.ind=T) # Adjust here the threshold value for                                   
  # mating.
  #  ij<-which(matr!=0,arr.ind=T)
  matrb<-matrix(0,dim(matr)[1],dim(matr)[2])
  matrb[ij]<-matr[ij]
  #  matrb[ij]<-1
  colnames(matrb)=cn
  rownames(matrb)=rn
  return(matrb)
}
sortmatr<-function(m){
  rsum<-rowSums(m)
  csum<-colSums(m)
  totsum<-sum(m)
  mext<-rbind(c(totsum,csum),cbind(rsum,m))
  i<-dim(mext)[1]
  j<-dim(mext)[2]
  mexts<-t(mext[order(-mext[,1]),])
  mexts<-t(mexts[order(-mexts[,1]),])
  ms<-mexts[2:i,2:j]
  return(ms)
}
# Plot matrix
plotmat<-function(matr,xlabel="Maternal trees",ylabel="Pollen donors",cexmin=0.2,cexmax=4,sortm=TRUE,sqroot=TRUE){
  r=nrow(matr)
  c=ncol(matr)
  if (sortm==T) matr=sortmatr(matr)
  if (sqroot==T) matr=sqrt(matr) 
  matplot(1,c,type="p",pch=21,col= rgb(0, 0, 1, 0.2),bg=rgb(0, 0, 1, 0.2),cex=cexmin+(cexmax*matr[1,1]/max(matr)),xlim=c(1,r),
          ylim=c(1,c),xlab=xlabel,ylab=ylabel,bty="n",xaxt="n",yaxt="n")
  for (i in 1:r){
    for (j in 1:c){
      if(matr[i,j]!=0){
        matplot(i,1+(c-j),type="p",pch=21,col= rgb(0, 0, 1, 0.1),bg=rgb(0, 0, 1, 0.3),cex=cexmin+(cexmax*matr[i,j]/max(matr)),add=T)
      }
    }
  }
}
#---------------------------------------------------------------------------

#############################
#####     JUAN VELA     #####
#############################

# Data. Input the population tables as matrix here.
matingJVq<-as.matrix(read.table("pollenNetwork_matrices/JV_pollenNetwork_median.txt",header=T,sep="\t",dec=".",na.strings="NA", row.names=1))
                                #col.names=c("JV0853","JV0854","JV0852","JV0851","JV0855","JV0624","JV0622","JV0619","JV0620","JV0621","JV0618","JV0617","JV0616","JV0615","JV0614","JV0613","JV0612","JV0611","JV0610","JV0609","JV0608","JV0607","JV0604","JV0606","JV0605","JV0603","JV0602","JV0601","JV0800","JV0799","JV0798","JV0797","JV0796","JV0793","JV0792","JV0791","JV0788","JV0789","JV0787","JV0786")))

# JV labels. All trees.
#row.names(matingJVq)<-c("JV0853","JV0854","JV0852","JV0851","JV0855","JV0624","JV0622","JV0619","JV0620","JV0621","JV0618","JV0617","JV0616","JV0615","JV0614","JV0613","JV0612","JV0611","JV0610","JV0609","JV0608","JV0607","JV0604","JV0606","JV0605","JV0603","JV0602","JV0601","JV0800","JV0799","JV0798","JV0797","JV0796","JV0793","JV0792","JV0791","JV0788","JV0789","JV0787","JV0786")

# JV. colnames for estimated MEMM proportions on sampled maternal trees
#col.names("JV0605","JV0607","JV0612","JV0624","JV0786","JV0787","JV0798","JV0853")

# Set mating links according to threshold values.
# Adjust threshold values in function {quant2binPJ}.
matingJV<-quant2bin(matingJVq)     # Binary matrix with threshold value
matingJVPJ<-quant2binPJ(matingJVq) # Quant. matrix with threshold value

# Assing the matrix of the focal population
mymat<-matingJVPJ
# Initialize network
mat<-network.initialize(dim(mymat)[1], directed=F)
# mat<-network.bipartite(as.matrix(mymat),mat)  # Just for bipartite graphs
mat<-as.network.matrix(mymat,matrix.type="adjacency", directed=F)
#---------------------------------------------------------------------------
plotmat(matingJVq,sortm=F)
plotweb(matingJVPJ, method="normal")
plot.network(mat,
    usearrows=F,jitter= T,
#   mode= "circle",
    mode= "fruchtermanreingold",
#   mode= "kamadakawai", 
    label=network.vertex.names(mat),displaylabels = T,
    boxed.labels= F,label.pad=0,label.pos= 5,label.cex= 0.8,
    vertex.col= rgb(0,0,1,0.6),vertex.cex= 1.5, vertex.sides= 6, vertex.lty= 0,
    edge.lty=0.8,edge.col= rgb(0,0,1,0.8),
    label.lty=NULL,usecurve = F)

# The thresholded matrix
plot.sociomatrix(mat,cex.lab=0.4)
par(mfrow=c(1,2))
heatmap(as.matrix(matingJVq),col=grey((64:0)/64),margins=c(15,15))
heatmap(as.matrix(mat),col=grey((32:0)/32),margins=c(15,15))
par(mfrow=c(1,1))

# Network analsyses
visweb(matingJVq, type="nested",text="interaction",plotsize=40,textsize=0.8)
visweb(mat, type="nested",plotsize=40,textsize=0.5, labsize=0.5)
#---------------------------------------------------------------------------

###############################
#####     LA SAUCEDA 4    #####
###############################

# Data. Input the population tables as matrix here.
matingSauc4q<-as.matrix(read.table("pollenNetwork_matrices/SAUC_pollenNetwork_median.txt",header=T,sep="\t",dec=".",na.strings="NA",row.names=1)) #col.names=c("4S0238","4S0239","4S0237","4S0236","4S0235","4S0234","4S0233","4S0232","4S0231","4S0230","4S0229","4S0228","4S0227","4S0225","4S0226","4S0224","4S0223","4S0222","4S0221","4S0220","4S0219","4S0218","4S0217","4S0216","4S0215","4S0214","4S0213","4S0211","4S0210","4S0882","4S0883","4S0884","4S0207","4S0203","4S0202","4S0206","4S0206bis","4S0208","4S0204","4S0205","4S0201","4S0200","4S0199","4S0198","4S0196","4S0881","4S0195","4S0194","4S0197","4S0193","4S0192","4S0191","4S0190","4S0188","4S0187","4S0186","4S0184","4S0189","4S0182","4S0185","4S0183","4S0181","4S0180","4S0178","4S0179","4S0879","4S0880","4S0176","4S0175","4S0169","4S0168","4S0167","4S0170","4S0166","4S0165","4S0164","4S0878","4S0163","4S0162","4S0160","4S0159","4S0161","4S0877","4S0171","4S0172","4S0173","4S0174","4S0158")))

# Sauc4 labels. All trees.
#row.names(matingSauc4q)<-c("4S0238","4S0239","4S0237","4S0236","4S0235","4S0234","4S0233","4S0232","4S0231","4S0230","4S0229","4S0228","4S0227","4S0225","4S0226","4S0224","4S0223","4S0222","4S0221","4S0220","4S0219","4S0218","4S0217","4S0216","4S0215","4S0214","4S0213","4S0211","4S0210","4S0882","4S0883","4S0884","4S0207","4S0203","4S0202","4S0206","4S0206bis","4S0208","4S0204","4S0205","4S0201","4S0200","4S0199","4S0198","4S0196","4S0881","4S0195","4S0194","4S0197","4S0193","4S0192","4S0191","4S0190","4S0188","4S0187","4S0186","4S0184","4S0189","4S0182","4S0185","4S0183","4S0181","4S0180","4S0178","4S0179","4S0879","4S0880","4S0176","4S0175","4S0169","4S0168","4S0167","4S0170","4S0166","4S0165","4S0164","4S0878","4S0163","4S0162","4S0160","4S0159","4S0161","4S0877","4S0171","4S0172","4S0173","4S0174","4S0158")

# Sauc4. colnames for estimated MEMM proportions on sampled maternal trees
#col.names("JV0605","JV0607","JV0612","JV0624","JV0786","JV0787","JV0798","JV0853")

# Set mating links according to threshold values.
# Adjust threshold values in function {quant2binPJ}.
matingSauc4<-quant2bin(matingSauc4q)     # Binary matrix with threshold value
matingSauc4PJ<-quant2binPJ(matingSauc4q) # Quant. matrix with threshold value

# Assing the matrix of the focal population
mymat<-matingSauc4PJ
# Initialize network
mat<-network.initialize(dim(mymat)[1], directed=F)
# mat<-network.bipartite(as.matrix(mymat),mat)  # Just for bipartite graphs
mat<-as.network.matrix(mymat,matrix.type="adjacency", directed=F)
#---------------------------------------------------------------------------
plotmat(matingSauc4q,sortm=F)
plotweb(matingSauc4PJ, method="normal")
plot.network(mat,
             usearrows=F,jitter= T,
             #   mode= "circle",
             mode= "fruchtermanreingold",
             #   mode= "kamadakawai", 
             label=network.vertex.names(mat),displaylabels = T,
             boxed.labels= F,label.pad=0,label.pos= 5,label.cex= 0.8,
             vertex.col= rgb(0,0,1,0.6),vertex.cex= 1.5, vertex.sides= 6, vertex.lty= 0,
             edge.lty=0.8,edge.col= rgb(0,0,1,0.8),
             label.lty=NULL,usecurve = F)

# The thresholded matrix
plot.sociomatrix(mat,cex.lab=0.4)
par(mfrow=c(1,2))
heatmap(as.matrix(matingSauc4q),col=grey((64:0)/64),margins=c(15,15))
heatmap(as.matrix(mat),col=grey((32:0)/32),margins=c(15,15))
par(mfrow=c(1,1))

# Network analsyses
visweb(matingSauc4q, type="nested",text="interaction",plotsize=40,textsize=0.8)
visweb(mat, type="nested",plotsize=40,textsize=0.5, labsize=0.5)
#---------------------------------------------------------------------------

##############################
#####     Pasada Llana     #####
##############################

# Data. Input the population tables as matrix here.
matingSaucq<-as.matrix(read.table("pollenNetwork_matrices/PL_pollenNetwork_median.txt",header=T,sep="\t",dec=".",na.strings="NA",row.names=1)) #col.names=c("S0038","S0039","S0040","S0041","S0042","S0043","S0044","S0842","S0843","S0045","S0046","S0844","S0846","S0848","S0847","S0047","S0048","S0049","S0050","S0051","S0052","S0849","S0054","S0056","S0055","S0053","S0060","S0059","S0058","S0057","S0062","S0061","S0063","S0064","S0066","S0065","S0857","S0858","S0856","S0067","S0859","S0068","S0069","S0070","S0860","S0071","S0072","S0073","S0074","S0075","S0077","S000A","S000B","S0079","S0076","S0078","S0865","S0081","S0082","S0083","S0084","S0864","S0085","S0148","S0147","S0146","S0145","S0144","S0143","S0149","S0141","S0142","S0140","S0139","S0137","S0136","S0138","S0135","S0087","S0088","S0133","S0089","S0090","S0132","S0131","S0130","S0091","S0092","S0867","S0866","S0868","S0095","S0094","S0093","S0129","S0128","S0127","S0096","S0097","S0098","S0099","S0100","S0101","S0126","S0102","S0103","S0105","S0104","S0125","S0106","S0107","S0124","S0123","S0876","S0108","S0109","S0110","S0122","S0111","S0869","S0119","S0112","S0113","S0114","S0871","S0121","S0870","S0120","S0875","S0116","S0115","S0872","S0873","S0874","S0118","S0117","S0155","S0154","S0153","S0152","S0151","S0150","S0156","S0157")))

# Sauc labels. All trees.
#row.names(matingSaucq)<-c("S0038","S0039","S0040","S0041","S0042","S0043","S0044","S0842","S0843","S0045","S0046","S0844","S0846","S0848","S0847","S0047","S0048","S0049","S0050","S0051","S0052","S0849","S0054","S0056","S0055","S0053","S0060","S0059","S0058","S0057","S0062","S0061","S0063","S0064","S0066","S0065","S0857","S0858","S0856","S0067","S0859","S0068","S0069","S0070","S0860","S0071","S0072","S0073","S0074","S0075","S0077","S000A","S000B","S0079","S0076","S0078","S0865","S0081","S0082","S0083","S0084","S0864","S0085","S0148","S0147","S0146","S0145","S0144","S0143","S0149","S0141","S0142","S0140","S0139","S0137","S0136","S0138","S0135","S0087","S0088","S0133","S0089","S0090","S0132","S0131","S0130","S0091","S0092","S0867","S0866","S0868","S0095","S0094","S0093","S0129","S0128","S0127","S0096","S0097","S0098","S0099","S0100","S0101","S0126","S0102","S0103","S0105","S0104","S0125","S0106","S0107","S0124","S0123","S0876","S0108","S0109","S0110","S0122","S0111","S0869","S0119","S0112","S0113","S0114","S0871","S0121","S0870","S0120","S0875","S0116","S0115","S0872","S0873","S0874","S0118","S0117","S0155","S0154","S0153","S0152","S0151","S0150","S0156","S0157")

# Sauc. colnames for estimated MEMM proportions on sampled maternal trees
#col.names("JV0605","JV0607","JV0612","JV0624","JV0786","JV0787","JV0798","JV0853")

# Set mating links according to threshold values.
# Adjust threshold values in function {quant2binPJ}.
matingSauc<-quant2bin(matingSaucq)     # Binary matrix with threshold value
matingSaucPJ<-quant2binPJ(matingSaucq) # Quant. matrix with threshold value

# Assing the matrix of the focal population
mymat<-matingSaucPJ
# Initialize network
mat<-network.initialize(dim(mymat)[1], directed=F)
# mat<-network.bipartite(as.matrix(mymat),mat)  # Just for bipartite graphs
mat<-as.network.matrix(mymat,matrix.type="adjacency", directed=F)
#---------------------------------------------------------------------------
plotmat(matingSaucq,sortm=F)
plotweb(matingSaucPJ, method="normal")
plot.network(mat,
             usearrows=F,jitter= T,
             #   mode= "circle",
             mode= "fruchtermanreingold",
             #   mode= "kamadakawai", 
             label=network.vertex.names(mat),displaylabels = T,
             boxed.labels= F,label.pad=0,label.pos= 5,label.cex= 0.8,
             vertex.col= rgb(0,0,1,0.6),vertex.cex= 1.5, vertex.sides= 6, vertex.lty= 0,
             edge.lty=0.8,edge.col= rgb(0,0,1,0.8),
             label.lty=NULL,usecurve = F)

# The thresholded matrix
plot.sociomatrix(mat,cex.lab=0.4)
heatmap(as.matrix(matingSaucq),col=grey((64:0)/64),margins=c(15,15))
heatmap(as.matrix(mat),col=grey((32:0)/32),margins=c(15,15))

# Network analsyses
visweb(matingSaucq, type="nested",text="interaction",plotsize=144,textsize=0.8)
visweb(mat, type="nested",plotsize=40,textsize=0.5, labsize=0.4)
#---------------------------------------------------------------------------

#############################
#####     EL ZAPATO     #####
#############################

# Data. Input the population tables as matrix here.
matingZAPq<-as.matrix(read.table("pollenNetwork_matrices/ZP_pollenNetwork_median.txt",header=T,sep="\t",dec=".",na.strings="NA",row.names=1))
#col.names=c("Z0557","Z0558","Z0559","Z0356","Z0357","Z0358","Z0560","Z0561","Z0562","Z0562mother","Z0565","Z0563","Z0564","Z0569","Z0568","Z0567","Z0566","Z0570","Z0571","Z0572","Z0573","Z0574","Z0575","Z0577","Z0576","Z0578","Z0579","Z0580","Z0581","Z0582","Z0583","Z0585","Z0586","Z0587","Z0588","Z0591","Z0589","Z0590","Z0592","Z0784","Z0594","Z0595","Z0596","Z0783","Z0597","Z0598","Z0599","Z0600","Z0301","Z0862","Z0861","Z0302","Z0303","Z0307","Z0306","Z0305","Z0304","Z0311","Z0308","Z0310","Z0309","Z0313","Z0317","Z0312","Z0318","Z0316","Z0314","Z0315","Z0319","Z0320","Z0321","Z0324","Z0322","Z0323","Z0325","Z0326","Z0327","Z0328","Z0329","Z0330","Z0331","Z0332","Z0333","Z0785","Z0334","Z0335","Z0336","Z0338","Z0337","Z0339","Z0340","Z0343","Z0341","Z0342","Z0344","Z0349","Z0359","Z0348","Z0347","Z0346","Z0345","Z0350","Z0353","Z0352","Z0351","Z0354","Z0355","Z0360","Z0361","Z0362","Z0364","Z0365","Z0366","Z0367","Z0368","Z0369","Z0370","Z0371","Z0372","Z0373","Z0374","Z0375","Z0376","Z0377","Z0378","Z0379","Z0380","Z0381","Z0382","Z0383","Z0384","Z0385","Z0386","Z0387","Z0388","Z0389","Z0390","Z0391","Z0392","Z0393","Z0863","Z0394","Z0395","Z0396","Z0398","Z0397","Z0399","Z0400","Z0701","Z0702","Z0703","Z0704","Z0705","Z0706","Z0707","Z0709","Z0708","Z0710","Z0711","Z0713","Z0712","Z0715","Z0714","Z0717","Z0716","Z0718","Z0719","Z0720","Z0721","Z0722","Z0731","Z0730","Z0723","Z0725","Z0724","Z0727","Z0726","Z0728","Z0729","Z0733","Z0734","Z0732","Z0735","Z0736","Z0737","Z0738","Z0739","Z0740","Z0741","Z0742","Z0743","Z0746","Z0748","Z0747","Z0744","Z0745","Z0751","Z0749","Z0750","Z0752","Z0753","Z0754","Z0755","Z0756","Z0757","Z0758","Z0759","Z0760","Z0761","Z0763","Z0762","Z0764","Z0768","Z0767","Z0765","Z0769","Z0771","Z0772","Z0774","Z0775","Z0776","Z0777","Z0781","Z0779","Z0780","Z0782")))

# ZAP labels. All trees.
#row.names(matingZAPq)<-c("Z0557","Z0558","Z0559","Z0356","Z0357","Z0358","Z0560","Z0561","Z0562","Z0562mother","Z0565","Z0563","Z0564","Z0569","Z0568","Z0567","Z0566","Z0570","Z0571","Z0572","Z0573","Z0574","Z0575","Z0577","Z0576","Z0578","Z0579","Z0580","Z0581","Z0582","Z0583","Z0585","Z0586","Z0587","Z0588","Z0591","Z0589","Z0590","Z0592","Z0784","Z0594","Z0595","Z0596","Z0783","Z0597","Z0598","Z0599","Z0600","Z0301","Z0862","Z0861","Z0302","Z0303","Z0307","Z0306","Z0305","Z0304","Z0311","Z0308","Z0310","Z0309","Z0313","Z0317","Z0312","Z0318","Z0316","Z0314","Z0315","Z0319","Z0320","Z0321","Z0324","Z0322","Z0323","Z0325","Z0326","Z0327","Z0328","Z0329","Z0330","Z0331","Z0332","Z0333","Z0785","Z0334","Z0335","Z0336","Z0338","Z0337","Z0339","Z0340","Z0343","Z0341","Z0342","Z0344","Z0349","Z0359","Z0348","Z0347","Z0346","Z0345","Z0350","Z0353","Z0352","Z0351","Z0354","Z0355","Z0360","Z0361","Z0362","Z0364","Z0365","Z0366","Z0367","Z0368","Z0369","Z0370","Z0371","Z0372","Z0373","Z0374","Z0375","Z0376","Z0377","Z0378","Z0379","Z0380","Z0381","Z0382","Z0383","Z0384","Z0385","Z0386","Z0387","Z0388","Z0389","Z0390","Z0391","Z0392","Z0393","Z0863","Z0394","Z0395","Z0396","Z0398","Z0397","Z0399","Z0400","Z0701","Z0702","Z0703","Z0704","Z0705","Z0706","Z0707","Z0709","Z0708","Z0710","Z0711","Z0713","Z0712","Z0715","Z0714","Z0717","Z0716","Z0718","Z0719","Z0720","Z0721","Z0722","Z0731","Z0730","Z0723","Z0725","Z0724","Z0727","Z0726","Z0728","Z0729","Z0733","Z0734","Z0732","Z0735","Z0736","Z0737","Z0738","Z0739","Z0740","Z0741","Z0742","Z0743","Z0746","Z0748","Z0747","Z0744","Z0745","Z0751","Z0749","Z0750","Z0752","Z0753","Z0754","Z0755","Z0756","Z0757","Z0758","Z0759","Z0760","Z0761","Z0763","Z0762","Z0764","Z0768","Z0767","Z0765","Z0769","Z0771","Z0772","Z0774","Z0775","Z0776","Z0777","Z0781","Z0779","Z0780","Z0782")

# ZAP. colnames for estimated MEMM proportions on sampled maternal trees
#col.names("JV0605","JV0607","JV0612","JV0624","JV0786","JV0787","JV0798","JV0853")

# Set mating links according to threshold values.
# Adjust threshold values in function {quant2binPJ}.
matingZAP<-quant2bin(matingZAPq)     # Binary matrix with threshold value
matingZAPPJ<-quant2binPJ(matingZAPq) # Quant. matrix with threshold value

# Assing the matrix of the focal population
mymat<-matingZAPPJ
# Initialize network
mat<-network.initialize(dim(mymat)[1], directed=F)
# mat<-network.bipartite(as.matrix(mymat),mat)  # Just for bipartite graphs
mat<-as.network.matrix(mymat,matrix.type="adjacency", directed=F)
#---------------------------------------------------------------------------
plotmat(matingZAPq,sortm=F)
plotweb(matingZAPPJ, method="normal")
plot.network(mat,
             usearrows=F,jitter= T,
             #   mode= "circle",
             mode= "fruchtermanreingold",
             #   mode= "kamadakawai", 
             label=network.vertex.names(mat),displaylabels = T,
             boxed.labels= F,label.pad=0,label.pos= 5,label.cex= 0.8,
             vertex.col= rgb(0,0,1,0.6),vertex.cex= 1.5, vertex.sides= 6, vertex.lty= 0,
             edge.lty=0.8,edge.col= rgb(0,0,1,0.8),
             label.lty=NULL,usecurve = F)

# The thresholded matrix
plot.sociomatrix(mat,cex.lab=0.4)
par(mfrow=c(1,2))
heatmap(as.matrix(matingZAPq),col=grey((64:0)/64),margins=c(15,15))
heatmap(as.matrix(mat),col=grey((32:0)/32),margins=c(15,15))

# Network analsyses
visweb(matingZAPq, type="nested",text="interaction",plotsize=40,textsize=0.8)
visweb(mat, type="nested",plotsize=40,textsize=0.5, labsize=0.5)
#---------------------------------------------------------------------------

