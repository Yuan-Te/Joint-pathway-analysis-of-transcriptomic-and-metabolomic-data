#PCA using MixOmics from Bioconductor
# https://mixomicsteam.github.io/Bookdown/plsda.html#tuning:sPLSDA

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("mixOmics")

library(mixOmics)
library(ggplot2)
library(dplyr)

setwd("C:\\Users\\mayke\\OneDrive\\??????\\Jeff Research\\2. Raw_data")
data <- read.csv("Batch_norm_name.csv",check.names=FALSE)

metabolites <- data[,-2]
metabolites$PATIENT_SAMPLE_NAME <-gsub("DH-META-","",as.character(metabolites$PATIENT_SAMPLE_NAME))
metabolites$PATIENT_SAMPLE_NAME
rownames(metabolites) <- metabolites$PATIENT_SAMPLE_NAME
metabolites <- metabolites[,-1]

colnames(metabolites)[1001]
colnames(metabolites)[1002]
metabolites <- metabolites[,-1002:-1189]
length(colnames(metabolites))

MyResult.pca <- pca(metabolites)
plotIndiv(MyResult.pca)
plotVar(MyResult.pca, cutoff = 0.8) 
plotIndiv(MyResult.pca, group = data$GROUP_NAME, legend = TRUE)
#Outlier: 31123,31176,31196,31171,31211

#Remove outlier
which(rownames(metabolites)=="31123") #19
which(rownames(metabolites)=="31176") #72
which(rownames(metabolites)=="31196") #92
which(rownames(metabolites)=="31171") #67
which(rownames(metabolites)=="31211") #107
metabolites <- metabolites[-c(19,72,92,67,107),]
dim(metabolites)

data <- data[-c(19,72,92,67,107),]
dim(data)

MyResult.pca <- pca(metabolites)
plotIndiv(MyResult.pca)
plotVar(MyResult.pca, cutoff = 0.8) 
plotIndiv(MyResult.pca, group = data$GROUP_NAME, legend = TRUE)

#------------------------------------------------------------------
#PLS-DA
library(mixOmics)

data <- read.csv("Batch_norm_name.csv",check.names=FALSE)
data <- data[-(61:90),]

metabolites <- data[,-2]
metabolites$PATIENT_SAMPLE_NAME <-gsub("DH-META-","",as.character(metabolites$PATIENT_SAMPLE_NAME))
metabolites$PATIENT_SAMPLE_NAME
rownames(metabolites) <- metabolites$PATIENT_SAMPLE_NAME
metabolites <- metabolites[,-1]
colnames(metabolites)[1001]
colnames(metabolites)[1002]
metabolites <- metabolites[,-1002:-1189]
dim(metabolites)

X <- metabolites
Y <- data$GROUP_NAME
dim(X)
length(Y)

MyResult.splsda <- splsda(X, Y) #Use all 1001 metabolites
plotIndiv(MyResult.splsda)

MyResult.splsda <- splsda(X, Y, keepX = c(50,50)) # Use most 50 important variables
plotIndiv(MyResult.splsda)

plotVar(MyResult.splsda) 
selectVar(MyResult.splsda, comp=1)$name #50 important variables
selectVar(MyResult.splsda, comp=2)$name #50 important variables
Plsda_50_comp1 <- data.frame(selectVar(MyResult.splsda, comp=1)$name)
colnames(Plsda_50_comp1) <- "Comp_1"
Plsda_50_comp2 <- data.frame(selectVar(MyResult.splsda, comp=2)$name)
colnames(Plsda_50_comp2) <- "Comp_2"
Plsda_50 <- bind_cols(Plsda_50_comp1, Plsda_50_comp2)
#write.csv(Plsda_50, file = "/Users/LaiYuanTe/Desktop/Term2/Research/3. Volcanol and Pathway analysis/PCA and PLS_DA/Plsda_50.csv",row.names = FALSE)


#Total Pls-da
MyResult.splsda <- splsda(X, Y) #total metabolites
plotIndiv(MyResult.splsda, ind.names = FALSE, legend=TRUE, #Group with circle
          ellipse = TRUE, style="ggplot2", cex=0.3,col=c("indianred1","green4","royalblue1"), title = "PLS-DA on metabolomics data")
linn.vip <- vip(MyResult.splsda)
barplot(linn.vip,
        beside = TRUE, col = c("lightblue", "mistyrose", "lightcyan"),
        ylim = c(0, 3), legend = colnames(linn.vip),
        main = "Variable Importance in the Projection", font.main = 4)
#write.csv(linn.vip, "/Users/LaiYuanTe/Desktop/Term2/Research/3. Volcanol and Pathway analysis/PCA and PLS_DA/VIP_total.csv")

#Top 50 metabolites Pls-da
MyResult.splsda <- splsda(X, Y, keepX = c(50,50)) # Top 50 metabolites
plotIndiv(MyResult.splsda, ind.names = FALSE, legend=TRUE, #Group with circle
          ellipse = TRUE, style="ggplot2", cex=0.3,col=c("indianred1","green4","royalblue1"), title = "PLS-DA on metabolomics data")


plotVar(MyResult.splsda,var.names=FALSE) #No cutoff
plotVar(MyResult.splsda,var.names=FALSE,cutoff=0.5) #With cutoff

auc.plsda <- auroc(MyResult.splsda) #AUC courve

#Variable selection outputs
MyResult.splsda2 <- splsda(X,Y, ncomp=3)
selectVar(MyResult.splsda2, comp=1)$value #list of variables
selectVar(MyResult.splsda2, comp=2)$value #list of variables
selectVar(MyResult.splsda2, comp=3)$value #list of variables

plotLoadings(MyResult.splsda2, contrib = 'max', method = 'mean') #Name of metabolites are too long

#PLot 3D
#install.packages("rgl")
library("rgl")
plotIndiv(MyResult.splsda2, ind.names = FALSE, legend=TRUE, ellipse = TRUE,style="3d", cex=0.3, col=c("indianred1","green4","royalblue1"))
#HC: red
#SSC-CON: green
#PAH: blue


#-------------------------------------------------------------------------------------- 
#Decide how many component
MyResult.plsda2 <- splsda(X,Y, ncomp=10)
set.seed(30) # for reproducbility in this vignette, otherwise increase nrepeat
MyPerf.plsda <- perf(MyResult.plsda2, validation = "Mfold", folds = 3, 
                     progressBar = FALSE, nrepeat = 50) # we suggest nrepeat = 50
plot(MyPerf.plsda, col = color.mixo(5:7), sd = TRUE, legend.position = "horizontal") 
MyPerf.plsda

#Decide how many variables
list.keepX <- c(5:10,  seq(20, 100, 10))
set.seed(30) # for reproducbility in this vignette, otherwise increase nrepeat
tune.splsda.srbct <- tune.splsda(X, Y, ncomp = 4, # we suggest to push ncomp a bit more, e.g. 4
                                 validation = 'Mfold',
                                 folds = 3, dist = 'max.dist', progressBar = FALSE,
                                 measure = "BER", test.keepX = list.keepX,
                                 nrepeat = 50)   # we suggest nrepeat = 50
error <- tune.splsda.srbct$error.rate
ncomp <- tune.splsda.srbct$choice.ncomp$ncomp # optimal number of components based on t-tests on the error rate
ncomp
select.keepX <- tune.splsda.srbct$choice.keepX[1:ncomp]  # optimal number of variables to select
select.keepX 

#plot(tune.splsda.srbct, col = color.jet(ncomp))

#FINAL PLS-DA
MyResult.splsda.final <- splsda(X, Y, ncomp = ncomp, keepX = select.keepX)

plotIndiv(MyResult.splsda.final, ind.names = FALSE, legend=TRUE,
          ellipse = TRUE, title="PLS-DA on metabolomics data",cex=0.3, col=c("indianred1","green4","royalblue1"))






