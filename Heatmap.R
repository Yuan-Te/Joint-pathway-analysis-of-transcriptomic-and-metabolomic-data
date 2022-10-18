# Heatmap
library("reshape") 
library("ggplot2")        
library(pheatmap)

#ggplot 246 metabolites-------------------------------------------
metabolite_246 <- read.csv("/Users/LaiYuanTe/Desktop/Term2/Research/Heatmap/246_metabolites/heatmap_246.csv",check.names=FALSE)
class(metabolite_246) 
dim(metabolite_246)#120 x 248
sapply(metabolite_246, class) # Examinate every value is numaric
metabolite_246[1:5,1:5]

metabolite_246 <- scale(metabolite_246[,3:248]) # Scaling every metabolites, data will become matrix
class(metabolite_246)
dim(metabolite_246) # Matrix 120 x 246
metabolite_246 <- as.data.frame(metabolite_246) #Turn Matrix to Dataframe
class(metabolite_246) #dataframe 120 x 246
metabolite_246[1:5,1:5]

#Add back patient name and group name
name <- read.csv("/Users/LaiYuanTe/Desktop/Term2/Research/Heatmap/246_metabolites/heatmap_246.csv",check.names=FALSE)
metabolite_246 <- cbind(Group = name$GROUP_NAME, metabolite_246) #Add to first column
metabolite_246 <- cbind(Sample = name$PATIENT_SAMPLE_NAME, metabolite_246) #Add to very first column 
dim(metabolite_246) #120 x 248
sapply(metabolite_246, class)  # Examinate again
metabolite_246[1:5,1:5]

#https://jcoliver.github.io/learn-r/006-heatmaps.html
#transform to table
#install.packages("tidyr")
library("tidyr")
metabolite.long <- pivot_longer(data = metabolite_246, 
                                cols = -c(1:2), 
                                names_to = "Metabolites", 
                                values_to = "Abundance")
head(metabolite.long)


#plot with ggplot
library(ggplot2)
metabolite.heatmap <- ggplot(data = metabolite.long, mapping = aes(x = Sample,
                                                                   y = Metabolites, fill = Abundance))  +
  geom_tile() + ggtitle("246 Metabolites") + theme(plot.title = element_text(size=12, hjust = 0.5), 
                                                    axis.title = element_text(size = 10), 
                                                    axis.text = element_text(size = 1),
                                                    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                                                    legend.text = element_text(size = 8)) + 
                                              xlab(label = "Sample") + 
                                              scale_fill_gradient2(low = "#132F5D",
                                                                    mid = "#FFFFFF",
                                                                    high = "#66001F", lim=c(-6, 6)) +
                                              facet_grid(~ Group, switch = "x", scales = "free_x", space = "free_x") +
                                              theme(strip.text.x = element_text(size=5, hjust = 0.5))
metabolite.heatmap
---------------------------------------------------------------------------
#Average
average_246 <- read.csv("/Users/LaiYuanTe/Desktop/Term2/Research/Heatmap/246_metabolites/heatmap_average_246.csv",check.names=FALSE)
class(average_246) 
dim(average_246)#120 x 248
sapply(average_246, class) # Examinate every value is numaric
average_246[,1:5]

average_246 <- scale(average_246[,2:247]) # Scaling every metabolites, data will become matrix
class(average_246)
dim(average_246) # Matrix 4 x 246
average_246 <- as.data.frame(average_246) #Turn Matrix to Dataframe
class(average_246) #dataframe 4 x 246
metabolite_246[,1:5]

#Add back patient name and group name
name <- read.csv("/Users/LaiYuanTe/Desktop/Term2/Research/Heatmap/246_metabolites/heatmap_average_246.csv",check.names=FALSE)
average_246 <- cbind(Group = name$GROUP_NAME, average_246) #Add to first column
dim(average_246) #4 x 247
sapply(average_246, class)  # Examinate again
average_246[,1:5]

#https://jcoliver.github.io/learn-r/006-heatmaps.html
#transform to table
#install.packages("tidyr")
library("tidyr")
metabolite.long <- pivot_longer(data = average_246, 
                                cols = -c(1), 
                                names_to = "Metabolites", 
                                values_to = "Abundance")
head(metabolite.long)


#plot with ggplot
library(ggplot2)
metabolite.heatmap <- ggplot(data = metabolite.long, mapping = aes(x = Group,
                                                                   y = Metabolites, fill = Abundance))  +
  geom_tile() + ggtitle("246 Metabolites") + theme(plot.title = element_text(size=12, hjust = 0.5), 
                                                   axis.title = element_text(size = 10), 
                                                   axis.text = element_text(size = 1),
                                                   axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                                                   legend.text = element_text(size = 8)) + 
  xlab(label = "Group") + 
  scale_fill_gradient2(low = "#132F5D",
                       mid = "#FFFFFF",
                       high = "#66001F", lim=c(-3, 3)) +
  facet_grid(~ Group, switch = "x", scales = "free_x", space = "free_x") +
  theme(strip.text.x = element_text(size=5, hjust = 0.5))
metabolite.heatmap
----------------------------------------------------------------------------
# Different melting method with euclidean cluster
# Distance
data <-  metabolite_246
data[1:5,1:5]
colnames(data)
dim(data)
data <- scale(data[,3:248])
ord <- hclust( dist(data, method = "euclidean"), method = "ward.D" )$order
ord

pd <- as.data.frame( data )
pd[1:5,1:5]
dim(pd)
#Add back patient name and group name
name <- read.csv("/Users/LaiYuanTe/Desktop/Term2/Research/Heatmap/246_metabolites/heatmap_246.csv",check.names=FALSE)
pd <- cbind(Group = name$GROUP_NAME, pd) #Add to first column
pd <- cbind(Sample = name$PATIENT_SAMPLE_NAME, pd) #Add to very first column 
dim(pd) #120 x 248
sapply(metabolite_246, class)  # Examinate again
pd[1:5,1:5]

pd.m <- melt( pd, id.vars = c("Sample", "Group"), variable = "Metabolites" )

pd.m$Metabolites <- factor( pd.m$Metabolites, levels = colnames(data), labels = seq_along( colnames(data) ) )
pd.m$Sample <- factor( pd.m$Sample, levels = rownames(data)[ord])

#plot with ggplot
library(ggplot2)
metabolite.heatmap <- ggplot( pd.m, aes(Sample, Metabolites) ) +
                      geom_tile(aes(fill = value)) +
                      scale_fill_gradient2(low = "#132F5D", high = "#66001F")
metabolite.heatmap

metabolite.heatmap <- ggplot(data = pd.m, mapping = aes(x = Sample,
                                                        y = Metabolites, fill = value))  +
  geom_tile() + ggtitle("246 Metabolites") + theme(plot.title = element_text(size=12, hjust = 0.5), 
                                                   axis.title = element_text(size = 10), 
                                                   axis.text = element_text(size = 1),
                                                   axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                                                   legend.text = element_text(size = 8)) + 
  xlab(label = "Sample") + 
  scale_fill_gradient2(low = "#132F5D",
                       mid = "#FFFFFF",
                       high = "#66001F", lim=c(-10, 10)) +
  facet_grid(~ Group, switch = "x", scales = "free_x", space = "free_x") +
  theme(strip.text.x = element_text(size=5, hjust = 0.5))
metabolite.heatmap


#ggplot 393 metabolites-------------------------------------------
delim = ","  # or is it "\t" ?
dec = "."    # or is it "," ?
metabolite_393 <- read.csv("/Users/LaiYuanTe/Desktop/Term2/Research/Heatmap/heatmap_393.csv", header=TRUE, sep=delim, dec=dec, stringsAsFactors=FALSE)
class(metabolite_393) 
dim(metabolite_393)#120 x 395
sapply(metabolite_393, class) # Examinate every value is numaric
metabolite_393$palmitoleamide..16.1.. # Strangely high

metabolite_393 <- scale(metabolite_393[,3:395]) # Scaling every metabolites, data will become matrix
class(metabolite_393)
dim(metabolite_393) # Matrix 120 x 393
metabolite_393 <- as.data.frame(metabolite_393) #Turn Matrix to Dataframe
class(metabolite_393) #dataframe 120 x 393

#Add back patient name and group name
name <- read.csv("/Users/LaiYuanTe/Desktop/Term2/Research/Heatmap/heatmap_393.csv", header=TRUE, sep=delim, dec=dec, stringsAsFactors=FALSE)
metabolite_393 <- data.frame(Group = name$GROUP_NAME, metabolite_393) #Add to first column
metabolite_393 <- data.frame(Sample = name$PATIENT_SAMPLE_NAME, metabolite_393) #Add to very first column 
dim(metabolite_393) #120 x 393
sapply(metabolite_393, class)  # Examinate again
metabolite_393$palmitoleamide..16.1.. # Range from -3 to 3

#https://jcoliver.github.io/learn-r/006-heatmaps.html
#transform to table
#install.packages("tidyr")
library("tidyr")
metabolite.long <- pivot_longer(data = metabolite_393, 
                          cols = -c(1:2), 
                          names_to = "Metabolites", 
                          values_to = "Abundance")
head(metabolite.long)

#plot with ggplot
library(ggplot2)
metabolite.heatmap <- ggplot(data = metabolite.long, mapping = aes(x = Sample,
                                                       y = Metabolites, fill = Abundance))  +
                              geom_tile() + 
                              ggtitle("393 Metabolites") +
                              theme(plot.title = element_text(size=12, hjust = 0.5), 
                                    axis.title = element_text(size = 10), 
                                    axis.text = element_text(size = 1),
                                    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                                    legend.text = element_text(size = 8)) + 
                              xlab(label = "Sample") + 
                              scale_fill_gradient2(low = "#08264D",
                                                   mid = "#FFFFFF",
                                                   high = "#FF0000") +
                              facet_grid(~ Group, switch = "x", scales = "free_x", space = "free_x") +
                              theme(strip.text.x = element_text(size=5, hjust = 0.5))
metabolite.heatmap











#Pheatmap-------------------------------------
#heatmap_393 <- read.csv("/Users/LaiYuanTe/Desktop/Term2/Research/Heatmap/heatmap_393_transpose.csv",check.names=FALSE)
#class(heatmap_393)
#heatmap_393 <- heatmap_393[-1,] #delete the first row
#rownames(heatmap_393) <- heatmap_393[,1]
#heatmap_393 <- heatmap_393[,-1]  #delete the second row
#which(is.na(heatmap_393))

#sapply(heatmap_393, class) 
#heatmap_393_num <- as.data.frame(apply(heatmap_393, 2, as.numeric)) #change value from character to numeric
#sapply(heatmap_393_num, class) 

# Plot the heatmap
#pheatmap(heatmap_393_num, scale = "row", color = c("blue", "light blue", "white", "white", "yellow", "orange"), breaks = c(-2, -1, -0.5, 0, 0.5, 1, 2))

#pheatmap(heatmap_393_num, legend_breaks = c(5, 0, -5))



# Heatmap extraction(246)
Batch_norm <- read.csv("/Users/LaiYuanTe/Desktop/Term2/Research/Raw_data/Batch_norm_name.csv",check.names=FALSE)
nrow(Batch_norm)
ncol(Batch_norm)
Batch_norm[1:5, 1:5]
class(Batch_norm)

vol <- read.csv("/Users/LaiYuanTe/Desktop/Term2/Research/Heatmap/List of 246 metabolites.csv")
nrow(vol)
class(vol)
length(vol$Metabolites)

heatmap <- data.frame(Patient = Batch_norm[,1], Group = Batch_norm[,2])

a <- which(colnames(Batch_norm) == vol$Metabolites[1])

total_index <- c()
total_name <- c()
for (i in 1:length(vol$Metabolites)){
  index <- which(colnames(Batch_norm) == vol$Metabolites[i])
  total_index <- c(total_index, index)
  total_name <- c(total_name, vol$Metabolites[i])
}

Batch_norm <- Batch_norm[,total_index]
nrow(Batch_norm)
ncol(Batch_norm)

rownames(Batch_norm) <- c(1:120)

Batch_norm$PARENT_SAMPLE_NAME=heatmap[,1]
Batch_norm$GROUP_NAME=heatmap[,2]
nrow(Batch_norm)
ncol(Batch_norm)

write.csv(Batch_norm, file = "/Users/LaiYuanTe/Desktop/Term2/Research/Heatmap/heatmap_246.csv",row.names = FALSE)



# Heatmap extraction (393 metabolites)
Batch_norm <- read.csv("/Users/LaiYuanTe/Desktop/Term2/Research/Raw_data/Batch_norm_name.csv",check.names=FALSE)
nrow(Batch_norm)
ncol(Batch_norm)
Batch_norm[1:5, 1:5]


vol <- read.csv("/Users/LaiYuanTe/Desktop/Term2/Research/Heatmap/Name of repetitive removed metabolites.csv")
nrow(vol)
class(vol)
length(vol$All_metabolites)

#heatmap <- data.frame(PATIENT_SAMPLE_NAME=Batch_norm$PATIENT_SAMPLE_NAME,GROUP_NAME=Batch_norm$GROUP_NAME)
#rownames(heatmap)

a <- which(colnames(Batch_norm) == vol$All_metabolites[1])

total_index <- c()
total_name <- c()

for (i in 1:length(vol$All_metabolites)){
  index <- which(colnames(Batch_norm) == vol$All_metabolites[i])
  total_index <- c(total_index, index)
  total_name <- c(total_name, vol$All_metabolites[i])
}

Batch_norm <- Batch_norm[,total_index]
nrow(Batch_norm)
ncol(Batch_norm)

rownames(Batch_norm) <- c(1:120)

Batch_norm$PARENT_SAMPLE_NAME=heatmap$PARENT_SAMPLE_NAME
Batch_norm$GROUP_NAME=heatmap$GROUP_NAME
nrow(Batch_norm)
ncol(Batch_norm)

write.csv(Batch_norm, file = "/Users/LaiYuanTe/Desktop/Term2/Research/Heatmap/heatmap_393_1.csv",row.names = FALSE)


