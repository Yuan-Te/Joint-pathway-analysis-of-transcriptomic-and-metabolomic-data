# Multiple regression_three groups
library(ggplot2)
library(ggpubr)
library(fastDummies)
library(tidyverse)
library(caret)


data <- read.csv("/Users/LaiYuanTe/Desktop/Term2/Research/7. Patient information/Charateristic & 189 metabolites.csv", check.names = FALSE)
data[1:5,1:10]
dim(data) #18 variables and 189 metabolites and 1 disease and 1 patient's code

vital_metabolites <- read.csv("/Users/LaiYuanTe/Desktop/Term2/Research/7. Patient information/35 important metabolites.csv", check.names = FALSE)

#modify categorical variables in patient data
# Condition
data <- dummy_cols(data, select_columns = "Condition") #new cols are at last
# Gender
data$New_Gender <- ifelse(data$Gender == "M", 1, 0)
# Disease subset
data <- dummy_cols(data, select_columns = "Disease_subset")

# change name to be readable
colnames(data)[which(colnames(data)=="Condition_SSC CON")] = "Condition_SSC_CON"
colnames(data)[which(colnames(data)=="Condition_SSC PAH")] = "Condition_SSC_PAH"

#------------------------------------------------------------------------------------------
#Visualization
index <- which(colnames(data)=="1-methyl-4-imidazoleacetate")
#Condition
ggplot(data, aes(x=Condition, y=data[,index])) + 
  geom_boxplot()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  labs(title="1-methyl-4-imidazoleacetate",x="Condition", y = "Relative Abundance")
#Gender
ggplot(data, aes(x=Gender, y=data[,index])) + 
  geom_boxplot()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  labs(title="1-methyl-4-imidazoleacetate",x="Gender", y = "Relative Abundance")
#Age
ggscatter(data, x = "Age", y = "1-methyl-4-imidazoleacetate",
          add = "reg.line",                                 # Add regression line
          conf.int = TRUE,                                  # Add confidence interval
          add.params = list(color = "blue",
                            fill = "lightgray"))+
  stat_cor(method = "spearman", label.x = 18, label.y = 4) 


#Multiple (Linear) regression with age, gender and condition
# 1-methyl-4-imidazoleacetate
index <- which(colnames(data)=="1-methyl-4-imidazoleacetate")
index
df <- data.frame(
  metaboltite_name = data[,index],
  Age = data$Age,
  New_Gender = data$New_Gender,
  Condition_HC = data$Condition_HC,
  Condition_SSC_CON = data$Condition_SSC_CON,
  Condition_SSC_PAH = data$Condition_SSC_PAH
)
mean(df$metaboltite_name)    

model <- lm(metaboltite_name ~ Age + New_Gender + Condition_SSC_CON + Condition_SSC_PAH, data=df)
summary(model)
model <- lm(metaboltite_name ~ Age + Condition_SSC_PAH + Condition_SSC_CON, data=df)
summary(model)
model <- lm(metaboltite_name ~ Age + Condition_SSC_PAH, data=df)
summary(model)
summary(model)[4] #table
summary(model)[8] #Multiple R-squared
summary(model)[9] #adj.r.squared

#Age
#Age_Total
ggscatter(data, x = "Age", y = "1-methyl-4-imidazoleacetate",
          add = "reg.line",                                 # Add regression line
          conf.int = TRUE,                                  # Add confidence interval
          add.params = list(color = "blue",
                            fill = "lightgray"))+
  stat_cor(method = "spearman", label.x = 18, label.y = 4) 
#Age_By group
ggscatter(data, x = "Age", y = "1-methyl-4-imidazoleacetate",
          add = "reg.line",                         # Add regression line
          #conf.int = TRUE,                          # Add confidence interval
          color = "Condition", palette = c("lightcoral","darkgoldenrod2", "steelblue1"),                       # Color by groups "Condition"
          shape = "Condition"                             # Change point shape by groups "Condition"
)+
  stat_cor(aes(color = Condition)) 


#Condition
which(colnames(data)=="1-methyl-4-imidazoleacetate")
colnames(data)[which(colnames(data)=="1-methyl-4-imidazoleacetate")] = "imidazoleacetate"
colnames(data)[90]

my_comparisons <- list( c("HC", "SSC CON"), c("SSC CON", "SSC PAH"), c("HC", "SSC PAH") )

data %>% 
  ggplot(aes(x=Condition,y=imidazoleacetate,fill=Condition)) +
  geom_boxplot(lwd=0.5,fatten=1.05,outlier.shape = NA) + geom_jitter(width=0.11,alpha=0.3) + theme_bw() +
  ggtitle("1-methyl-4-imidazoleacetate") +  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 10)) + 
  xlab("Condition") + ylab("Relative abundance") + 
  labs(fill = "Condition") + scale_fill_manual(values=c("lightpink1", "lightgoldenrod", "lightskyblue")) +
  stat_compare_means(comparisons = my_comparisons, size=4)




# 2-aminoadipate
index <- which(colnames(data)=="2-aminoadipate")
index
df <- data.frame(
  metaboltite_name = data[,index],
  Age = data$Age,
  New_Gender = data$New_Gender,
  Condition_HC = data$Condition_HC,
  Condition_SSC_CON = data$Condition_SSC_CON,
  Condition_SSC_PAH = data$Condition_SSC_PAH
)
mean(df$metaboltite_name)    

model <- lm(metaboltite_name ~ Age + New_Gender + Condition_SSC_CON + Condition_SSC_PAH, data=df)
summary(model)
model <- lm(metaboltite_name ~ Age + Condition_SSC_CON + Condition_SSC_PAH, data=df)
summary(model)
summary(model)[8] 
summary(model)[9] 


#2-hydroxyglutarate
index <- which(colnames(data)=="2-hydroxyglutarate")
index
df <- data.frame(
  metaboltite_name = data[,index],
  Age = data$Age,
  New_Gender = data$New_Gender,
  Condition_HC = data$Condition_HC,
  Condition_SSC_CON = data$Condition_SSC_CON,
  Condition_SSC_PAH = data$Condition_SSC_PAH
)
mean(df$metaboltite_name)    

model <- lm(metaboltite_name ~ Age + New_Gender + Condition_SSC_CON + Condition_SSC_PAH, data=df)
summary(model)
model <- lm(metaboltite_name ~ Age + Condition_SSC_CON, data=df)
summary(model)
# Only Condition_SSC_CON is significant
summary(model)[8] 
summary(model)[9] 


#3beta,7alpha-dihydroxy-5-cholestenoate
index <- which(colnames(data)=="3beta,7alpha-dihydroxy-5-cholestenoate")
index
df <- data.frame(
  metaboltite_name = data[,index],
  Age = data$Age,
  New_Gender = data$New_Gender,
  Condition_HC = data$Condition_HC,
  Condition_SSC_CON = data$Condition_SSC_CON,
  Condition_SSC_PAH = data$Condition_SSC_PAH
)
mean(df$metaboltite_name)    

model <- lm(metaboltite_name ~ Age + New_Gender + Condition_SSC_CON + Condition_SSC_PAH, data=df)
summary(model)
model <- lm(metaboltite_name ~ New_Gender + Condition_SSC_CON + Condition_SSC_PAH, data=df)
summary(model)
model <- lm(metaboltite_name ~ Condition_SSC_PAH, data=df)
summary(model)


#4-acetamidobutanoate
index <- which(colnames(data)=="4-acetamidobutanoate")
index
df <- data.frame(
  metaboltite_name = data[,index],
  Age = data$Age,
  New_Gender = data$New_Gender,
  Condition_HC = data$Condition_HC,
  Condition_SSC_CON = data$Condition_SSC_CON,
  Condition_SSC_PAH = data$Condition_SSC_PAH
)
mean(df$metaboltite_name)    

model <- lm(metaboltite_name ~ Age + New_Gender + Condition_SSC_CON + Condition_SSC_PAH, data=df)
summary(model)
model <- lm(metaboltite_name ~ Age + Condition_SSC_PAH, data=df)
summary(model)
#Age is significant but has less effect than Condition

#Age
#Age_Total
ggscatter(data, x = "Age", y = "4-acetamidobutanoate",
          add = "reg.line",                                 # Add regression line
          conf.int = TRUE,                                  # Add confidence interval
          add.params = list(color = "blue",
                            fill = "lightgray"))+
  stat_cor(method = "spearman", label.x = 18, label.y = 4) 
#Age_By group
ggscatter(data, x = "Age", y = "4-acetamidobutanoate",
          add = "reg.line",                         # Add regression line
          #conf.int = TRUE,                          # Add confidence interval
          color = "Condition", palette = c("lightcoral","darkgoldenrod2", "steelblue1"),                       # Color by groups "Condition"
          shape = "Condition"                             # Change point shape by groups "Condition"
)+
  stat_cor(aes(color = Condition)) 


#Condition
which(colnames(data)=="4-acetamidobutanoate")
colnames(data)[which(colnames(data)=="4-acetamidobutanoate")] = "acetamidobutanoate"
colnames(data)[95]

my_comparisons <- list( c("HC", "SSC CON"), c("SSC CON", "SSC PAH"), c("HC", "SSC PAH") )

data %>% 
  ggplot(aes(x=Condition,y=acetamidobutanoate,fill=Condition)) +
  geom_boxplot(lwd=0.5,fatten=1.05,outlier.shape = NA) + geom_jitter(width=0.11,alpha=0.3) + theme_bw() +
  ggtitle("4-acetamidobutanoate") +  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 10)) + 
  xlab("Condition") + ylab("Relative abundance") + 
  labs(fill = "Condition") + scale_fill_manual(values=c("lightpink1", "lightgoldenrod", "lightskyblue")) +
  stat_compare_means(comparisons = my_comparisons, size=4)





#5,6-dihydrouracil
index <- which(colnames(data)=="5,6-dihydrouracil")
index
df <- data.frame(
  metaboltite_name = data[,index],
  Age = data$Age,
  New_Gender = data$New_Gender,
  Condition_HC = data$Condition_HC,
  Condition_SSC_CON = data$Condition_SSC_CON,
  Condition_SSC_PAH = data$Condition_SSC_PAH
)
mean(df$metaboltite_name)    

model <- lm(metaboltite_name ~ Age + New_Gender + Condition_SSC_CON + Condition_SSC_PAH, data=df)
summary(model)
model <- lm(metaboltite_name ~ Age + New_Gender + Condition_SSC_PAH, data=df)
summary(model)
model <- lm(metaboltite_name ~ Condition_SSC_PAH, data=df)
summary(model)
# Age and New_Gender are not significant

#7-methylurate
index <- which(colnames(data)=="7-methylurate")
index
df <- data.frame(
  metaboltite_name = data[,index],
  Age = data$Age,
  New_Gender = data$New_Gender,
  Condition_HC = data$Condition_HC,
  Condition_SSC_CON = data$Condition_SSC_CON,
  Condition_SSC_PAH = data$Condition_SSC_PAH
)
mean(df$metaboltite_name) 
model <- lm(metaboltite_name ~ Age + New_Gender + Condition_SSC_CON + Condition_SSC_PAH, data=df)
summary(model)
model <- lm(metaboltite_name ~ Condition_SSC_PAH, data=df)
summary(model)


#7-methylxanthine
index <- which(colnames(data)=="7-methylxanthine")
index
df <- data.frame(
  metaboltite_name = data[,index],
  Age = data$Age,
  New_Gender = data$New_Gender,
  Condition_HC = data$Condition_HC,
  Condition_SSC_CON = data$Condition_SSC_CON,
  Condition_SSC_PAH = data$Condition_SSC_PAH
)
mean(df$metaboltite_name) 
model <- lm(metaboltite_name ~ Age + New_Gender + Condition_SSC_CON + Condition_SSC_PAH, data=df)
summary(model)
model <- lm(metaboltite_name ~ Condition_SSC_PAH, data=df)
summary(model)


#adenosine 5'-monophosphate (AMP)
index <- which(colnames(data)=="adenosine 5'-monophosphate (AMP)")
index
df <- data.frame(
  metaboltite_name = data[,index],
  Age = data$Age,
  New_Gender = data$New_Gender,
  Condition_HC = data$Condition_HC,
  Condition_SSC_CON = data$Condition_SSC_CON,
  Condition_SSC_PAH = data$Condition_SSC_PAH
)
mean(df$metaboltite_name) 
model <- lm(metaboltite_name ~ Age + New_Gender + Condition_SSC_CON + Condition_SSC_PAH, data=df)
summary(model)
model <- lm(metaboltite_name ~ Condition_SSC_PAH, data=df)
summary(model)


#arabitol/xylitol
index <- which(colnames(data)=="arabitol/xylitol")
index
df <- data.frame(
  metaboltite_name = data[,index],
  Age = data$Age,
  New_Gender = data$New_Gender,
  Condition_HC = data$Condition_HC,
  Condition_SSC_CON = data$Condition_SSC_CON,
  Condition_SSC_PAH = data$Condition_SSC_PAH
)
mean(df$metaboltite_name) 
model <- lm(metaboltite_name ~ Age + New_Gender + Condition_SSC_CON + Condition_SSC_PAH, data=df)
summary(model)
model <- lm(metaboltite_name ~ Age + New_Gender + Condition_SSC_PAH, data=df)
summary(model)



#arginine
index <- which(colnames(data)=="arginine")
index
df <- data.frame(
  metaboltite_name = data[,index],
  Age = data$Age,
  New_Gender = data$New_Gender,
  Condition_HC = data$Condition_HC,
  Condition_SSC_CON = data$Condition_SSC_CON,
  Condition_SSC_PAH = data$Condition_SSC_PAH
)
mean(df$metaboltite_name) 
model <- lm(metaboltite_name ~ Age + New_Gender + Condition_SSC_CON + Condition_SSC_PAH, data=df)
summary(model)
model <- lm(metaboltite_name ~ Age + Condition_SSC_CON + Condition_SSC_PAH, data=df)
summary(model)


#carnitine
index <- which(colnames(data)=="carnitine")
index
df <- data.frame(
  metaboltite_name = data[,index],
  Age = data$Age,
  New_Gender = data$New_Gender,
  Condition_HC = data$Condition_HC,
  Condition_SSC_CON = data$Condition_SSC_CON,
  Condition_SSC_PAH = data$Condition_SSC_PAH
)
mean(df$metaboltite_name) 
model <- lm(metaboltite_name ~ Age + New_Gender + Condition_SSC_CON + Condition_SSC_PAH, data=df)
summary(model)
model <- lm(metaboltite_name ~ Age + Condition_SSC_CON + Condition_SSC_PAH, data=df)
summary(model)
model <- lm(metaboltite_name ~ Age + Condition_SSC_CON, data=df)
summary(model)
model <- lm(metaboltite_name ~ Condition_SSC_CON, data=df)
summary(model)


#dimethylglycine
index <- which(colnames(data)=="dimethylglycine")
index
df <- data.frame(
  metaboltite_name = data[,index],
  Age = data$Age,
  New_Gender = data$New_Gender,
  Condition_HC = data$Condition_HC,
  Condition_SSC_CON = data$Condition_SSC_CON,
  Condition_SSC_PAH = data$Condition_SSC_PAH
)
mean(df$metaboltite_name) 
model <- lm(metaboltite_name ~ Age + New_Gender + Condition_SSC_CON + Condition_SSC_PAH, data=df)
summary(model)
model <- lm(metaboltite_name ~ Age + Condition_SSC_CON + Condition_SSC_PAH, data=df)
summary(model)
model <- lm(metaboltite_name ~ Age + Condition_SSC_CON, data=df)
summary(model)



#fumarate
index <- which(colnames(data)=="fumarate")
index
df <- data.frame(
  metaboltite_name = data[,index],
  Age = data$Age,
  New_Gender = data$New_Gender,
  Condition_HC = data$Condition_HC,
  Condition_SSC_CON = data$Condition_SSC_CON,
  Condition_SSC_PAH = data$Condition_SSC_PAH
)
mean(df$metaboltite_name) 
model <- lm(metaboltite_name ~ Age + New_Gender + Condition_SSC_CON + Condition_SSC_PAH, data=df)
summary(model)
model <- lm(metaboltite_name ~ Condition_SSC_PAH, data=df)
summary(model)




#glycine
index <- which(colnames(data)=="glycine")
index
df <- data.frame(
  metaboltite_name = data[,index],
  Age = data$Age,
  New_Gender = data$New_Gender,
  Condition_HC = data$Condition_HC,
  Condition_SSC_CON = data$Condition_SSC_CON,
  Condition_SSC_PAH = data$Condition_SSC_PAH
)
mean(df$metaboltite_name) 
model <- lm(metaboltite_name ~ Age + New_Gender + Condition_SSC_CON + Condition_SSC_PAH, data=df)
summary(model)
model <- lm(metaboltite_name ~ New_Gender + Condition_SSC_CON + Condition_SSC_PAH, data=df)
summary(model)
model <- lm(metaboltite_name ~ New_Gender + Condition_SSC_PAH, data=df)
summary(model)
model <- lm(metaboltite_name ~ Condition_SSC_PAH, data=df)
summary(model)





#	gulonate*
index <- which(colnames(data)=="gulonate*")
index
df <- data.frame(
  metaboltite_name = data[,index],
  Age = data$Age,
  New_Gender = data$New_Gender,
  Condition_HC = data$Condition_HC,
  Condition_SSC_CON = data$Condition_SSC_CON,
  Condition_SSC_PAH = data$Condition_SSC_PAH
)
mean(df$metaboltite_name) 
model <- lm(metaboltite_name ~ Age + New_Gender + Condition_SSC_CON + Condition_SSC_PAH, data=df)
summary(model)
model <- lm(metaboltite_name ~ Age + New_Gender + Condition_SSC_PAH, data=df)
summary(model)
model <- lm(metaboltite_name ~ Age + Condition_SSC_PAH, data=df)
summary(model)




#histidine
index <- which(colnames(data)=="histidine")
index
df <- data.frame(
  metaboltite_name = data[,index],
  Age = data$Age,
  New_Gender = data$New_Gender,
  Condition_HC = data$Condition_HC,
  Condition_SSC_CON = data$Condition_SSC_CON,
  Condition_SSC_PAH = data$Condition_SSC_PAH
)
mean(df$metaboltite_name) 
model <- lm(metaboltite_name ~ Age + New_Gender + Condition_SSC_CON + Condition_SSC_PAH, data=df)
summary(model)
model <- lm(metaboltite_name ~ New_Gender + Condition_SSC_CON + Condition_SSC_PAH, data=df)
summary(model)
model <- lm(metaboltite_name ~ Condition_SSC_CON + Condition_SSC_PAH, data=df)
summary(model)
model <- lm(metaboltite_name ~ Condition_SSC_PAH, data=df)
summary(model)



#kynurenine
index <- which(colnames(data)=="kynurenine")
index
df <- data.frame(
  metaboltite_name = data[,index],
  Age = data$Age,
  New_Gender = data$New_Gender,
  Condition_HC = data$Condition_HC,
  Condition_SSC_CON = data$Condition_SSC_CON,
  Condition_SSC_PAH = data$Condition_SSC_PAH
)
mean(df$metaboltite_name) 
model <- lm(metaboltite_name ~ Age + New_Gender + Condition_SSC_CON + Condition_SSC_PAH, data=df)
summary(model)
model <- lm(metaboltite_name ~ Age + New_Gender + Condition_SSC_PAH, data=df)
summary(model)
model <- lm(metaboltite_name ~ Condition_SSC_PAH, data=df)
summary(model)


#lactate
index <- which(colnames(data)=="lactate")
index
df <- data.frame(
  metaboltite_name = data[,index],
  Age = data$Age,
  New_Gender = data$New_Gender,
  Condition_HC = data$Condition_HC,
  Condition_SSC_CON = data$Condition_SSC_CON,
  Condition_SSC_PAH = data$Condition_SSC_PAH
)
mean(df$metaboltite_name) 
model <- lm(metaboltite_name ~ Age + New_Gender + Condition_SSC_CON + Condition_SSC_PAH, data=df)
summary(model)
model <- lm(metaboltite_name ~ New_Gender + Condition_SSC_CON + Condition_SSC_PAH, data=df)
summary(model)
model <- lm(metaboltite_name ~ Condition_SSC_CON + Condition_SSC_PAH, data=df)
summary(model)
model <- lm(metaboltite_name ~ Condition_SSC_PAH, data=df)
summary(model)


#leucine
index <- which(colnames(data)=="leucine")
index
df <- data.frame(
  metaboltite_name = data[,index],
  Age = data$Age,
  New_Gender = data$New_Gender,
  Condition_HC = data$Condition_HC,
  Condition_SSC_CON = data$Condition_SSC_CON,
  Condition_SSC_PAH = data$Condition_SSC_PAH
)
mean(df$metaboltite_name) 
model <- lm(metaboltite_name ~ Age + New_Gender + Condition_SSC_CON + Condition_SSC_PAH, data=df)
summary(model)
model <- lm(metaboltite_name ~ Age + New_Gender + Condition_SSC_PAH, data=df)
summary(model)
model <- lm(metaboltite_name ~ Age + New_Gender, data=df)
summary(model)
model <- lm(metaboltite_name ~ New_Gender, data=df)
summary(model)


#N1-Methyl-2-pyridone-5-carboxamide
index <- which(colnames(data)=="N1-Methyl-2-pyridone-5-carboxamide")
index
df <- data.frame(
  metaboltite_name = data[,index],
  Age = data$Age,
  New_Gender = data$New_Gender,
  Condition_HC = data$Condition_HC,
  Condition_SSC_CON = data$Condition_SSC_CON,
  Condition_SSC_PAH = data$Condition_SSC_PAH
)
mean(df$metaboltite_name) 
model <- lm(metaboltite_name ~ Age + New_Gender + Condition_SSC_CON + Condition_SSC_PAH, data=df)
summary(model)
model <- lm(metaboltite_name ~ Age + New_Gender + Condition_SSC_PAH, data=df)
summary(model)
model <- lm(metaboltite_name ~ Age + Condition_SSC_PAH, data=df)
summary(model)
model <- lm(metaboltite_name ~ Condition_SSC_PAH, data=df)
summary(model)



#N-acetylneuraminate
index <- which(colnames(data)=="N-acetylneuraminate")
index
df <- data.frame(
  metaboltite_name = data[,index],
  Age = data$Age,
  New_Gender = data$New_Gender,
  Condition_HC = data$Condition_HC,
  Condition_SSC_CON = data$Condition_SSC_CON,
  Condition_SSC_PAH = data$Condition_SSC_PAH
)
mean(df$metaboltite_name) 
model <- lm(metaboltite_name ~ Age + New_Gender + Condition_SSC_CON + Condition_SSC_PAH, data=df)
summary(model)
model <- lm(metaboltite_name ~ Age + Condition_SSC_CON + Condition_SSC_PAH, data=df)
summary(model)
model <- lm(metaboltite_name ~ Age + Condition_SSC_PAH, data=df)
summary(model)

#Age
#Age_Total
ggscatter(data, x = "Age", y = "N-acetylneuraminate",
          add = "reg.line",                                 # Add regression line
          conf.int = TRUE,                                  # Add confidence interval
          add.params = list(color = "blue",
                            fill = "lightgray"))+
  stat_cor(method = "spearman", label.x = 18, label.y = 4) 
#Age_By group
ggscatter(data, x = "Age", y = "N-acetylneuraminate",
          add = "reg.line",                         # Add regression line
          #conf.int = TRUE,                          # Add confidence interval
          color = "Condition", palette = c("lightcoral","darkgoldenrod2", "steelblue1"),                       # Color by groups "Condition"
          shape = "Condition"                             # Change point shape by groups "Condition"
)+
  stat_cor(aes(color = Condition)) 


#Condition
which(colnames(data)=="N-acetylneuraminate")
colnames(data)[which(colnames(data)=="N-acetylneuraminate")] = "acetylneuraminate"
colnames(data)[32]

my_comparisons <- list( c("HC", "SSC CON"), c("SSC CON", "SSC PAH"), c("HC", "SSC PAH") )

data %>% 
  ggplot(aes(x=Condition,y=acetylneuraminate,fill=Condition)) +
  geom_boxplot(lwd=0.5,fatten=1.05,outlier.shape = NA) + geom_jitter(width=0.11,alpha=0.3) + theme_bw() +
  ggtitle("N-acetylneuraminate") +  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 10)) + 
  xlab("Condition") + ylab("Relative abundance") + 
  labs(fill = "Condition") + scale_fill_manual(values=c("lightpink1", "lightgoldenrod", "lightskyblue")) +
  stat_compare_means(comparisons = my_comparisons, size=4)



#N-carbamoylaspartate
index <- which(colnames(data)=="N-carbamoylaspartate")
index
df <- data.frame(
  metaboltite_name = data[,index],
  Age = data$Age,
  New_Gender = data$New_Gender,
  Condition_HC = data$Condition_HC,
  Condition_SSC_CON = data$Condition_SSC_CON,
  Condition_SSC_PAH = data$Condition_SSC_PAH
)
mean(df$metaboltite_name) 
model <- lm(metaboltite_name ~ Age + New_Gender + Condition_SSC_CON + Condition_SSC_PAH, data=df)
summary(model)
model <- lm(metaboltite_name ~ Age + Condition_SSC_CON + Condition_SSC_PAH, data=df)
summary(model)
model <- lm(metaboltite_name ~ Age + Condition_SSC_PAH, data=df)
summary(model)
model <- lm(metaboltite_name ~ Condition_SSC_PAH, data=df)
summary(model)






#nicotinamide
index <- which(colnames(data)=="nicotinamide")
index
df <- data.frame(
  metaboltite_name = data[,index],
  Age = data$Age,
  New_Gender = data$New_Gender,
  Condition_HC = data$Condition_HC,
  Condition_SSC_CON = data$Condition_SSC_CON,
  Condition_SSC_PAH = data$Condition_SSC_PAH
)
mean(df$metaboltite_name) 
model <- lm(metaboltite_name ~ Age + New_Gender + Condition_SSC_CON + Condition_SSC_PAH, data=df)
summary(model)
model <- lm(metaboltite_name ~ Age + New_Gender + Condition_SSC_PAH, data=df)
summary(model)
model <- lm(metaboltite_name ~ Age + Condition_SSC_PAH, data=df)
summary(model)
model <- lm(metaboltite_name ~ Condition_SSC_PAH, data=df)
summary(model)




#phenylacetate
index <- which(colnames(data)=="phenylacetate")
index
df <- data.frame(
  metaboltite_name = data[,index],
  Age = data$Age,
  New_Gender = data$New_Gender,
  Condition_HC = data$Condition_HC,
  Condition_SSC_CON = data$Condition_SSC_CON,
  Condition_SSC_PAH = data$Condition_SSC_PAH
)
mean(df$metaboltite_name) 
model <- lm(metaboltite_name ~ Age + New_Gender + Condition_SSC_CON + Condition_SSC_PAH, data=df)
summary(model)
model <- lm(metaboltite_name ~ Age + Condition_SSC_CON + Condition_SSC_PAH, data=df)
summary(model)
model <- lm(metaboltite_name ~ Age + Condition_SSC_PAH, data=df)
summary(model)




# pyridoxate
index <- which(colnames(data)=="pyridoxate")
index
df <- data.frame(
  metaboltite_name = data[,index],
  Age = data$Age,
  New_Gender = data$New_Gender,
  Condition_HC = data$Condition_HC,
  Condition_SSC_CON = data$Condition_SSC_CON,
  Condition_SSC_PAH = data$Condition_SSC_PAH
)
mean(df$metaboltite_name) 
model <- lm(metaboltite_name ~ Age + New_Gender + Condition_SSC_CON + Condition_SSC_PAH, data=df)
summary(model)
model <- lm(metaboltite_name ~ New_Gender + Condition_SSC_CON + Condition_SSC_PAH, data=df)
summary(model)
model <- lm(metaboltite_name ~ Condition_SSC_CON + Condition_SSC_PAH, data=df)
summary(model)
model <- lm(metaboltite_name ~ Condition_SSC_PAH, data=df)
summary(model)






#quinolinate
index <- which(colnames(data)=="quinolinate")
index
df <- data.frame(
  metaboltite_name = data[,index],
  Age = data$Age,
  New_Gender = data$New_Gender,
  Condition_HC = data$Condition_HC,
  Condition_SSC_CON = data$Condition_SSC_CON,
  Condition_SSC_PAH = data$Condition_SSC_PAH
)
mean(df$metaboltite_name) 
model <- lm(metaboltite_name ~ Age + New_Gender + Condition_SSC_CON + Condition_SSC_PAH, data=df)
summary(model)
model <- lm(metaboltite_name ~ Age + New_Gender + Condition_SSC_PAH, data=df)
summary(model)
model <- lm(metaboltite_name ~ Age + Condition_SSC_PAH, data=df)
summary(model)
model <- lm(metaboltite_name ~ Condition_SSC_PAH, data=df)
summary(model)




#retinal
index <- which(colnames(data)=="retinal")
index
df <- data.frame(
  metaboltite_name = data[,index],
  Age = data$Age,
  New_Gender = data$New_Gender,
  Condition_HC = data$Condition_HC,
  Condition_SSC_CON = data$Condition_SSC_CON,
  Condition_SSC_PAH = data$Condition_SSC_PAH
)
mean(df$metaboltite_name) 
model <- lm(metaboltite_name ~ Age + New_Gender + Condition_SSC_CON + Condition_SSC_PAH, data=df)
summary(model)
model <- lm(metaboltite_name ~ Age + Condition_SSC_CON + Condition_SSC_PAH, data=df)
summary(model)
model <- lm(metaboltite_name ~ Age + Condition_SSC_PAH, data=df)
summary(model)
model <- lm(metaboltite_name ~ Condition_SSC_PAH, data=df)
summary(model)


#serine
index <- which(colnames(data)=="serine")
index
df <- data.frame(
  metaboltite_name = data[,index],
  Age = data$Age,
  New_Gender = data$New_Gender,
  Condition_HC = data$Condition_HC,
  Condition_SSC_CON = data$Condition_SSC_CON,
  Condition_SSC_PAH = data$Condition_SSC_PAH
)
mean(df$metaboltite_name) 
model <- lm(metaboltite_name ~ Age + New_Gender + Condition_SSC_CON + Condition_SSC_PAH, data=df)
summary(model)
model <- lm(metaboltite_name ~ Age + New_Gender + Condition_SSC_PAH, data=df)
summary(model)
model <- lm(metaboltite_name ~ New_Gender + Condition_SSC_PAH, data=df)
summary(model)
model <- lm(metaboltite_name ~ Condition_SSC_PAH, data=df)
summary(model)







#sphinganine
index <- which(colnames(data)=="sphinganine")
index
df <- data.frame(
  metaboltite_name = data[,index],
  Age = data$Age,
  New_Gender = data$New_Gender,
  Condition_HC = data$Condition_HC,
  Condition_SSC_CON = data$Condition_SSC_CON,
  Condition_SSC_PAH = data$Condition_SSC_PAH
)
mean(df$metaboltite_name) 
model <- lm(metaboltite_name ~ Age + New_Gender + Condition_SSC_CON + Condition_SSC_PAH, data=df)
summary(model)
model <- lm(metaboltite_name ~ New_Gender + Condition_SSC_CON + Condition_SSC_PAH, data=df)
summary(model)
model <- lm(metaboltite_name ~ Condition_SSC_CON + Condition_SSC_PAH, data=df)
summary(model)


# sphingomyelin (d18:2/18:1)*
index <- which(colnames(data)=="sphingomyelin (d18:2/18:1)*")
index
df <- data.frame(
  metaboltite_name = data[,index],
  Age = data$Age,
  New_Gender = data$New_Gender,
  Condition_HC = data$Condition_HC,
  Condition_SSC_CON = data$Condition_SSC_CON,
  Condition_SSC_PAH = data$Condition_SSC_PAH
)
mean(df$metaboltite_name) 
model <- lm(metaboltite_name ~ Age + New_Gender + Condition_SSC_CON + Condition_SSC_PAH, data=df)
summary(model)
model <- lm(metaboltite_name ~ Age + New_Gender + Condition_SSC_CON, data=df)
summary(model)
model <- lm(metaboltite_name ~ New_Gender + Condition_SSC_CON, data=df)
summary(model)
model <- lm(metaboltite_name ~ New_Gender, data=df)
summary(model)








#sphingosine
index <- which(colnames(data)=="sphingosine")
index
df <- data.frame(
  metaboltite_name = data[,index],
  Age = data$Age,
  New_Gender = data$New_Gender,
  Condition_HC = data$Condition_HC,
  Condition_SSC_CON = data$Condition_SSC_CON,
  Condition_SSC_PAH = data$Condition_SSC_PAH
)
mean(df$metaboltite_name) 
model <- lm(metaboltite_name ~ Age + New_Gender + Condition_SSC_CON + Condition_SSC_PAH, data=df)
summary(model)
model <- lm(metaboltite_name ~ New_Gender + Condition_SSC_CON + Condition_SSC_PAH, data=df)
summary(model)
model <- lm(metaboltite_name ~ Condition_SSC_CON + Condition_SSC_PAH, data=df)
summary(model)
model <- lm(metaboltite_name ~ Condition_SSC_PAH, data=df)
summary(model)







#sucrose
index <- which(colnames(data)=="sucrose")
index
df <- data.frame(
  metaboltite_name = data[,index],
  Age = data$Age,
  New_Gender = data$New_Gender,
  Condition_HC = data$Condition_HC,
  Condition_SSC_CON = data$Condition_SSC_CON,
  Condition_SSC_PAH = data$Condition_SSC_PAH
)
mean(df$metaboltite_name) 
model <- lm(metaboltite_name ~ Age + New_Gender + Condition_SSC_CON + Condition_SSC_PAH, data=df)
summary(model)
model <- lm(metaboltite_name ~ Age + New_Gender + Condition_SSC_PAH, data=df)
summary(model)
model <- lm(metaboltite_name ~ New_Gender + Condition_SSC_PAH, data=df)
summary(model)
model <- lm(metaboltite_name ~ Condition_SSC_PAH, data=df)
summary(model)








#uracil
index <- which(colnames(data)=="uracil")
index
df <- data.frame(
  metaboltite_name = data[,index],
  Age = data$Age,
  New_Gender = data$New_Gender,
  Condition_HC = data$Condition_HC,
  Condition_SSC_CON = data$Condition_SSC_CON,
  Condition_SSC_PAH = data$Condition_SSC_PAH
)
mean(df$metaboltite_name) 
model <- lm(metaboltite_name ~ Age + New_Gender + Condition_SSC_CON + Condition_SSC_PAH, data=df)
summary(model)
model <- lm(metaboltite_name ~ New_Gender + Condition_SSC_CON + Condition_SSC_PAH, data=df)
summary(model)
model <- lm(metaboltite_name ~ Condition_SSC_CON + Condition_SSC_PAH, data=df)
summary(model)
model <- lm(metaboltite_name ~ Condition_SSC_PAH, data=df)
summary(model)








#vanillylmandelate (VMA)
index <- which(colnames(data)=="vanillylmandelate (VMA)")
index
df <- data.frame(
  metaboltite_name = data[,index],
  Age = data$Age,
  New_Gender = data$New_Gender,
  Condition_HC = data$Condition_HC,
  Condition_SSC_CON = data$Condition_SSC_CON,
  Condition_SSC_PAH = data$Condition_SSC_PAH
)
mean(df$metaboltite_name) 
model <- lm(metaboltite_name ~ Age + New_Gender + Condition_SSC_CON + Condition_SSC_PAH, data=df)
summary(model)
model <- lm(metaboltite_name ~ Age + New_Gender + Condition_SSC_PAH, data=df)
summary(model)
model <- lm(metaboltite_name ~ Age + Condition_SSC_PAH, data=df)
summary(model)


#Age
#Age_Total
ggscatter(data, x = "Age", y = "vanillylmandelate (VMA)",
          add = "reg.line",                                 # Add regression line
          conf.int = TRUE,                                  # Add confidence interval
          add.params = list(color = "blue",
                            fill = "lightgray"))+
  stat_cor(method = "spearman", label.x = 18, label.y = 4) 
#Age_By group
ggscatter(data, x = "Age", y = "vanillylmandelate (VMA)",
          add = "reg.line",                         # Add regression line
          #conf.int = TRUE,                          # Add confidence interval
          color = "Condition", palette = c("lightcoral","darkgoldenrod2", "steelblue1"),                       # Color by groups "Condition"
          shape = "Condition"                             # Change point shape by groups "Condition"
)+
  stat_cor(aes(color = Condition)) 


#Condition
which(colnames(data)=="vanillylmandelate (VMA)")
colnames(data)[which(colnames(data)=="vanillylmandelate (VMA)")] = "vanillylmandelate"
colnames(data)[39]

my_comparisons <- list( c("HC", "SSC CON"), c("SSC CON", "SSC PAH"), c("HC", "SSC PAH") )

data %>% 
  ggplot(aes(x=Condition,y=vanillylmandelate,fill=Condition)) +
  geom_boxplot(lwd=0.5,fatten=1.05,outlier.shape = NA) + geom_jitter(width=0.11,alpha=0.3) + theme_bw() +
  ggtitle("Vanillylmandelate (VMA)") +  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 10)) + 
  xlab("Condition") + ylab("Relative abundance") + 
  labs(fill = "Condition") + scale_fill_manual(values=c("lightpink1", "lightgoldenrod", "lightskyblue")) +
  stat_compare_means(comparisons = my_comparisons, size=4)

























#Scatter Plot


#Age_Total
ggscatter(data, x = "Age", y = "leucine",
          add = "reg.line",                                 # Add regression line
          conf.int = TRUE,                                  # Add confidence interval
          add.params = list(color = "blue",
                            fill = "lightgray"))+
          stat_cor(method = "spearman", label.x = 18, label.y = 4) 
#Age_By group
ggscatter(data, x = "Age", y = "leucine",
          add = "reg.line",                         # Add regression line
          #conf.int = TRUE,                          # Add confidence interval
          color = "Condition", palette = c("lightcoral", "darkgoldenrod2", "steelblue1"),                       # Color by groups "Condition"
          shape = "Condition"                             # Change point shape by groups "Condition"
)+
  stat_cor(aes(color = Condition), label.x = 3)   

#Gender
compare_means(leucine ~ Gender, method = "t.test",  data = data)

data %>% 
  ggplot(aes(x=Gender, y=leucine, fill=Gender)) +
  geom_boxplot(lwd=0.5,fatten=1.05,outlier.shape = NA) + geom_jitter(width=0.11,alpha=0.3) + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 10)) + 
  xlab("Group") + ylab("Relative Abundance") + 
  labs(fill = "Group") + scale_fill_manual(values=c("lightpink1", "lightskyblue")) +
  stat_compare_means(method = "t.test", size=4)

varImp(model)



