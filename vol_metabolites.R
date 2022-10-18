
setwd('D:/Research/3. Volcanol and Pathway analysis')
metabolites <- read.csv('D:/Research/2. Raw_data/Batch_norm_name.csv',check.names=FALSE)
#metabolites <- metabolites[-c(61:90),] #delete PF group
dim(metabolites) # 120x1191

mean(metabolites[1:30,3]) #PAH first metabolite
mean(metabolites[31:60, 3]) #HC first metabolite
mean(metabolites[61:90, 3]) #PF first metabolite
mean(metabolites[91:120, 3]) #SSC first metabolite


table <- data.frame(metabolites=colnames(metabolites)[3:1191])

#table$LC_PAH_HC[1] <- mean(metabolites[1:30,3])/mean(metabolites[31:60, 3])
#table$LC_PAH_HC[2] <- mean(metabolites[1:30,4])/mean(metabolites[31:60, 4])

for (i in 1:1189){
  table$LC_PAH_HC[i] <- mean(metabolites[1:30,i+2])/mean(metabolites[31:60, i+2]) #PAH vs HC
  table$LC_SSC_HC[i] <- mean(metabolites[91:120,i+2])/mean(metabolites[31:60, i+2]) #SSC vs HC
  table$LC_PAH_SSC[i] <- mean(metabolites[1:30,i+2])/mean(metabolites[91:120, i+2]) #PAH vs SSC
  table$LC_PF_HC[i] <- mean(metabolites[61:90,i+2])/mean(metabolites[31:60, i+2])#PF vs HC
  table$LC_PF_SSC[i] <- mean(metabolites[61:90,i+2])/mean(metabolites[91:120, i+2])#PF vs HC
  table$LC_PAH_PF[i] <- mean(metabolites[1:30,i+2])/mean(metabolites[61:90, i+2])#PAH vs PF
  
  table$log2FC_PAH_HC[i] <- log2(mean(metabolites[1:30,i+2])/mean(metabolites[31:60, i+2]))
  table$log2FC_SSC_HC[i] <- log2(mean(metabolites[91:120,i+2])/mean(metabolites[31:60, i+2]))
  table$log2FC_PAH_SSC[i] <- log2(mean(metabolites[1:30,i+2])/mean(metabolites[91:120, i+2]))
  table$log2FC_PF_HC[i] <- log2(mean(metabolites[61:90,i+2])/mean(metabolites[31:60, i+2]))
  table$log2FC_PF_SSC[i] <- log2(mean(metabolites[61:90,i+2])/mean(metabolites[91:120, i+2]))
  table$log2FC_PAH_PF[i] <- log2(mean(metabolites[1:30,i+2])/mean(metabolites[61:90, i+2]))
  
}


library(matrixTests)

#table$pvalue_PAH_HC[1] <- row_t_welch(metabolites[1:30,3], metabolites[31:60,3]) #跳過無法產生p value的值

for (i in 1:1189){
  
  table$pvalue_PAH_HC[i] <- row_t_welch(metabolites[1:30,i+2], metabolites[31:60,i+2])[12]
  table$pvalue_SSC_HC[i] <- row_t_welch(metabolites[91:120,i+2], metabolites[31:60,i+2])[12]
  table$pvalue_PAH_SSC[i] <- row_t_welch(metabolites[1:30,i+2], metabolites[91:120,i+2])[12]
  table$pvalue_PF_HC[i] <- row_t_welch(metabolites[61:90,i+2], metabolites[31:60,i+2])[12]
  table$pvalue_PF_SSC[i] <- row_t_welch(metabolites[61:90,i+2], metabolites[91:120,i+2])[12]
  table$pvalue_PAH_PF[i] <- row_t_welch(metabolites[1:30,i+2], metabolites[61:90,i+2])[12]
  
  table$log10_pvalue_PAH_HC[i] <- log10(row_t_welch(metabolites[1:30,i+2], metabolites[31:60,i+2])[12])
  table$log10_pvalue_SSC_HC[i] <- log10(row_t_welch(metabolites[91:120,i+2], metabolites[31:60,i+2])[12])
  table$log10_pvalue_PAH_SSC[i] <- log10(row_t_welch(metabolites[1:30,i+2], metabolites[91:120,i+2])[12])
  table$log10_pvalue_PF_HC[i] <- log10(row_t_welch(metabolites[61:90,i+2], metabolites[31:60,i+2])[12])
  table$log10_pvalue_PF_SSC[i] <- log10(row_t_welch(metabolites[61:90,i+2], metabolites[91:120,i+2])[12])
  table$log10_pvalue_PAH_PF[i] <- log10(row_t_welch(metabolites[1:30,i+2], metabolites[61:90,i+2])[12])

}

dim(metabolites)
dim(table)

library(dplyr)

table <- table %>%
  mutate(PAH_vs_HC_threshold = factor(case_when(log2FC_PAH_HC > 0 & pvalue_PAH_HC < 0.01 ~ "UP",
                                      log2FC_PAH_HC < 0 & pvalue_PAH_HC < 0.01 ~ "DOWN",
                                      TRUE ~ "NOT SIGNIFICANTLY DIFFERENT")))
table <- table %>%
  mutate(SSC_vs_HC_threshold = factor(case_when(log2FC_SSC_HC > 0 & pvalue_SSC_HC < 0.01 ~ "UP",
                                                log2FC_SSC_HC < 0 & pvalue_SSC_HC < 0.01 ~ "DOWN",
                                                TRUE ~ "NOT SIGNIFICANTLY DIFFERENT")))
table <- table %>%
  mutate(PAH_vs_SSC_threshold = factor(case_when(log2FC_PAH_SSC > 0 & pvalue_PAH_SSC < 0.01 ~ "UP",
                                                 log2FC_PAH_SSC < 0 & pvalue_PAH_SSC < 0.01 ~ "DOWN",
                                                TRUE ~ "NOT SIGNIFICANTLY DIFFERENT")))
table <- table %>%
  mutate(PF_vs_HC_threshold = factor(case_when(log2FC_PF_HC > 0 & pvalue_PF_HC < 0.01 ~ "UP",
                                                log2FC_PF_HC < 0 & pvalue_PF_HC < 0.01 ~ "DOWN",
                                                TRUE ~ "NOT SIGNIFICANTLY DIFFERENT")))
table <- table %>%
  mutate(PF_vs_SSC_threshold = factor(case_when(log2FC_PF_SSC > 0 & pvalue_PF_SSC < 0.01 ~ "UP",
                                               log2FC_PF_SSC < 0 & pvalue_PF_SSC < 0.01 ~ "DOWN",
                                               TRUE ~ "NOT SIGNIFICANTLY DIFFERENT")))
table <- table %>%
  mutate(PAH_vs_PF_threshold = factor(case_when(log2FC_PAH_PF > 0 & pvalue_PAH_PF < 0.01 ~ "UP",
                                                log2FC_PAH_PF < 0 & pvalue_PAH_PF < 0.01 ~ "DOWN",
                                                TRUE ~ "NOT SIGNIFICANTLY DIFFERENT")))

summary(table$PAH_vs_HC_threshold)
summary(table$SSC_vs_HC_threshold)
summary(table$PAH_vs_SSC_threshold)
summary(table$PF_vs_HC_threshold)
summary(table$PF_vs_SSC_threshold)
summary(table$PAH_vs_PF_threshold)


#Volcanol
library("BiocManager")
library("ggrepel")

table <- as.data.frame(lapply(table, unlist))

#PAH vs HC
table <- table[order(table$pvalue_PAH_HC),]

ggplot(table, aes(x=log2FC_PAH_HC, y=-log10(pvalue_PAH_HC))) +
  ggtitle("SSc-PAH vs HC") +
  xlab("log2 fold change") + 
  ylab("-log10 p-value") +
  geom_point(aes(colour = PAH_vs_HC_threshold),size=1) +
  scale_color_manual(values=c("steelblue2", "gray55","indianred1")) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) +
  geom_text_repel(data=head(table, 20), aes(label=metabolites), max.overlaps = Inf, size=4) +
  ylim(0,8) + xlim(-9.5,9.5)

#SSC vs HC
table <- table[order(table$pvalue_SSC_HC),]

ggplot(table, aes(x=log2FC_SSC_HC, y=-log10(pvalue_SSC_HC))) +
  ggtitle("SSc-Con vs HC") +
  xlab("log2 fold change") + 
  ylab("-log10 p-value") +
  geom_point(aes(colour = SSC_vs_HC_threshold),size=1) +
  scale_color_manual(values=c("steelblue2", "gray55","indianred1")) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) +
  geom_text_repel(data=head(table, 20), aes(label=metabolites), max.overlaps = Inf, size=4) +
  ylim(0,8) + xlim(-9.5,9.5)


#PAH vs SSC
table <- table[order(table$pvalue_PAH_SSC),]

ggplot(table, aes(x=log2FC_PAH_SSC, y=-log10(pvalue_PAH_SSC))) +
  ggtitle("SSc-PAH vs SSc-Con") +
  xlab("log2 fold change") + 
  ylab("-log10 p-value") +
  geom_point(aes(colour = PAH_vs_SSC_threshold),size=1) +
  scale_color_manual(values=c("steelblue2", "gray55","indianred1")) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) +
  geom_text_repel(data=head(table, 20), aes(label=metabolites), max.overlaps = Inf, size=4) +
  ylim(0,8) + xlim(-9.5,9.5)

#PF vs HC
table <- table[order(table$pvalue_PF_HC),]

ggplot(table, aes(x=log2FC_PF_HC, y=-log10(pvalue_PF_HC))) +
  ggtitle("SSc-PF vs HC") +
  xlab("log2 fold change") + 
  ylab("-log10 p-value") +
  geom_point(aes(colour = PF_vs_HC_threshold),size=1) +
  scale_color_manual(values=c("steelblue2", "gray55","indianred1")) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) +
  geom_text_repel(data=head(table, 20), aes(label=metabolites), max.overlaps = Inf, size=4) +
  ylim(0,8) + xlim(-9.5,9.5)

#PF vs SSC
table <- table[order(table$pvalue_PF_SSC),]

ggplot(table, aes(x=log2FC_PF_SSC, y=-log10(pvalue_PF_SSC))) +
  ggtitle("SSc-PF vs SSC") +
  xlab("log2 fold change") + 
  ylab("-log10 p-value") +
  geom_point(aes(colour = PF_vs_SSC_threshold),size=1) +
  scale_color_manual(values=c("steelblue2", "gray55","indianred1")) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) +
  geom_text_repel(data=head(table, 20), aes(label=metabolites), max.overlaps = Inf, size=4) +
  ylim(0,8) + xlim(-9.5,9.5)
#PAH vs PF
table <- table[order(table$pvalue_PAH_PF),]

ggplot(table, aes(x=log2FC_PAH_PF, y=-log10(pvalue_PAH_PF))) +
  ggtitle("SSc-PAH vs PF") +
  xlab("log2 fold change") + 
  ylab("-log10 p-value") +
  geom_point(aes(colour = PAH_vs_PF_threshold),size=1) +
  scale_color_manual(values=c("steelblue2", "gray55","indianred1")) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) +
  geom_text_repel(data=head(table, 20), aes(label=metabolites), max.overlaps = Inf, size=4) +
  ylim(0,8) + xlim(-9.5,9.5)



