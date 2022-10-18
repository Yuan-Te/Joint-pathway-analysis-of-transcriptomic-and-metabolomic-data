#Boxplot
#http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/76-add-p-values-and-significance-levels-to-ggplots/
library(tidyverse)
library(ggplot2)
library(dplyr)
library(gapminder)
library(ggpubr)
library(ggsignif)
library(grid)
setwd("/Users/LaiYuanTe/Desktop/Term2/Research/5. Enrichment analysis/Box plot_enrichment/HC_SSc_PAH") 


#----------------------------------------------------------------
metabolites = read.csv('OROTIDINE.csv', check.names=FALSE)

compare_means(orotidine ~ GROUP_NAME, method = 't.test',  data = metabolites)
my_comparisons <- list( c("HC", "SSC CON"), c("SSC CON", "SSC PAH"), c("HC", "SSC PAH") )

tiff("Orotidine.tiff", units="in", width=5, height=5, res=1000)
metabolites %>% 
  ggplot(aes(x=GROUP_NAME,y=metabolites[,3],fill=GROUP_NAME)) +
  geom_boxplot(lwd=0.5,fatten=1.05,outlier.shape = NA) + geom_jitter(width=0.11,alpha=0.3) + theme_bw() +
  ggtitle("Orotidine") +  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 10)) + 
  xlab("Group") + ylab("Relative Abundance") + 
  labs(fill = "Group") + scale_fill_manual(values=c("lightpink1", "lightgoldenrod", "lightskyblue")) +
  stat_compare_means(comparisons = my_comparisons, size=4, label.y=c(3.5,3.8,4.2))+ylim(0,5.5)
dev.off()

#----------------------------------------------------------------
metabolites = read.csv('2-aminoadipate.csv') #, check.names=FALSE

compare_means(X2.aminoadipate ~ GROUP_NAME, method = 't.test',  data = metabolites)
my_comparisons <- list( c("HC", "SSC CON"), c("SSC CON", "SSC PAH"), c("HC", "SSC PAH") )

tiff("2-aminoadipate.tiff", units="in", width=5, height=5, res=1000)
metabolites %>% 
  ggplot(aes(x=GROUP_NAME,y=metabolites[,3],fill=GROUP_NAME)) +
  geom_boxplot(lwd=0.5,fatten=1.05,outlier.shape = NA) + geom_jitter(width=0.11,alpha=0.3) + theme_bw() +
  ggtitle("2-aminoadipate") +  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 10)) + 
  xlab("Group") + ylab("Relative Abundance") + 
  labs(fill = "Group") + scale_fill_manual(values=c("lightpink1", "lightgoldenrod", "lightskyblue")) +
  stat_compare_means(comparisons = my_comparisons, size=4, label.y=c(3.5,3.8,4.2))+ylim(0,5.5)
  
dev.off()
#----------------------------------------------------------------
metabolites = read.csv('3beta,7alpha-dihydroxy-5-cholestenoate.csv') #, check.names=FALSE
colnames(metabolites)[3]

compare_means(X3beta.7alpha.dihydroxy.5.cholestenoate ~ GROUP_NAME, method = 't.test',  data = metabolites)
my_comparisons <- list( c("HC", "SSC CON"), c("SSC CON", "SSC PAH"), c("HC", "SSC PAH") )

tiff("3beta,7alpha-dihydroxy-5-cholestenoate.tiff", units="in", width=5, height=5, res=1000)
metabolites %>% 
  ggplot(aes(x=GROUP_NAME,y=metabolites[,3],fill=GROUP_NAME)) +
  geom_boxplot(lwd=0.5,fatten=1.05,outlier.shape = NA) + geom_jitter(width=0.11,alpha=0.3) + theme_bw() +
  ggtitle('3beta,7alpha-dihydroxy-5-cholestenoate') +  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 10)) + 
  xlab("Group") + ylab("Relative Abundance") + 
  labs(fill = "Group") + scale_fill_manual(values=c("lightpink1", "lightgoldenrod", "lightskyblue")) +
  stat_compare_means(comparisons = my_comparisons, size=4, label.y=c(3.5,3.8,4.2))+ylim(0,5.5)

dev.off()
  
#----------------------------------------------------------------
metabolites = read.csv('5,6-dihydrouracil.csv') #, check.names=FALSE
colnames(metabolites)[3]

compare_means(X5.6.dihydrouracil ~ GROUP_NAME, method = 't.test',  data = metabolites)
my_comparisons <- list( c("HC", "SSC CON"), c("SSC CON", "SSC PAH"), c("HC", "SSC PAH") )

tiff("5,6-dihydrouracil.tiff", units="in", width=5, height=5, res=1000)
metabolites %>% 
  ggplot(aes(x=GROUP_NAME,y=metabolites[,3],fill=GROUP_NAME)) +
  geom_boxplot(lwd=0.5,fatten=1.05,outlier.shape = NA) + geom_jitter(width=0.11,alpha=0.3) + theme_bw() +
  ggtitle('5,6-dihydrouracil') +  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 10)) + 
  xlab("Group") + ylab("Relative Abundance") + 
  labs(fill = "Group") + scale_fill_manual(values=c("lightpink1", "lightgoldenrod", "lightskyblue")) +
  stat_compare_means(comparisons = my_comparisons, size=4, label.y=c(3.5,3.8,4.2))+ylim(0,5.5)

dev.off()

#----------------------------------------------------------------
metabolites = read.csv('adenosine 5-monophosphate(AMP).csv') #, check.names=FALSE
colnames(metabolites)[3]

compare_means(adenosine.5..monophosphate..AMP. ~ GROUP_NAME, method = 't.test',  data = metabolites)
my_comparisons <- list( c("HC", "SSC CON"), c("SSC CON", "SSC PAH"), c("HC", "SSC PAH") )

tiff("adenosine 5-monophosphate(AMP).tiff", units="in", width=5, height=5, res=1000)
metabolites %>% 
  ggplot(aes(x=GROUP_NAME,y=metabolites[,3],fill=GROUP_NAME)) +
  geom_boxplot(lwd=0.5,fatten=1.05,outlier.shape = NA) + geom_jitter(width=0.11,alpha=0.3) + theme_bw() +
  ggtitle('adenosine 5-monophosphate(AMP)') +  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 10)) + 
  xlab("Group") + ylab("Relative Abundance") + 
  labs(fill = "Group") + scale_fill_manual(values=c("lightpink1", "lightgoldenrod", "lightskyblue")) +
  stat_compare_means(comparisons = my_comparisons, size=4, label.y=c(3.5,3.8,4.2))+ylim(0,5.5)

dev.off()


#----------------------------------------------------------------
metabolites = read.csv('arginine.csv') #, check.names=FALSE
colnames(metabolites)[3]

compare_means(arginine ~ GROUP_NAME, method = 't.test',  data = metabolites)
my_comparisons <- list( c("HC", "SSC CON"), c("SSC CON", "SSC PAH"), c("HC", "SSC PAH") )

tiff("arginine.tiff", units="in", width=5, height=5, res=1000)
metabolites %>% 
  ggplot(aes(x=GROUP_NAME,y=metabolites[,3],fill=GROUP_NAME)) +
  geom_boxplot(lwd=0.5,fatten=1.05,outlier.shape = NA) + geom_jitter(width=0.11,alpha=0.3) + theme_bw() +
  ggtitle('arginine') +  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 10)) + 
  xlab("Group") + ylab("Relative Abundance") + 
  labs(fill = "Group") + scale_fill_manual(values=c("lightpink1", "lightgoldenrod", "lightskyblue")) +
  stat_compare_means(comparisons = my_comparisons, size=4, label.y=c(3.5,3.8,4.2))+ylim(0,5.5)

dev.off()
  
#----------------------------------------------------------------
metabolites = read.csv('fumarate.csv') #, check.names=FALSE
colnames(metabolites)[3]

compare_means(fumarate ~ GROUP_NAME, method = 't.test',  data = metabolites)
my_comparisons <- list( c("HC", "SSC CON"), c("SSC CON", "SSC PAH"), c("HC", "SSC PAH") )

tiff("fumarate.tiff", units="in", width=5, height=5, res=1000)
metabolites %>% 
  ggplot(aes(x=GROUP_NAME,y=metabolites[,3],fill=GROUP_NAME)) +
  geom_boxplot(lwd=0.5,fatten=1.05,outlier.shape = NA) + geom_jitter(width=0.11,alpha=0.3) + theme_bw() +
  ggtitle('fumarate') +  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 10)) + 
  xlab("Group") + ylab("Relative Abundance") + 
  labs(fill = "Group") + scale_fill_manual(values=c("lightpink1", "lightgoldenrod", "lightskyblue")) +
  stat_compare_means(comparisons = my_comparisons, size=4, label.y=c(3.5,3.8,4.2))+ylim(0,5.5)

dev.off()

#----------------------------------------------------------------
metabolites = read.csv('glycine.csv') #, check.names=FALSE
colnames(metabolites)[3]

compare_means(glycine ~ GROUP_NAME, method = 't.test',  data = metabolites)
my_comparisons <- list( c("HC", "SSC CON"), c("SSC CON", "SSC PAH"), c("HC", "SSC PAH") )

tiff("glycine.tiff", units="in", width=5, height=5, res=1000)
metabolites %>% 
  ggplot(aes(x=GROUP_NAME,y=metabolites[,3],fill=GROUP_NAME)) +
  geom_boxplot(lwd=0.5,fatten=1.05,outlier.shape = NA) + geom_jitter(width=0.11,alpha=0.3) + theme_bw() +
  ggtitle('glycine') +  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 10)) + 
  xlab("Group") + ylab("Relative Abundance") + 
  labs(fill = "Group") + scale_fill_manual(values=c("lightpink1", "lightgoldenrod", "lightskyblue")) +
  stat_compare_means(comparisons = my_comparisons, size=4, label.y=c(3.5,3.8,4.2))+ylim(0,5.5)

dev.off()

#----------------------------------------------------------------
metabolites = read.csv('histidine.csv') #, check.names=FALSE
colnames(metabolites)[3]

compare_means(histidine ~ GROUP_NAME, method = 't.test',  data = metabolites)
my_comparisons <- list( c("HC", "SSC CON"), c("SSC CON", "SSC PAH"), c("HC", "SSC PAH") )

tiff("histidine.tiff", units="in", width=5, height=5, res=1000)
metabolites %>% 
  ggplot(aes(x=GROUP_NAME,y=metabolites[,3],fill=GROUP_NAME)) +
  geom_boxplot(lwd=0.5,fatten=1.05,outlier.shape = NA) + geom_jitter(width=0.11,alpha=0.3) + theme_bw() +
  ggtitle('histidine') +  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 10)) + 
  xlab("Group") + ylab("Relative Abundance") + 
  labs(fill = "Group") + scale_fill_manual(values=c("lightpink1", "lightgoldenrod", "lightskyblue")) +
  stat_compare_means(comparisons = my_comparisons, size=4, label.y=c(3.5,3.8,4.2))+ylim(0,5.5)

dev.off()

#----------------------------------------------------------------
metabolites = read.csv('lactate.csv') #, check.names=FALSE
colnames(metabolites)[3]

compare_means(lactate ~ GROUP_NAME, method = 't.test',  data = metabolites)
my_comparisons <- list( c("HC", "SSC CON"), c("SSC CON", "SSC PAH"), c("HC", "SSC PAH") )

tiff("lactate.tiff", units="in", width=5, height=5, res=1000)
metabolites %>% 
  ggplot(aes(x=GROUP_NAME,y=metabolites[,3],fill=GROUP_NAME)) +
  geom_boxplot(lwd=0.5,fatten=1.05,outlier.shape = NA) + geom_jitter(width=0.11,alpha=0.3) + theme_bw() +
  ggtitle('lactate') +  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 10)) + 
  xlab("Group") + ylab("Relative Abundance") + 
  labs(fill = "Group") + scale_fill_manual(values=c("lightpink1", "lightgoldenrod", "lightskyblue")) +
  stat_compare_means(comparisons = my_comparisons, size=4, label.y=c(3.5,3.8,4.2))+ylim(0,5.5)

dev.off()

#----------------------------------------------------------------
metabolites = read.csv('leucine.csv') #, check.names=FALSE
colnames(metabolites)[3]

compare_means(leucine ~ GROUP_NAME, method = 't.test',  data = metabolites)
my_comparisons <- list( c("HC", "SSC CON"), c("SSC CON", "SSC PAH"), c("HC", "SSC PAH") )

tiff("leucine.tiff", units="in", width=5, height=5, res=1000)
metabolites %>% 
  ggplot(aes(x=GROUP_NAME,y=metabolites[,3],fill=GROUP_NAME)) +
  geom_boxplot(lwd=0.5,fatten=1.05,outlier.shape = NA) + geom_jitter(width=0.11,alpha=0.3) + theme_bw() +
  ggtitle('leucine') +  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 10)) + 
  xlab("Group") + ylab("Relative Abundance") + 
  labs(fill = "Group") + scale_fill_manual(values=c("lightpink1", "lightgoldenrod", "lightskyblue")) +
  stat_compare_means(comparisons = my_comparisons, size=4, label.y=c(3.5,3.8,4.2))+ylim(0,5.5)

dev.off()

#----------------------------------------------------------------
metabolites = read.csv('N-carbamoylaspartate.csv') #, check.names=FALSE
colnames(metabolites)[3]

compare_means(N.carbamoylaspartate ~ GROUP_NAME, method = 't.test',  data = metabolites)
my_comparisons <- list( c("HC", "SSC CON"), c("SSC CON", "SSC PAH"), c("HC", "SSC PAH") )

tiff("N-carbamoylaspartate.tiff", units="in", width=5, height=5, res=1000)
metabolites %>% 
  ggplot(aes(x=GROUP_NAME,y=metabolites[,3],fill=GROUP_NAME)) +
  geom_boxplot(lwd=0.5,fatten=1.05,outlier.shape = NA) + geom_jitter(width=0.11,alpha=0.3) + theme_bw() +
  ggtitle('N-carbamoylaspartate') +  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 10)) + 
  xlab("Group") + ylab("Relative Abundance") + 
  labs(fill = "Group") + scale_fill_manual(values=c("lightpink1", "lightgoldenrod", "lightskyblue")) +
  stat_compare_means(comparisons = my_comparisons, size=4, label.y=c(3.5,3.8,4.2))+ylim(0,5.5)

dev.off()

#----------------------------------------------------------------
metabolites = read.csv('phenylacetate.csv') #, check.names=FALSE
colnames(metabolites)[3]

compare_means(phenylacetate ~ GROUP_NAME, method = 't.test',  data = metabolites)
my_comparisons <- list( c("HC", "SSC CON"), c("SSC CON", "SSC PAH"), c("HC", "SSC PAH") )

tiff("phenylacetate.tiff", units="in", width=5, height=5, res=1000)
metabolites %>% 
  ggplot(aes(x=GROUP_NAME,y=metabolites[,3],fill=GROUP_NAME)) +
  geom_boxplot(lwd=0.5,fatten=1.05,outlier.shape = NA) + geom_jitter(width=0.11,alpha=0.3) + theme_bw() +
  ggtitle('phenylacetate') +  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 10)) + 
  xlab("Group") + ylab("Relative Abundance") + 
  labs(fill = "Group") + scale_fill_manual(values=c("lightpink1", "lightgoldenrod", "lightskyblue")) +
  stat_compare_means(comparisons = my_comparisons, size=4, label.y=c(3.5,3.8,4.2))+ylim(0,5.5)

dev.off()

#----------------------------------------------------------------
metabolites = read.csv('pyridoxate.csv') #, check.names=FALSE
colnames(metabolites)[3]

compare_means(pyridoxate ~ GROUP_NAME, method = 't.test',  data = metabolites)
my_comparisons <- list( c("HC", "SSC CON"), c("SSC CON", "SSC PAH"), c("HC", "SSC PAH") )

tiff("pyridoxate.tiff", units="in", width=5, height=5, res=1000)
metabolites %>% 
  ggplot(aes(x=GROUP_NAME,y=metabolites[,3],fill=GROUP_NAME)) +
  geom_boxplot(lwd=0.5,fatten=1.05,outlier.shape = NA) + geom_jitter(width=0.11,alpha=0.3) + theme_bw() +
  ggtitle('pyridoxate') +  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 10)) + 
  xlab("Group") + ylab("Relative Abundance") + 
  labs(fill = "Group") + scale_fill_manual(values=c("lightpink1", "lightgoldenrod", "lightskyblue")) +
  stat_compare_means(comparisons = my_comparisons, size=4, label.y=c(3.5,3.8,4.2))+ylim(0,5.5)

dev.off()

#----------------------------------------------------------------
metabolites = read.csv('retinal.csv') #, check.names=FALSE
colnames(metabolites)[3]

compare_means(retinal ~ GROUP_NAME, method = 't.test',  data = metabolites)
my_comparisons <- list( c("HC", "SSC CON"), c("SSC CON", "SSC PAH"), c("HC", "SSC PAH") )

tiff("retinal.tiff", units="in", width=5, height=5, res=1000)
metabolites %>% 
  ggplot(aes(x=GROUP_NAME,y=metabolites[,3],fill=GROUP_NAME)) +
  geom_boxplot(lwd=0.5,fatten=1.05,outlier.shape = NA) + geom_jitter(width=0.11,alpha=0.3) + theme_bw() +
  ggtitle('retinal') +  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 10)) + 
  xlab("Group") + ylab("Relative Abundance") + 
  labs(fill = "Group") + scale_fill_manual(values=c("lightpink1", "lightgoldenrod", "lightskyblue")) +
  stat_compare_means(comparisons = my_comparisons, size=4, label.y=c(3.5,3.8,4.2))+ylim(0,5.5)

dev.off()

#----------------------------------------------------------------
metabolites = read.csv('serine.csv') #, check.names=FALSE
colnames(metabolites)[3]

compare_means(serine ~ GROUP_NAME, method = 't.test',  data = metabolites)
my_comparisons <- list( c("HC", "SSC CON"), c("SSC CON", "SSC PAH"), c("HC", "SSC PAH") )

tiff("serine.tiff", units="in", width=5, height=5, res=1000)
metabolites %>% 
  ggplot(aes(x=GROUP_NAME,y=metabolites[,3],fill=GROUP_NAME)) +
  geom_boxplot(lwd=0.5,fatten=1.05,outlier.shape = NA) + geom_jitter(width=0.11,alpha=0.3) + theme_bw() +
  ggtitle('serine') +  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 10)) + 
  xlab("Group") + ylab("Relative Abundance") + 
  labs(fill = "Group") + scale_fill_manual(values=c("lightpink1", "lightgoldenrod", "lightskyblue")) +
  stat_compare_means(comparisons = my_comparisons, size=4, label.y=c(3.5,3.8,4.2))+ylim(0,5.5)

dev.off()

#----------------------------------------------------------------
metabolites = read.csv('sphinganine.csv') #, check.names=FALSE
colnames(metabolites)[3]

compare_means(sphinganine ~ GROUP_NAME, method = 't.test',  data = metabolites)
my_comparisons <- list( c("HC", "SSC CON"), c("SSC CON", "SSC PAH"), c("HC", "SSC PAH") )

tiff("sphinganine.tiff", units="in", width=5, height=5, res=1000)
metabolites %>% 
  ggplot(aes(x=GROUP_NAME,y=metabolites[,3],fill=GROUP_NAME)) +
  geom_boxplot(lwd=0.5,fatten=1.05,outlier.shape = NA) + geom_jitter(width=0.11,alpha=0.3) + theme_bw() +
  ggtitle('sphinganine') +  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 10)) + 
  xlab("Group") + ylab("Relative Abundance") + 
  labs(fill = "Group") + scale_fill_manual(values=c("lightpink1", "lightgoldenrod", "lightskyblue")) +
  stat_compare_means(comparisons = my_comparisons, size=4, label.y=c(3.5,3.8,4.2))+ylim(0,5.5)

dev.off()

#----------------------------------------------------------------
metabolites = read.csv('sphingomyelin (d18_2_18_1).csv') #, check.names=FALSE
colnames(metabolites)[3]

compare_means(sphingomyelin..d18.2.18.1.. ~ GROUP_NAME, method = 't.test',  data = metabolites)
my_comparisons <- list( c("HC", "SSC CON"), c("SSC CON", "SSC PAH"), c("HC", "SSC PAH") )

tiff("sphingomyelin (d18_2_18_1).tiff", units="in", width=5, height=5, res=1000)
metabolites %>% 
  ggplot(aes(x=GROUP_NAME,y=metabolites[,3],fill=GROUP_NAME)) +
  geom_boxplot(lwd=0.5,fatten=1.05,outlier.shape = NA) + geom_jitter(width=0.11,alpha=0.3) + theme_bw() +
  ggtitle('sphingomyelin (d18:2/18:1)*') +  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 10)) + 
  xlab("Group") + ylab("Relative Abundance") + 
  labs(fill = "Group") + scale_fill_manual(values=c("lightpink1", "lightgoldenrod", "lightskyblue")) +
  stat_compare_means(comparisons = my_comparisons, size=4, label.y=c(3.5,3.8,4.2))+ylim(0,5.5)

dev.off()

#----------------------------------------------------------------
metabolites = read.csv('sphingosine.csv') #, check.names=FALSE
colnames(metabolites)[3]

compare_means(sphingosine ~ GROUP_NAME, method = 't.test',  data = metabolites)
my_comparisons <- list( c("HC", "SSC CON"), c("SSC CON", "SSC PAH"), c("HC", "SSC PAH") )

tiff("sphingosine.tiff", units="in", width=5, height=5, res=1000)
metabolites %>% 
  ggplot(aes(x=GROUP_NAME,y=metabolites[,3],fill=GROUP_NAME)) +
  geom_boxplot(lwd=0.5,fatten=1.05,outlier.shape = NA) + geom_jitter(width=0.11,alpha=0.3) + theme_bw() +
  ggtitle('sphingosine') +  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 10)) + 
  xlab("Group") + ylab("Relative Abundance") + 
  labs(fill = "Group") + scale_fill_manual(values=c("lightpink1", "lightgoldenrod", "lightskyblue")) +
  stat_compare_means(comparisons = my_comparisons, size=4, label.y=c(3.5,3.8,4.2))+ylim(0,5.5)

dev.off()

#----------------------------------------------------------------
metabolites = read.csv('uracil.csv') #, check.names=FALSE
colnames(metabolites)[3]

compare_means(uracil ~ GROUP_NAME, method = 't.test',  data = metabolites)
my_comparisons <- list( c("HC", "SSC CON"), c("SSC CON", "SSC PAH"), c("HC", "SSC PAH") )

tiff("uracil.tiff", units="in", width=5, height=5, res=1000)
metabolites %>% 
  ggplot(aes(x=GROUP_NAME,y=metabolites[,3],fill=GROUP_NAME)) +
  geom_boxplot(lwd=0.5,fatten=1.05,outlier.shape = NA) + geom_jitter(width=0.11,alpha=0.3) + theme_bw() +
  ggtitle('uracil') +  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 10)) + 
  xlab("Group") + ylab("Relative Abundance") + 
  labs(fill = "Group") + scale_fill_manual(values=c("lightpink1", "lightgoldenrod", "lightskyblue")) +
  stat_compare_means(comparisons = my_comparisons, size=4, label.y=c(3.5,3.8,4.2))+ylim(0,5.5)

dev.off()

#----------------------------------------------------------------
metabolites = read.csv('vanillylmandelate (VMA).csv') #, check.names=FALSE
colnames(metabolites)[3]

compare_means(vanillylmandelate..VMA. ~ GROUP_NAME, method = 't.test',  data = metabolites)
my_comparisons <- list( c("HC", "SSC CON"), c("SSC CON", "SSC PAH"), c("HC", "SSC PAH") )

tiff("vanillylmandelate (VMA).tiff", units="in", width=5, height=5, res=1000)
metabolites %>% 
  ggplot(aes(x=GROUP_NAME,y=metabolites[,3],fill=GROUP_NAME)) +
  geom_boxplot(lwd=0.5,fatten=1.05,outlier.shape = NA) + geom_jitter(width=0.11,alpha=0.3) + theme_bw() +
  ggtitle('vanillylmandelate (VMA)') +  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 10)) + 
  xlab("Group") + ylab("Relative Abundance") + 
  labs(fill = "Group") + scale_fill_manual(values=c("lightpink1", "lightgoldenrod", "lightskyblue")) +
  stat_compare_means(comparisons = my_comparisons, size=4, label.y=c(3.5,3.8,4.2))+ylim(0,5.5)

dev.off()

#----------------------------------------------------------------
metabolites = read.csv('1-methyl-4-imidazoleacetate.csv') #, check.names=FALSE
colnames(metabolites)[3]

compare_means(X1.methyl.4.imidazoleacetate ~ GROUP_NAME, method = 't.test',  data = metabolites)
my_comparisons <- list( c("HC", "SSC CON"), c("SSC CON", "SSC PAH"), c("HC", "SSC PAH") )

tiff("1-methyl-4-imidazoleacetate.tiff", units="in", width=5, height=5, res=1000)
metabolites %>% 
  ggplot(aes(x=GROUP_NAME,y=metabolites[,3],fill=GROUP_NAME)) +
  geom_boxplot(lwd=0.5,fatten=1.05,outlier.shape = NA) + geom_jitter(width=0.11,alpha=0.3) + theme_bw() +
  ggtitle('1-methyl-4-imidazoleacetate') +  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 10)) + 
  xlab("Group") + ylab("Relative Abundance") + 
  labs(fill = "Group") + scale_fill_manual(values=c("lightpink1", "lightgoldenrod", "lightskyblue")) +
  stat_compare_means(comparisons = my_comparisons, size=4, label.y=c(3.5,3.8,4.2))+ylim(0,5.5)

dev.off()

#----------------------------------------------------------------
metabolites = read.csv('4-acetamidobutanoate.csv') #, check.names=FALSE
colnames(metabolites)[3]

compare_means(X4.acetamidobutanoate ~ GROUP_NAME, method = 't.test',  data = metabolites)
my_comparisons <- list( c("HC", "SSC CON"), c("SSC CON", "SSC PAH"), c("HC", "SSC PAH") )

tiff("4-acetamidobutanoate.tiff", units="in", width=5, height=5, res=1000)
metabolites %>% 
  ggplot(aes(x=GROUP_NAME,y=metabolites[,3],fill=GROUP_NAME)) +
  geom_boxplot(lwd=0.5,fatten=1.05,outlier.shape = NA) + geom_jitter(width=0.11,alpha=0.3) + theme_bw() +
  ggtitle('4-acetamidobutanoate') +  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 10)) + 
  xlab("Group") + ylab("Relative Abundance") + 
  labs(fill = "Group") + scale_fill_manual(values=c("lightpink1", "lightgoldenrod", "lightskyblue")) +
  stat_compare_means(comparisons = my_comparisons, size=4, label.y=c(3.5,3.8,4.2))+ylim(0,5.5)

dev.off()

#----------------------------------------------------------------
metabolites = read.csv('7-methylurate.csv') #, check.names=FALSE
colnames(metabolites)[3]

compare_means(X7.methylurate ~ GROUP_NAME, method = 't.test',  data = metabolites)
my_comparisons <- list( c("HC", "SSC CON"), c("SSC CON", "SSC PAH"), c("HC", "SSC PAH") )

tiff("7-methylurate.tiff", units="in", width=5, height=5, res=1000)
metabolites %>% 
  ggplot(aes(x=GROUP_NAME,y=metabolites[,3],fill=GROUP_NAME)) +
  geom_boxplot(lwd=0.5,fatten=1.05,outlier.shape = NA) + geom_jitter(width=0.11,alpha=0.3) + theme_bw() +
  ggtitle('7-methylurate') +  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 10)) + 
  xlab("Group") + ylab("Relative Abundance") + 
  labs(fill = "Group") + scale_fill_manual(values=c("lightpink1", "lightgoldenrod", "lightskyblue")) +
  stat_compare_means(comparisons = my_comparisons, size=4, label.y=c(3.5,3.8,4.2))+ylim(0,5.5)

dev.off()

#----------------------------------------------------------------
metabolites = read.csv('kynurenine.csv') #, check.names=FALSE
colnames(metabolites)[3]

compare_means(kynurenine ~ GROUP_NAME, method = 't.test',  data = metabolites)
my_comparisons <- list( c("HC", "SSC CON"), c("SSC CON", "SSC PAH"), c("HC", "SSC PAH") )

tiff("kynurenine.tiff", units="in", width=5, height=5, res=1000)
metabolites %>% 
  ggplot(aes(x=GROUP_NAME,y=metabolites[,3],fill=GROUP_NAME)) +
  geom_boxplot(lwd=0.5,fatten=1.05,outlier.shape = NA) + geom_jitter(width=0.11,alpha=0.3) + theme_bw() +
  ggtitle('kynurenine') +  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 10)) + 
  xlab("Group") + ylab("Relative Abundance") + 
  labs(fill = "Group") + scale_fill_manual(values=c("lightpink1", "lightgoldenrod", "lightskyblue")) +
  stat_compare_means(comparisons = my_comparisons, size=4, label.y=c(3.5,3.8,4.2))+ylim(0,5.5)

dev.off()

#----------------------------------------------------------------
metabolites = read.csv('N1-Methyl-2-pyridone-5-carboxamide.csv') #, check.names=FALSE
colnames(metabolites)[3]

compare_means(N1.Methyl.2.pyridone.5.carboxamide ~ GROUP_NAME, method = 't.test',  data = metabolites)
my_comparisons <- list( c("HC", "SSC CON"), c("SSC CON", "SSC PAH"), c("HC", "SSC PAH") )

tiff("N1-Methyl-2-pyridone-5-carboxamide.tiff", units="in", width=5, height=5, res=1000)
metabolites %>% 
  ggplot(aes(x=GROUP_NAME,y=metabolites[,3],fill=GROUP_NAME)) +
  geom_boxplot(lwd=0.5,fatten=1.05,outlier.shape = NA) + geom_jitter(width=0.11,alpha=0.3) + theme_bw() +
  ggtitle('N1-Methyl-2-pyridone-5-carboxamide') +  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 10)) + 
  xlab("Group") + ylab("Relative Abundance") + 
  labs(fill = "Group") + scale_fill_manual(values=c("lightpink1", "lightgoldenrod", "lightskyblue")) +
  stat_compare_means(comparisons = my_comparisons, size=4, label.y=c(3.5,3.8,4.2))+ylim(0,5.5)

dev.off()

#----------------------------------------------------------------
metabolites = read.csv('N-acetylneuraminate.csv') #, check.names=FALSE
colnames(metabolites)[3]

compare_means(N.acetylneuraminate ~ GROUP_NAME, method = 't.test',  data = metabolites)
my_comparisons <- list( c("HC", "SSC CON"), c("SSC CON", "SSC PAH"), c("HC", "SSC PAH") )

tiff("N-acetylneuraminate.tiff", units="in", width=5, height=5, res=1000)
metabolites %>% 
  ggplot(aes(x=GROUP_NAME,y=metabolites[,3],fill=GROUP_NAME)) +
  geom_boxplot(lwd=0.5,fatten=1.05,outlier.shape = NA) + geom_jitter(width=0.11,alpha=0.3) + theme_bw() +
  ggtitle('N-acetylneuraminate') +  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 10)) + 
  xlab("Group") + ylab("Relative Abundance") + 
  labs(fill = "Group") + scale_fill_manual(values=c("lightpink1", "lightgoldenrod", "lightskyblue")) +
  stat_compare_means(comparisons = my_comparisons, size=4, label.y=c(3.5,3.8,4.2))+ylim(0,5.5)

dev.off()

#----------------------------------------------------------------
metabolites = read.csv('nicotinamide.csv') #, check.names=FALSE
colnames(metabolites)[3]

compare_means(nicotinamide ~ GROUP_NAME, method = 't.test',  data = metabolites)
my_comparisons <- list( c("HC", "SSC CON"), c("SSC CON", "SSC PAH"), c("HC", "SSC PAH") )

tiff("nicotinamide.tiff", units="in", width=5, height=5, res=1000)
metabolites %>% 
  ggplot(aes(x=GROUP_NAME,y=metabolites[,3],fill=GROUP_NAME)) +
  geom_boxplot(lwd=0.5,fatten=1.05,outlier.shape = NA) + geom_jitter(width=0.11,alpha=0.3) + theme_bw() +
  ggtitle('nicotinamide') +  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 10)) + 
  xlab("Group") + ylab("Relative Abundance") + 
  labs(fill = "Group") + scale_fill_manual(values=c("lightpink1", "lightgoldenrod", "lightskyblue")) +
  stat_compare_means(comparisons = my_comparisons, size=4, label.y=c(3.5,3.8,4.2))+ylim(0,5.5)

dev.off()

#----------------------------------------------------------------
metabolites = read.csv('quinolinate.csv') #, check.names=FALSE
colnames(metabolites)[3]

compare_means(quinolinate ~ GROUP_NAME, method = 't.test',  data = metabolites)
my_comparisons <- list( c("HC", "SSC CON"), c("SSC CON", "SSC PAH"), c("HC", "SSC PAH") )

tiff("quinolinate.tiff", units="in", width=5, height=5, res=1000)
metabolites %>% 
  ggplot(aes(x=GROUP_NAME,y=metabolites[,3],fill=GROUP_NAME)) +
  geom_boxplot(lwd=0.5,fatten=1.05,outlier.shape = NA) + geom_jitter(width=0.11,alpha=0.3) + theme_bw() +
  ggtitle('quinolinate') +  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 10)) + 
  xlab("Group") + ylab("Relative Abundance") + 
  labs(fill = "Group") + scale_fill_manual(values=c("lightpink1", "lightgoldenrod", "lightskyblue")) +
  stat_compare_means(comparisons = my_comparisons, size=4, label.y=c(3.5,3.8,4.2))+ylim(0,5.5)

dev.off()

#----------------------------------------------------------------
metabolites = read.csv('sucrose.csv') #, check.names=FALSE
colnames(metabolites)[3]

compare_means(sucrose ~ GROUP_NAME, method = 't.test',  data = metabolites)
my_comparisons <- list( c("HC", "SSC CON"), c("SSC CON", "SSC PAH"), c("HC", "SSC PAH") )

tiff("sucrose.tiff", units="in", width=5, height=5, res=1000)
metabolites %>% 
  ggplot(aes(x=GROUP_NAME,y=metabolites[,3],fill=GROUP_NAME)) +
  geom_boxplot(lwd=0.5,fatten=1.05,outlier.shape = NA) + geom_jitter(width=0.11,alpha=0.3) + theme_bw() +
  ggtitle('sucrose') +  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 10)) + 
  xlab("Group") + ylab("Relative Abundance") + 
  labs(fill = "Group") + scale_fill_manual(values=c("lightpink1", "lightgoldenrod", "lightskyblue")) +
  stat_compare_means(comparisons = my_comparisons, size=4)

dev.off()

#----------------------------------------------------------------
metabolites = read.csv('gulonate.csv') #, check.names=FALSE
colnames(metabolites)[3]

compare_means(gulonate. ~ GROUP_NAME, method = 't.test',  data = metabolites)
my_comparisons <- list( c("HC", "SSC CON"), c("SSC CON", "SSC PAH"), c("HC", "SSC PAH") )

tiff("gulonate.tiff", units="in", width=5, height=5, res=1000)
metabolites %>% 
  ggplot(aes(x=GROUP_NAME,y=metabolites[,3],fill=GROUP_NAME)) +
  geom_boxplot(lwd=0.5,fatten=1.05,outlier.shape = NA) + geom_jitter(width=0.11,alpha=0.3) + theme_bw() +
  ggtitle('gulonate*') +  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 10)) + 
  xlab("Group") + ylab("Relative Abundance") + 
  labs(fill = "Group") + scale_fill_manual(values=c("lightpink1", "lightgoldenrod", "lightskyblue")) +
  stat_compare_means(comparisons = my_comparisons, size=4, label.y=c(3.5,3.8,4.2))+ylim(0,5.5)

dev.off()

#----------------------------------------------------------------
metabolites = read.csv('dimethylglycine.csv') #, check.names=FALSE
colnames(metabolites)[3]

compare_means(dimethylglycine ~ GROUP_NAME, method = 't.test',  data = metabolites)
my_comparisons <- list( c("HC", "SSC CON"), c("SSC CON", "SSC PAH"), c("HC", "SSC PAH") )

tiff("dimethylglycine.tiff", units="in", width=5, height=5, res=1000)
metabolites %>% 
  ggplot(aes(x=GROUP_NAME,y=metabolites[,3],fill=GROUP_NAME)) +
  geom_boxplot(lwd=0.5,fatten=1.05,outlier.shape = NA) + geom_jitter(width=0.11,alpha=0.3) + theme_bw() +
  ggtitle('dimethylglycine') +  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 10)) + 
  xlab("Group") + ylab("Relative Abundance") + 
  labs(fill = "Group") + scale_fill_manual(values=c("lightpink1", "lightgoldenrod", "lightskyblue")) +
  stat_compare_means(comparisons = my_comparisons, size=4, label.y=c(3.5,3.8,4.2))+ylim(0,5.5)

dev.off()

#----------------------------------------------------------------
metabolites = read.csv('carnitine.csv') #, check.names=FALSE
colnames(metabolites)[3]

compare_means(carnitine ~ GROUP_NAME, method = 't.test',  data = metabolites)
my_comparisons <- list( c("HC", "SSC CON"), c("SSC CON", "SSC PAH"), c("HC", "SSC PAH") )

tiff("carnitine.tiff", units="in", width=5, height=5, res=1000)
metabolites %>% 
  ggplot(aes(x=GROUP_NAME,y=metabolites[,3],fill=GROUP_NAME)) +
  geom_boxplot(lwd=0.5,fatten=1.05,outlier.shape = NA) + geom_jitter(width=0.11,alpha=0.3) + theme_bw() +
  ggtitle('carnitine') +  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 10)) + 
  xlab("Group") + ylab("Relative Abundance") + 
  labs(fill = "Group") + scale_fill_manual(values=c("lightpink1", "lightgoldenrod", "lightskyblue")) +
  stat_compare_means(comparisons = my_comparisons, size=4, label.y=c(3.5,3.8,4.2))+ylim(0,5.5)

dev.off()

#----------------------------------------------------------------
metabolites = read.csv('arabitol_xylitol.csv') #, check.names=FALSE
colnames(metabolites)[3]

compare_means(arabitol.xylitol ~ GROUP_NAME, method = 't.test',  data = metabolites)
my_comparisons <- list( c("HC", "SSC CON"), c("SSC CON", "SSC PAH"), c("HC", "SSC PAH") )

tiff("arabitol_xylitol.tiff", units="in", width=5, height=5, res=1000)
metabolites %>% 
  ggplot(aes(x=GROUP_NAME,y=metabolites[,3],fill=GROUP_NAME)) +
  geom_boxplot(lwd=0.5,fatten=1.05,outlier.shape = NA) + geom_jitter(width=0.11,alpha=0.3) + theme_bw() +
  ggtitle('arabitol/xylitol') +  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 10)) + 
  xlab("Group") + ylab("Relative Abundance") + 
  labs(fill = "Group") + scale_fill_manual(values=c("lightpink1", "lightgoldenrod", "lightskyblue")) +
  stat_compare_means(comparisons = my_comparisons, size=4, label.y=c(3.5,3.8,4.2))+ylim(0,5.5)

dev.off()

#----------------------------------------------------------------
metabolites = read.csv('7-methylxanthine.csv') #, check.names=FALSE
colnames(metabolites)[3]

compare_means(X7.methylxanthine ~ GROUP_NAME, method = 't.test',  data = metabolites)
my_comparisons <- list( c("HC", "SSC CON"), c("SSC CON", "SSC PAH"), c("HC", "SSC PAH") )

tiff("7-methylxanthine.tiff", units="in", width=5, height=5, res=1000)
metabolites %>% 
  ggplot(aes(x=GROUP_NAME,y=metabolites[,3],fill=GROUP_NAME)) +
  geom_boxplot(lwd=0.5,fatten=1.05,outlier.shape = NA) + geom_jitter(width=0.11,alpha=0.3) + theme_bw() +
  ggtitle('7-methylxanthine') +  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 10)) + 
  xlab("Group") + ylab("Relative Abundance") + 
  labs(fill = "Group") + scale_fill_manual(values=c("lightpink1", "lightgoldenrod", "lightskyblue")) +
  stat_compare_means(comparisons = my_comparisons, size=4, label.y=c(3.5,3.8,4.2))+ylim(0,5.5)

dev.off()

#----------------------------------------------------------------
metabolites = read.csv('2-hydroxyglutarate.csv') #, check.names=FALSE
colnames(metabolites)[3]

compare_means(X2.hydroxyglutarate ~ GROUP_NAME, method = 't.test',  data = metabolites)
my_comparisons <- list( c("HC", "SSC CON"), c("SSC CON", "SSC PAH"), c("HC", "SSC PAH") )

tiff("2-hydroxyglutarate.tiff", units="in", width=5, height=5, res=1000)
metabolites %>% 
  ggplot(aes(x=GROUP_NAME,y=metabolites[,3],fill=GROUP_NAME)) +
  geom_boxplot(lwd=0.5,fatten=1.05,outlier.shape = NA) + geom_jitter(width=0.11,alpha=0.3) + theme_bw() +
  ggtitle('2-hydroxyglutarate') +  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 10)) + 
  xlab("Group") + ylab("Relative Abundance") + 
  labs(fill = "Group") + scale_fill_manual(values=c("lightpink1", "lightgoldenrod", "lightskyblue")) +
  stat_compare_means(comparisons = my_comparisons, size=4, label.y=c(3.5,3.8,4.2))+ylim(0,5.5)

dev.off()
  
  
  
  
  
  
  
  
  
  
  


  
  
  
  




