#NADPH

#http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/76-add-p-values-and-significance-levels-to-ggplots/
library(tidyverse)
library(ggplot2)
library(plyr)
library(dplyr)
library(gapminder)
library(ggpubr)
library(ggsignif)
library(grid)

setwd("/Users/LaiYuanTe/Desktop/Term2/Research/8. Validation test/Final") 

NADPH_NADP = read.csv('NADPH_NADP.csv', check.names=FALSE)
NADPH <- NADPH_NADP[1:24,] #NADPH
NADP <- NADPH_NADP[25:48,] #NADP+


total <- ddply(NADPH_NADP, c("Group", "Type"), summarise,
               N    = length(Luminescence),
               mean = mean(Luminescence),
               sd   = sd(Luminescence),
               se   = sd / sqrt(N)
)
total
#NADPH
compare_means(Luminescence ~ Group, method = 't.test',  data = NADPH_NADP[1:24,]) #t-test
compare_means(Luminescence ~ Group,  data = NADPH_NADP[1:24,]) #Wilcoxon
#NADP
compare_means(Luminescence ~ Group, method = 't.test',  data = NADPH_NADP[25:48,]) #t-test
compare_means(Luminescence ~ Group,  data = NADPH_NADP[25:48,]) #Wilcoxon

my_comparisons <- list( c("HC", "SSC-CON"), c("SSC-CON", "SSC-PAH"), c("HC", "SSC-PAH") )


#Bar chart
total$Type <- factor(total$Type)
ggplot(total, aes(x=Type, y=mean, fill=Group)) + 
  geom_bar(width=0.7, position=position_dodge(0.85), stat="identity",
           colour="black", # Use black outlines,
           size=.3) +      # Thinner lines
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),
                size=.3,    # Thinner lines
                width=.2,
                position=position_dodge(.85)) + 
  scale_y_continuous(limits = c(0,5000), expand=c(0,0)) + 
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major.y = element_line(color = "black",
                                        size = 0.1,
                                        linetype = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5), text = element_text(size = 10), 
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5, size=12),
        axis.text.y = element_text(size=12)) + 
  xlab("Metabolite") +
  ylab("Luminescence") +
  labs(fill = "Group") + scale_fill_manual(values=c("gray87", "gray55", "gray36")) + 
  ggtitle("The fluctuation of NADPH and NADP+ in three patient groups") 


#Boxplot
NADPH %>% 
  ggplot(aes(x=Group,y=Luminescence,fill=Group)) +
  geom_boxplot(lwd=0.5,fatten=1.05,outlier.shape = NA) + geom_jitter(width=0.11,alpha=0.3) + theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5), text = element_text(size = 10), 
        axis.text.x = element_text(angle = 30, vjust = 0.5, hjust=0.5)) + 
  xlab("Group") + ylab("Luminescence") + 
  labs(fill = "Group") + scale_fill_manual(values=c("gray87", "gray55", "gray36")) +
  ggtitle("NADPH") +
  stat_compare_means(comparisons = my_comparisons) # Wilcoxon p-value
  
NADP %>% 
  ggplot(aes(x=Group,y=Luminescence,fill=Group)) +
  geom_boxplot(lwd=0.5,fatten=1.05,outlier.shape = NA) + geom_jitter(width=0.11,alpha=0.3) + theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5), text = element_text(size = 10), 
        axis.text.x = element_text(angle = 30, vjust = 0.5, hjust=0.5)) + 
  xlab("Group") + ylab("Luminescence") + 
  labs(fill = "Group") + scale_fill_manual(values=c("gray87", "gray55", "gray36")) +
  ggtitle("NADP+") +
  stat_compare_means(comparisons = my_comparisons) # Wilcoxon p-value

#NADPH ratio
NADPH_ratio = read.csv('NADPH ratio.csv', check.names=FALSE)
compare_means(Ratio ~ Group, method = 't.test',  data = NADPH_ratio) #t-test
compare_means(Ratio ~ Group,  data = NADPH_ratio) #Wilcoxon
my_comparisons <- list( c("HC", "SSC-CON"), c("SSC-CON", "SSC-PAH"), c("HC", "SSC-PAH") )

total <- ddply(NADPH_ratio, c("Group"), summarise,
               N    = length(Ratio),
               mean = mean(Ratio),
               sd   = sd(Ratio),
               se   = sd / sqrt(N)
)
total

#Bar chart
total$Group <- factor(total$Group)
ggplot(total, aes(x=Group, y=mean, fill=Group)) + 
  geom_bar(width=0.7, stat="identity",
           colour="black", # Use black outlines,
           size=.3) +      # Thinner lines
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),
                size=.3,    # Thinner lines
                width=.2,
                position=position_dodge(.8)) + 
  scale_y_continuous(limits = c(0,1.9), expand=c(0,0)) + 
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major.y = element_line(color = "black",
                                          size = 0.1,
                                          linetype = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5), text = element_text(size = 10), 
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5, size=12),
        axis.text.y = element_text(size=12)) + 
  xlab("Patient Group") +
  ylab("NADPH/NADP+ Ratio") +
  labs(fill = "Group") + scale_fill_manual(values=c("gray87", "gray55", "gray36")) + 
  ggtitle("NADPH/NADP+ ratio difference among three patient groups") 




#GSH
GSH = read.csv('GSH.csv', check.names=FALSE)
compare_means(Fluorescence ~ Group, method = 't.test',  data = GSH) #t-test
compare_means(Fluorescence ~ Group,  data = GSH) #Wilcoxon
my_comparisons <- list( c("HC", "SSC-CON"), c("SSC-CON", "SSC-PAH"), c("HC", "SSC-PAH") )

total <- ddply(GSH, c("Group"), summarise,
               N    = length(Fluorescence),
               mean = mean(Fluorescence),
               sd   = sd(Fluorescence),
               se   = sd / sqrt(N)
)
total

#Bar chart
total$Group <- factor(total$Group)

ggplot(total, aes(x=Group, y=mean, fill=Group)) + 
  geom_bar(width=0.75, stat="identity",
           colour="black", # Use black outlines,
           size=.3) +      # Thinner lines
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),
                size=.3,    # Thinner lines
                width=.2,
                position=position_dodge(.8)) + 
  scale_y_continuous(limits = c(0,4500000), expand=c(0,0)) + 
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major.y = element_line(color = "black",
                                          size = 0.1,
                                          linetype = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5), text = element_text(size = 10), 
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5, size=11.5),
        axis.text.y = element_text(size=11.5)) + 
  xlab("Patient Group") +
  ylab("Luminescence") +
  labs(fill = "Group") + scale_fill_manual(values=c("gray87", "gray55", "gray36")) + 
  ggtitle("Glutathione expression") 




#ROS
ROS = read.csv('ROS.csv', check.names=FALSE)
compare_means(Fluorescence ~ Group, method = 't.test',  data = ROS) #t-test
compare_means(Fluorescence ~ Group,  data = ROS) #Wilcoxon
my_comparisons <- list( c("HC", "SSC-CON"), c("SSC-CON", "SSC-PAH"), c("HC", "SSC-PAH") )

total <- ddply(ROS, c("Group"), summarise,
               N    = length(Fluorescence),
               mean = mean(Fluorescence),
               sd   = sd(Fluorescence),
               se   = sd / sqrt(N)
)
total

#Bar chart
total$Group <- factor(total$Group)

ggplot(total, aes(x=Group, y=mean, fill=Group)) + 
  geom_bar(width=0.75, stat="identity",
           colour="black", # Use black outlines,
           size=.3) +      # Thinner lines
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),
                size=.3,    # Thinner lines
                width=.2,
                position=position_dodge(.8)) + 
  scale_y_continuous(limits = c(0,180000), expand=c(0,0)) + 
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major.y = element_line(color = "black",
                                          size = 0.1,
                                          linetype = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5), text = element_text(size = 10), 
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5, size=11.5),
        axis.text.y = element_text(size=11.5)) + 
  xlab("Patient Group") +
  ylab("Luminescence") +
  labs(fill = "Group") + scale_fill_manual(values=c("gray87", "gray55", "gray36")) + 
  ggtitle("ROS expression") 


#ROS ratio
ROS_ratio = read.csv('ROS_ratio.csv', check.names=FALSE)
compare_means(Ratio ~ Group, method = 't.test',  data = ROS_ratio) #t-test
compare_means(Ratio ~ Group,  data = ROS_ratio) #Wilcoxon
my_comparisons <- list( c("HC", "SSC-CON"), c("SSC-CON", "SSC-PAH"), c("HC", "SSC-PAH") )

total <- ddply(ROS_ratio, c("Group"), summarise,
               N    = length(Ratio),
               mean = mean(Ratio),
               sd   = sd(Ratio),
               se   = sd / sqrt(N)
)
total

#Bar chart
total$Group <- factor(total$Group)

ggplot(total, aes(x=Group, y=mean, fill=Group)) + 
  geom_bar(width=0.75, stat="identity",
           colour="black", # Use black outlines,
           size=.3) +      # Thinner lines
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),
                size=.3,    # Thinner lines
                width=.2,
                position=position_dodge(.8)) + 
  scale_y_continuous(limits = c(0,0.18), expand=c(0,0)) + 
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major.y = element_line(color = "black",
                                          size = 0.1,
                                          linetype = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5), text = element_text(size = 10), 
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5, size=11.5),
        axis.text.y = element_text(size=11.5)) + 
  xlab("Patient Group") +
  ylab("530/590 on 355/460 ratio") +
  labs(fill = "Group") + scale_fill_manual(values=c("gray87", "gray55", "gray36")) + 
  ggtitle("530/590 on 355/460 ratio") 

