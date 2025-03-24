
rm(list = ls())

setwd('~/projects/testacea_female_drive/code/')

library(tidyverse)
library(grid)
library(gridExtra)
library(ggbeeswarm)
library(glue)



#---------------------------- load data ----------------------------
females <- read.csv("../data/female_drive_suppressor.csv", header=TRUE)
fathers <- read.csv("../data/female_drive_suppressor_fathers.csv", header=TRUE)

# remove father 5, as there are no daughters of father 5
fathers$father_ID
fathers <- fathers[fathers$father_ID != '5',]
# re-factor
fathers$father_ID <- factor(fathers$father_ID)
fathers$father_ID



#---------------------------- check data types ----------------------------
head(females)
summary(females) # female is numeric, needs to be character
females$female_ID <- as.character(females$female_ID)
females$father_ID <- as.character(females$father_ID)
females$father_genotype <- as.factor(females$father_genotype)
summary(females) # check

summary(fathers)
fathers <- fathers[order(-fathers$sex_ratio),]
fathers
fathers$father_genotype <- factor(fathers$father_genotype)
fathers$father_ID <- factor(fathers$father_ID, levels=fathers$father_ID)



#-----------------------------------------------------------------------
females$sum <- (females$XDXB + females$XDY + females$XBXB + females$XBY)

dim(females)
females
females <- females[females$sum >= 30, ] # one female produced fewer than 30 offspring
dim(females)

females$XD <- females$XDXB + females$XDY
females$XB <- females$XBXB + females$XBY
females$prop <- females$XD / females$sum



#-------------------------- fig. S4a: father sex ratio ---------------------------------
fathers_sex_ratio <- ggplot(fathers, aes(y = sex_ratio, x = father_ID, fill = father_genotype)) + 
  geom_bar(stat = "identity", alpha=0.6, colour="black", linewidth = 0.4, width = 0.6) +
  theme_classic(base_size = 8) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "black", linewidth = 0.3) +
  labs(x = "Father", y = "Proportion female offspring") + 
  scale_fill_manual(values=c("black", "grey")) + 
  theme(axis.line=element_line()) + 
  theme(axis.line.x = element_line(colour = 'black', linewidth=0.3, linetype='solid')) +
  theme(axis.line.y = element_line(colour = 'black', linewidth=0.3, linetype='solid')) +
  scale_y_continuous(expand = c(0,0)) + # expand brings the y axis line to the origin
  coord_cartesian(ylim = c(0, 1)) + # prevents black outline around bars from being cut off at the top
  theme(legend.key.size = unit(0.3, 'cm'), #change legend key size
        legend.key.height = unit(0.3, 'cm'), #change legend key height
        legend.key.width = unit(0.3, 'cm'), #change legend key width
        legend.title = element_blank(), #change legend title font size
        legend.text = element_text(size=5),
        legend.position = c(0.8, 0.9)) +
  labs(tag = "A")

fathers_sex_ratio



#--------------------- fig. S4b: bee swarm female drive ------------------------------------
female_sup <- ggplot(females, aes(x = father_genotype, y = prop)) +
  geom_jitter(aes(shape=father_ID), cex = 2, width = 0.2) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "black", linewidth = 0.5) +
  stat_summary(fun=mean, geom='crossbar', colour='slateblue', width=0.4, lwd=0.3)+
  labs(x = "Father's Genotype", y = expression("Proportion dark-eyed offspring")) +
  theme_bw(base_size=8) +
  scale_x_discrete(labels = c('suppressed', 'unsuppressed')) +
  guides(shape=guide_legend(title="Father")) +
  theme(legend.justification = "top") + 
  labs(tag = "B")

female_sup



#-----------------------------------------------------------------------
date <- Sys.Date()
pdf(glue("../figures/{date}_female_suppression.pdf"), width=5, height=2.5)
grid.arrange(fathers_sex_ratio, female_sup, ncol = 2)
dev.off()



#---------------- stats ---------------------------------------
head(females)
female_sup_glm <- glm(formula = cbind(XD, XB) ~ father_genotype, data = females, family = binomial)
summary(female_sup_glm)



# fin
