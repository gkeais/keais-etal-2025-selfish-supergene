
rm(list = ls())

setwd('~/projects/testacea_female_drive/code/')

library(tidyverse)
library(grid)
library(gridExtra)
library(ggbeeswarm)
library(glue)


embryos <- read.csv("../data/embryo_transmission.csv")
hatch <- read.csv("../data/egg_hatch.csv")
fitness <- read.csv("../data/female_fitness.csv") 

head(hatch)
table(hatch$female_genotype)



#-------------------------- prep data for plotting --------------------------------
head(embryos)

dim(embryos) # 284
embryos <- embryos[embryos$sample_type == 'embryo',]
dim(embryos) # 252
284-(8*4) # 8 replicates of the four adult genotypes removed

# filter to normal embryos
embryos_filt <- embryos %>%
  filter(embryo_genotype != 'failed') %>%
  filter(embryo_genotype != 'OY') %>%
  filter(embryo_genotype != 'XDXBY')

dim(embryos_filt) # 214

# groupby female
embryos_grouped <- embryos %>%
  group_by(female, embryo_genotype) %>%
  summarise(count=n())

embryos_grouped

# convert to wide format to calculate XB and XD numbers
embryos_wide <- spread(embryos_grouped, embryo_genotype, count)
embryos_wide


embryos_wide[is.na(embryos_wide)] <- 0
embryos_wide

embryos_wide$total_embryos <- rowSums(embryos_wide[2:8])
embryos_wide$total_success <- (embryos_wide$total_embryos - embryos_wide$failed)
embryos_wide$prop_success <- embryos_wide$total_success / embryos_wide$total_embryos 
embryos_wide

embryos_wide$XB <- (embryos_wide$XBXB + embryos_wide$XBY)
embryos_wide$XD <- (embryos_wide$XDXB + embryos_wide$XDY)


embryos_wide$normal <- embryos_wide$XB + embryos_wide$XD
embryos_wide$aneuploidies <- embryos_wide$total_success - embryos_wide$normal
embryos_wide
# remove genotype columns
embryos_wide_temp <- embryos_wide[c('female', 'XB', 'XD')]
embryos_wide_temp

# convert back to long for ggplot
embryos_XB_XD_long <- gather(embryos_wide_temp, embryo_genotype, count, XB:XD)
embryos_XB_XD_long

# factor females and order by drive strength
embryos_XB_XD_long$female <- factor(embryos_XB_XD_long$female, levels = c('female1', 'DB13', 'female6'))
embryos_XB_XD_long



#------------------------------ figure 2A: embryo transmission individual --------------------------------------------------
female_XB_XD_counts <- ggplot(embryos_XB_XD_long, aes(y = count, x = female, fill = embryo_genotype)) + 
  geom_bar(position = position_dodge(width=0.5), stat = "identity", alpha=0.6, colour="black", linewidth = 0.4, width = 0.5) +
  theme_classic(base_size = 8) +
  labs(x = "Female", y = "Count") + 
  scale_fill_manual(values=c('slategray', 'seashell2'), labels = c(expression(X^b), expression(X^D))) +
  # adding vjust = 1 to adjust for the extra height introduced by X^B (https://github.com/tidyverse/ggplot2/issues/3216)
  # expression() also distorts the key shape, so I manually set the key height and width below
  theme(legend.text = element_text(vjust = 1)) +
  theme(axis.line=element_line()) + 
  theme(axis.line.x = element_line(colour = 'black', linewidth=0.3, linetype='solid')) +
  theme(axis.line.y = element_line(colour = 'black', linewidth=0.3, linetype='solid')) +
  scale_y_continuous(expand = c(0,0)) + # expand brings the y axis line to the origin
  coord_cartesian(ylim = c(0, max(embryos_XB_XD_long$count) + 0.1)) + # prevents black outline around bars from being cut off at the top
  scale_x_discrete(labels = c('fe1', 'DB13', 'fe6')) +
  theme(legend.key.size = unit(0.3, 'cm'), #change legend key size
        legend.key.height = unit(0.3, 'cm'), #change legend key height
        legend.key.width = unit(0.3, 'cm'), #change legend key width
        legend.title = element_blank(), #change legend title font size
        legend.text = element_text(size=8)) +
  theme(legend.position = c(0.9, 0.9)) +
  labs(tag = "A")

female_XB_XD_counts



#---------------------------- fig 2B: embryo transmission pooled ----------------------------
pooled_embryo_genotypes <- embryos_grouped %>%
  group_by(embryo_genotype) %>%
  summarise(count = sum(count))
pooled_embryo_genotypes

pooled_embryo_genotypes_filt <- pooled_embryo_genotypes %>%
  filter(embryo_genotype != 'failed') %>%
  filter(embryo_genotype != 'OY') %>%
  filter(embryo_genotype != 'XDXBY')

pooled_embryo_genotypes_filt$embryo_genotype <- factor(pooled_embryo_genotypes_filt$embryo_genotype, 
                                                  levels = c("XBY", "XBXB", "XDY", "XDXB"))

pooled_embryo <- ggplot(pooled_embryo_genotypes_filt, aes(y = count, x = embryo_genotype, fill = embryo_genotype)) + 
  geom_bar(stat = "identity", alpha=0.6, colour="black", linewidth = 0.4, width = 0.6) +
  theme_classic(base_size = 8) +
  labs(x = "Embryo Genotype", y = "Count") + 
  scale_fill_manual(values=c("black", "grey", "steelblue", "darkorange")) + 
  theme(axis.line=element_line()) + 
  theme(axis.line.x = element_line(colour = 'black', linewidth=0.3, linetype='solid')) +
  theme(axis.line.y = element_line(colour = 'black', linewidth=0.3, linetype='solid')) +
  scale_y_continuous(expand = c(0,0)) + # expand brings the y axis line to the origin
  coord_cartesian(ylim = c(0, (max(pooled_embryo_genotypes$count) + 10))) + # prevents black outline around bars from being cut off at the top
  scale_x_discrete(labels = c(expression('X'^'b'*'Y'), 
                              expression('X'^'b'*'X'^'b'), 
                              expression('X'^'D'*'Y'),
                              expression('X'^'D'*'X'^'b'))) +
  theme(legend.position = "none") +
  labs(tag = "B")

pooled_embryo


# plot with aneuploidies and failed samples
pooled_embryo_genotypes
pooled_embryo_genotypes$embryo_genotype <- factor(pooled_embryo_genotypes$embryo_genotype, 
                                                       levels = c("XBY", "XBXB", "XDY", "XDXB", "XDXBY", "OY", "failed"))


pooled_embryo_aneuploidies <- ggplot(pooled_embryo_genotypes, aes(y = count, x = embryo_genotype, fill = embryo_genotype)) + 
  geom_bar(stat = "identity", alpha=0.6, colour="black", linewidth = 0.4, width = 0.6) +
  theme_classic(base_size = 8) +
  labs(x = "Embryo Genotype", y = "Count") + 
  scale_fill_manual(values=c("black", "grey", "steelblue", "darkorange", 'darkorchid', 'pink', 'red')) + 
  theme(axis.line=element_line()) + 
  theme(axis.line.x = element_line(colour = 'black', linewidth=0.3, linetype='solid')) +
  theme(axis.line.y = element_line(colour = 'black', linewidth=0.3, linetype='solid')) +
  scale_y_continuous(expand = c(0,0)) + # expand brings the y axis line to the origin
  coord_cartesian(ylim = c(0, (max(pooled_embryo_genotypes$count) + 10))) + # prevents black outline around bars from being cut off at the top
  scale_x_discrete(labels = c(expression('X'^'b'*'Y'), 
                              expression('X'^'b'*'X'^'b'), 
                              expression('X'^'D'*'Y'),
                              expression('X'^'D'*'X'^'b'),
                              expression('X'^'D'*'X'^'b'*'Y'),
                              expression('OY'^phantom('/')),
                              expression('failed'^phantom('/')))) +
  theme(legend.position = "none")

pooled_embryo_aneuploidies



#---------------------------- fig 2C: egg hatch plot ---------------------------
hatch$female_genotype <- factor(hatch$female_genotype, levels = c('XX', 'XdX'))
# figure 2C
prop_hatch <- ggplot(hatch, aes(x = female_genotype, y = prop_hatched)) +
  geom_beeswarm(cex = 3, size=0.5) +
  #geom_violin(alpha=0, width=0.5) +
  stat_summary(fun=mean, geom='crossbar', colour='slateblue', width=0.5, lwd=0.4) +
  labs(x = "Female Genotype", y = "Proportion hatched") +
  # adding phantom("/") so that XX aligns with XDX (https://github.com/tidyverse/ggplot2/issues/3216)
  scale_x_discrete(label = c(expression('XX'^phantom("/")), expression('X'^'D'*'X'), expression('X'^'D'*'X'^'D'))) +
  theme_bw(base_size=8) +
  labs(tag = "C")

prop_hatch



#---------------------------- fig 2D: female fitness ----------------------------
table(fitness$genotype)
fitness$genotype <- factor(fitness$genotype, levels = c('XBXB', 'XDXB', 'XDXD'))

female_fitness <- ggplot(fitness, aes(x = genotype, y = count)) +
  geom_beeswarm(cex = 3, size=0.5) + #aes(colour=geno)
  #geom_violin(alpha=0, width=0.5) +
  stat_summary(fun=mean, geom='crossbar', colour='slateblue', width=0.5, lwd=0.4)+
  labs(x = "Female Genotype", y = "Number of offspring") +
  scale_x_discrete(labels = c(expression('X'^'b'*'X'^'b'), expression('X'^'D'*'X'^'b'), expression('X'^'D'*'X'^'D'))) +
  theme_bw(base_size=8) +
  labs(tag = "D")

female_fitness



#------------------------------------- save plot to pdf ------------------------
date <- Sys.Date()


pdf(glue("../figures/{date}_embryo_transmission.pdf"), width=4.76*(2/3), height=4)
grid.arrange(female_XB_XD_counts, pooled_embryo, prop_hatch, female_fitness, ncol = 2)
dev.off()



#----------------------------------- stats -----------------------------------

#--------------------- fig 2A and B: embryo transmission -------------------
# fig 2b

embryos_wide_temp$total <- (embryos_wide_temp$XB + embryos_wide_temp$XD)
embryos_wide_temp$prop <- embryos_wide_temp$XD/embryos_wide_temp$total
embryos_wide_temp

pooled_embryo_genotypes

37+39 # XbXb + XbY = 76
65+73 # XDXb + XDY = 138

76+138 # total = 214

138/214 # proportion inheriting the XD
# binomial test of the number of embyros inheriting the XD
# where
# x = number of successes (138)
# n = number of trial (214)
# p = hypothesized probability of success (0.5)
binom.test(x = (138), n = (214), p = 0.5) # p-value = 2.701e-05



#--------------------- fig 2C: egg hatch ---------------------------
# fig 2c
# egg hatch averages by genotype
hatch %>%
  group_by(female_genotype) %>%
  summarise_at(vars(prop_hatched), list(hatch_mean = mean))

head(hatch)
hatch_glm <- glm(formula = cbind(hatched, unhatched) ~ female_genotype, data = hatch, family = binomial)
summary(hatch_glm)



#--------------------- fig 2D: female fitness ---------------------------
# female fitness averages by genotype

fitness %>%
  group_by(genotype) %>%
  summarise_at(vars(count), list(count_mean = mean))

51.5 / 136 # xdxd fitness = 0.38

# anova
head(fitness)

fitness_lm <- lm(count ~ genotype, data = fitness)
fitness_aov <- aov(fitness_lm)
summary(fitness_aov)
tukey_test <- TukeyHSD(fitness_aov)
tukey_test



# fin