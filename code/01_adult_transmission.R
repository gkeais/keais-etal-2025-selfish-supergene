
rm(list = ls())

setwd("~/projects/testacea_female_drive/code/")

library(tidyverse)
library(ggbeeswarm)
library(gridExtra)
library(grid)
library(glue)



#----------------------------------- import data -----------------------------------
df <- read.csv('../data/adult_transmission.csv')



#----------------------------------- prep data for plotting -------------------------
# split by female_genotype
drivers <- df[df$female_genotype == 'driver',]
length(unique(drivers$female))
controls <- df[df$female_genotype == 'control',]

# convert to wide format to calculate proportion Xd
df_prop_drivers <- spread(drivers, offspring_genotype, count)
df_prop_controls <- spread(controls, offspring_genotype, count)

# add proportion column for drivers (proportion non-Xb transmitted to offspring)
df_prop_drivers$nonXB <- (df_prop_drivers$XDXB + df_prop_drivers$XDY)
df_prop_drivers$XB <- (df_prop_drivers$XBXB + df_prop_drivers$XBY)
df_prop_drivers$sum <- (df_prop_drivers$XB + df_prop_drivers$nonXB)
df_prop_drivers$prop <- (df_prop_drivers$nonXB) / df_prop_drivers$sum

# add proportion column for controls (proportion non-Xb transmitted to offspring)
df_prop_controls$nonXB <- (df_prop_controls$XXB + df_prop_controls$XY)
df_prop_controls$XB <- (df_prop_controls$XBXB + df_prop_controls$XBY)
df_prop_controls$sum <- (df_prop_controls$XB + df_prop_controls$nonXB)
df_prop_controls$prop <- (df_prop_controls$nonXB) / df_prop_controls$sum

# pool data across all females
pooled_drivers <- drivers %>%
  group_by(offspring_genotype) %>%
  summarise(count = sum(count))

pooled_controls <- controls %>%
  group_by(offspring_genotype) %>%
  summarise(count = sum(count))

# factor offspring_genotype
pooled_drivers$offspring_genotype <- factor(pooled_drivers$offspring_genotype, levels = c('XBY', 'XBXB', 'XDY', 'XDXB'))
pooled_controls$offspring_genotype <- factor(pooled_controls$offspring_genotype, levels = c('XBY', 'XBXB', 'XY', 'XXB'))

# combine df_prop_drivers and df_prop_controls
df_prop <- rbind(df_prop_drivers[c('female', 'female_genotype', 'line', 'XB', 'nonXB', 'prop')], 
                 df_prop_controls[c('female', 'female_genotype', 'line', 'XB', 'nonXB', 'prop')])

df_prop <- df_prop[order(-df_prop$prop),]


# order females in the driver dataframe by their Xd proportion
drivers$female <- factor(drivers$female, levels = df_prop[df_prop$female_genotype == 'driver',]$female)
drivers$offspring_genotype <- factor(drivers$offspring_genotype, levels = c('XBY', 'XBXB', 'XDY', 'XDXB'))



#---------------------------- fig 1A: beeswarm plots of X transmission proportion ----------------------------
# figure 1A
df_prop$female_genotype <- factor(df_prop$female_genotype, levels = c("driver", "control"))


prop_female_genotype <- ggplot(df_prop, aes(x = female_genotype, y = prop)) +
  geom_beeswarm(cex = 3, size = 0.5) + 
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "black", linewidth = 0.5) +
  #stat_summary(fun.data = "mean_cl_normal", geom = "errorbar", colour = 'salmon', width = 0.25, lwd = 0.45) +
  stat_summary(fun=mean, geom='crossbar', colour='slateblue', width = 0.5, lwd = 0.4) +
#  geom_signif(comparisons = list(c("control", "driver")), annotations='***') +
  labs(x = "Female Genotype", y = "Proportion WT-eyed offspring") +
  scale_x_discrete(labels = c(expression('X'^'D'*'X'^'b'), expression('X'*'X'^'b'))) +
  theme_bw(base_size = 8) +
  labs(tag = "A")
prop_female_genotype

  
# figure S3
head(df_prop)
df_prop %>% group_by(line) %>% summarize(mean_prop = mean(prop))
# NAME CHANGES TO XD LINES:
# X^3 = X^D21 (X3 was collected in 2021)
# X^D = X^D12 (XD was collected in 2012)
# X^63 = X^D16 (X63 was collected in 2016)
# X3 > D > 63 > 12 (original names)
# D21 > D12 > D16 > 12 (new names)
df_prop$line <- factor(df_prop$line, levels = c("X3", "D", "63", "12"))
#prop_line <- ggplot(df_prop[df_prop$female_genotype != 'control',], aes(x = line, y = prop)) +
#  geom_beeswarm(cex = 3, size = 0.5) + 
#  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "black", linewidth = 0.5) +
#  stat_summary(fun=mean, geom='crossbar', colour='slateblue', width = 0.5, lwd = 0.4) +
#  labs(x = "Female Genotype", y = "Proportion WT-eyed offspring") +# y = expression("Proportion non-"*'X'^'B')) +
#  scale_x_discrete(labels = c(expression('X'^'D63'*'X'^'b'), expression('X'^'D12'*'X'^'b'), expression('X'^'D3'*'X'^'b'))) +
#  theme_bw(base_size = 8)
#prop_line

prop_line <- ggplot(df_prop, aes(x = line, y = prop)) +
  geom_beeswarm(cex = 3, size = 0.5) + 
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "black", linewidth = 0.5) +
  stat_summary(fun=mean, geom='crossbar', colour='slateblue', width = 0.5, lwd = 0.4) +
  labs(x = "Female Genotype", y = "Proportion WT-eyed offspring") +
  scale_x_discrete(labels = c(expression('X'^'D21'*'X'^'b'), 
                              expression('X'^'D12'*'X'^'b'), 
                              expression('X'^'D16'*'X'^'b'),
                              expression('XX'^'b'))) +
  theme_bw(base_size = 8)
prop_line



#----------------------------  plot pooled counts - drivers ----------------------------
# figure 1B
pooled_drivers_plot <- ggplot(pooled_drivers, aes(y = count, x = offspring_genotype, fill = offspring_genotype)) + 
  geom_bar(stat = "identity", alpha=0.6, colour="black", linewidth = 0.4, width = 0.6) +
  theme_classic(base_size = 8) +
  labs(x = "Offspring Genotype", y = "Count") + 
  scale_fill_manual(values=c("black", "grey", "steelblue", "darkorange")) + 
  theme(axis.line=element_line()) + 
  theme(axis.line.x = element_line(colour = 'black', linewidth=0.3, linetype='solid')) +
  theme(axis.line.y = element_line(colour = 'black', linewidth=0.3, linetype='solid')) +
  scale_y_continuous(expand = c(0,0)) + # expand brings the y axis line to the origin
  coord_cartesian(ylim = c(0, (max(pooled_drivers$count) + 10))) + # prevents black outline around bars from being cut off at the top
  scale_x_discrete(labels = c(expression('X'^'b'*'Y'), 
                              expression('X'^'b'*'X'^'b'), 
                              expression('X'^'D'*'Y'),
                              expression('X'^'D'*'X'^'b'))) +
  theme(legend.position = "none") +
  labs(tag = "B")
  
pooled_drivers_plot



#----------------------------  plot pooled counts - controls ----------------------------
# figure 1C
pooled_controls_plot <- ggplot(pooled_controls, aes(y = count, x = offspring_genotype, fill = offspring_genotype)) + 
  geom_bar(stat = "identity", alpha=0.6, colour="black", linewidth = 0.4, width = 0.6) +
  theme_classic(base_size = 8) +
  labs(x = "Offspring Genotype", y = "Count") + 
  scale_fill_manual(values=c("black", "grey", "#FCEDDA", "#EE4E34")) + 
  theme(axis.line=element_line()) + 
  theme(axis.line.x = element_line(colour = 'black', linewidth=0.3, linetype='solid')) +
  theme(axis.line.y = element_line(colour = 'black', linewidth=0.3, linetype='solid')) +
  scale_y_continuous(expand = c(0,0)) +
  coord_cartesian(ylim = c(0, (max(pooled_controls$count) + 10))) +
  scale_x_discrete(labels = c(expression('X'^'b'*'Y'), 
                              expression('X'^'b'*'X'^'b'), 
                              expression('XY'^phantom('|')),
                              expression('X'*'X'^'b'))) +
  theme(legend.position = "none") +
  labs(tag = "C")
pooled_controls_plot



#---------------------------- grid --------------------------------------------

# get date to append to figure output names
date <- Sys.Date()

# stacked bar
pdf(glue("../figures/{date}_genotype_proportions.pdf"), width=7.24, height=4)
grid.arrange(drivers_genotype_proportions, ncol=1)
dev.off()

# figure 1
pdf(glue("../figures/{date}_adult_transmission.pdf"), width=4.76, height=2)
grid.arrange(prop_female_genotype, pooled_drivers_plot, pooled_controls_plot, ncol=3)
dev.off()

# figure S2
pdf(glue("../figures/{date}_driver_lines.pdf"), width=(4.76/2), height=2)
grid.arrange(prop_line, ncol=1)
dev.off()



#------------------------- stats ------------------------------------------
# summarize the data
head(df_prop)
df_prop$total <- df_prop$XB + df_prop$nonXB
df_prop_summary <- df_prop %>%
  group_by(female_genotype) %>%
  summarize(mean = mean(prop, na.rm=TRUE), n = n(), 
            sum_XB = sum(XB), sum_nonXB = sum(nonXB), sum_total = sum(total),
            min = min(prop, na.rm=TRUE), max = max(prop, na.rm=TRUE))

df_prop_summary



#-------------------- figure 1a -----------------------------

# the number of WT-eyed offspring given n total offspring is binomially distributed.
# i.e., the probability of observing k WT-eyed offspring out of n total offspring given a success rate, p, 
# is given by the binomial probability mass function:
# nCk * p^k * (1-p)^(n-k)

# the proportion of WT-eyed offspring is: x/n
# note that x is binomially distributed, but x/n itself is not

# these data can be treated in at least three ways:
# (1) treat each proportion as a numerical value and perform a one sample t-test
# (2) pool offspring genotypes (XD and Xb) across all females and perform a binomial test (with large n, a Z-test is also appropriate)
# (3) test for an effect of female genotype (XDXb and XXb) on the transmission of the Xb using logistic regression
# (another approach would be a chi-squared test)

# Arguably, (1) is the easiest to understand, so that is what I've chosen.
# More details on this approach:
# The proportion of WT-eyed offspring produced by a given female is considered a sample proportion
# Each sample proportion is an estimate of the population proportion
# Therefore, the distribution of proportions in figure 1 is a sampling distribution of the sample proportions:
# (one for XDXb, and one for XXb)
hist(df_prop[df_prop$female_genotype == 'driver',]$prop)
hist(df_prop[df_prop$female_genotype == 'control',]$prop)
# These should approach a normal distribution as n increases, according to the CLT

# In this scenario, a one sample t-test tests the hypothesis that the mean of the 
# sampling distribution of the sample proportion is significantly different from 0.5
t.test(df_prop[df_prop$female_genotype == 'driver',]$prop, mu=0.5)
t.test(df_prop[df_prop$female_genotype == 'control',]$prop, mu=0.5)

out <- t.test(df_prop[df_prop$female_genotype == 'driver',]$prop, mu=0.5)
out$conf.int

# manually calculating the CIs
proportions <- df_prop[df_prop$female_genotype == 'driver',]$prop
n <- length(proportions)  # 200
mean_p <- mean(proportions)
sd_p <- sd(proportions)  
# t-critical value for 95% CI
t_crit <- qt(0.975, df = n - 1) 
# margin of Error
ME <- t_crit * (sd_p / sqrt(n))
# confidence interval
c(mean_p - ME, mean_p + ME)


#--- binomial test and Z-test on pooled data ---
n <- df_prop_summary[df_prop_summary$female_genotype == 'driver',]$sum_total
successes <- df_prop_summary[df_prop_summary$female_genotype == 'driver',]$sum_nonXB
binom.test(successes, n, p = 0.5, alternative = "two.sided")
prop.test(successes, n, p = 0.5, alternative = "two.sided", correct = TRUE)

n_control <- df_prop_summary[df_prop_summary$female_genotype == 'control',]$sum_total
successes_control <- df_prop_summary[df_prop_summary$female_genotype == 'control',]$sum_nonXB
binom.test(successes_control, n_control, p = 0.5, alternative = "two.sided")
prop.test(successes_control, n_control, p = 0.5, alternative = "two.sided", correct = TRUE)


#--- logistic regression ---
# https://stats.stackexchange.com/questions/26762/how-to-do-logistic-regression-in-r-when-outcome-is-fractional-a-ratio-of-two-co
# dependent variable is a 2 column matrix 
# the first column contains counts of 'successes' and the second column contains counts of 'failures'
# defining transmission of a nonXb chromosome as 'success'
head(df_prop)
head(cbind(df_prop$nonXB, df_prop$XB)) # note that this is the opposite column order to what is in df_prop
glm_female_genotype <- glm(cbind(nonXB, XB) ~ female_genotype, data = df_prop, family = binomial)
summary(glm_female_genotype)



#-------------------- figure S1 -----------------------------
df_prop$line
glm_female_genotype <- glm(cbind(nonXB, XB) ~ line, data = df_prop, family = binomial)
summary(glm_female_genotype)

head(df_prop)
lines_lm <- lm(prop ~ line, data = df_prop)
lines_aov <- aov(lines_lm)
summary(lines_aov)
tukey_test <- TukeyHSD(lines_aov)
tukey_test


glm_line <- glm(cbind(nonXB, XB) ~ line, data = df_prop[df_prop$female_genotype == 'driver',], family = binomial)
summary(glm_line)



# fin