
rm(list = ls())

setwd("~/projects/testacea_female_drive/code/")
library(tidyverse)
#library(viridis)
library(RColorBrewer)
library(gridExtra)
library(grid)
library(glue)



#------------------------------ import data -------------------------------
both <- read.csv('../data/model_output/2024-02-08_run1_both_nested_fitness_loop_ouput.csv')
both_parameters <- read.csv('../data/model_output/2024-02-08_run1_both_nested_fitness_loop_parameters.csv')
both_parameters
head(both)
dim(both)

male_only <- read.csv('../data/model_output/2024-02-08_run1_male_only_nested_fitness_loop_ouput.csv')
male_only_parameters <- read.csv('../data/model_output/2024-02-08_run1_male_only_nested_fitness_loop_parameters.csv')
male_only_parameters
head(male_only)
dim(male_only)


#------------------------------ colour palette -------------------------------
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")



#---------------------------- compute stability conditions ----------------------------
# why not to 1?
fW_vec = seq(from = 0, to = 0.999, length.out = 1000)
head(fW_vec)
length(fW_vec)

# PST rearranged in terms of fY
# fY computed as a function of fW (a vector ranging from 0 to 0.999 in increments of 0.001)
calc_fY_PST = function(n,m,fV,fU,fX){
  fY = (fX*((2*fU)-fV-(n*fV))) / ((1+m)*(1+n)*fV)
  return(fY)
}

# PSR rearranged in terms of fY
# fY computed as a function of fW (a vector ranging from 0 to 0.999 in increments of 0.001)
calc_fY_PSR = function(n,m,fV,fX,fW_vec){
  fY = (-fX*(1-n)*fV) / (((1+m)*(1-n)*fV)-(2*fW_vec*(1+m)))
  return(fY)
}

# three genotypes with constant fitness
fV = 1 # XdX
fU = 1 # XX
fX = 1 # XY



#----------------- make dataframes for PST and PSR under the "both" scenario ------------------------------

#----- compute a dataframe for PST -----
n = 0.317
m = 0.994
both_vec <- rep('both', 1000)
head(both_vec)
length(both_vec)

# XdY
fY_PST_both = calc_fY_PST(n,m,fV,fU,fX)
fY_PST_both

# fY_PST_both is a scalar, so rep to vector
fY_PST_both_vec <- rep(fY_PST_both, 1000)
head(fY_PST_both_vec)
length(fY_PST_both_vec)

PST_both <- data.frame(
  XdY = fY_PST_both_vec,
  XdXd = fW_vec,
  run = both_vec
)

dim(PST_both)
head(PST_both, n=10) # flat line at y=0.2753315
tail(PST_both, n=10)


#----- compute a dataframe for PSR -----
fY_PSR_both = calc_fY_PSR(n,m,fV,fX,fW_vec)
head(fY_PSR_both)
fY_PSR_both
length(fY_PSR_both)

PSR_both <- data.frame(
  XdY = fY_PSR_both,
  XdXd = fW_vec,
  run = both_vec
)

head(PSR_both, n=10)
tail(PSR_both, n=10)
dim(PSR_both)

PSR_both
PSR_both$XdY[PSR_both$XdY >= 1] <- NA
PSR_both
PSR_both$XdY[PSR_both$XdY < 0] <- NA
PSR_both


#----------------- compute dataframes for PST and PSR under the "male only" scenario ------------------------------
n = 0
m = 0.994

fY_PST_male = calc_fY_PST(n,m,fV,fU,fX)
fY_PST_male
fY_PSR_male = calc_fY_PSR(n,m,fV,fX,fW_vec)
fY_PSR_male

male_only_vec <- rep('male_only', 1000)
male_only_vec
fY_PST_male_vec <- rep(fY_PST_male, 1000)

PST_male <- data.frame(
  XdY = fY_PST_male_vec,
  XdXd = fW_vec,
  run = male_only_vec
)

PSR_male <- data.frame(
  XdY = fY_PSR_male,
  XdXd = fW_vec,
  run = male_only_vec
)

head(PSR_male)
PSR_male$XdY[PSR_male$XdY >= 1] <- NA
PSR_male
PSR_male$XdY[PSR_male$XdY < 0] <- NA
PSR_male



#---------------------------- rbind -------------------------
df_PST <- rbind(PST_both, PST_male)
table(df_PST$run)
df_PST$run = factor(df_PST$run, levels=c('male_only', 'both'))

df_PSR <- rbind(PSR_both, PSR_male)
table(df_PSR$run)
df_PSR$run = factor(df_PSR$run, levels=c('male_only', 'both'))



#------------------------------ plotting heatmap -------------------------------
XDXD_estimate <- 0.379

# Create a dataframe with labels and positions
labels <- data.frame(label = c("fixes", "polymorphism", "does not invade"),
                     x = c(0.881, 0.2, 0.5),
                     y = c(0.95, 0.75, 0.1))


#--- male only scenario (Fig. 4A) ---
head(male_only)
numerics_male_only <- ggplot() + 
  geom_raster(data=male_only, mapping = aes(x=fW, y=fY, fill=Xd_equil)) +
  scale_fill_gradientn(colours = myPalette(100), name=expression('X'^'D'*' frequency')) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_fixed() +
  labs(x = expression('X'^'D'*'X'^'D'*' fitness'), y = expression('X'^'D'*'Y'*' fitness')) +
  theme_bw(base_size=10) +
  theme(legend.justification = "top") + 
  geom_line(data = PSR_male, aes(x = XdXd, y = XdY)) +
  geom_line(data = PST_male, aes(x = XdXd, y = XdY)) +
  geom_text(data = labels, aes(x = x, y = y, label = label), size=2, colour='black') +
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(0.35, 'cm'), #change legend key height
        legend.key.width = unit(0.35, 'cm'), #change legend key width
        legend.title = element_text(size=8), #change legend title font size
        legend.text = element_text(size=8))  #change legend text font size
numerics_male_only


#--- both scenario (Fig. 4B) ---
head(both)
numerics_both <- ggplot() + 
  geom_raster(data=both, mapping = aes(x=fW, y=fY, fill=Xd_equil)) +
  scale_fill_gradientn(colours = myPalette(100), name=expression('X'^'D'*' frequency')) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_fixed() +
  labs(x = expression('X'^'D'*'X'^'D'*' fitness'), y = expression('X'^'D'*'Y'*' fitness')) +
  theme_bw(base_size=10) +
  theme(legend.justification = "top") + 
  geom_line(data = PSR_both, aes(x = XdXd, y = XdY)) +
  geom_line(data = PST_both, aes(x = XdXd, y = XdY)) +
  geom_text(data = labels, aes(x = x, y = y, label = label), size=2, colour='black') +
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(0.35, 'cm'), #change legend key height
        legend.key.width = unit(0.35, 'cm'), #change legend key width
        legend.title = element_text(size=8), #change legend title font size
        legend.text = element_text(size=8))  #change legend text font size
numerics_both


#------------------------------ plotting overlap (Fig. 4D) -------------------------------

#https://stackoverflow.com/questions/28586635/shade-region-between-two-lines-with-ggplot

# shading region of the male_only polymorphism region that is also polymorphic in the both scenario
# i.e. to the left of df_PSR where run == both and above df_PST where run == male_only
head(df_PSR)
a <- df_PSR[df_PSR$run == 'both', ]
b <- df_PST[df_PST$run == 'male_only', ]
c <- merge(a, b, by = "XdXd")

max(c$XdY.x, na.rm=TRUE)

stability_overlap <- ggplot(data=c, aes(x=XdXd)) +
  geom_line(data=c, aes(y=XdY.x), linewidth=1) +
  geom_line(data=c, aes(y=XdY.y), linewidth=1) +
  geom_ribbon(aes(ymin = XdY.y, ymax = pmax(XdY.x, XdY.y)), fill = "gray", alpha = 0.5) +
  # 0.514 is the x value where the "both" line reaches with 1 on the y axis
  geom_ribbon(data = c[c$XdXd<0.514,], aes(ymin = XdY.y, ymax=1), fill = "gray", alpha = 0.5) + 
  geom_line(data=df_PSR, aes(x=XdXd, y=XdY, group=run, colour=run), linewidth=1) +
  geom_line(data=df_PST, aes(x=XdXd, y=XdY, group=run, colour=run), linewidth=1) +
  scale_colour_manual(values=c('gray', 'plum3'), labels = c('male only', 'both')) +
  labs(x = expression('X'^'D'*'X'^'D'*' fitness'), y = expression('X'^'D'*'Y'*' fitness')) +
  scale_x_continuous(expand = c(0, 0), limits=c(0,1)) +
  scale_y_continuous(expand = c(0, 0), limits=c(0,1)) +
    theme_bw(base_size=10) +
  geom_vline(aes(xintercept = 0.379), linetype='dotted') + 
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(0.35, 'cm'), #change legend key height
        legend.key.width = unit(0.35, 'cm'), #change legend key width
        legend.title = element_text(size=8), #change legend title font size
        legend.text = element_text(size=8)) +  #change legend text font size
  coord_fixed() +
  theme(legend.justification = "top",
        legend.title = element_blank()) 
stability_overlap


#------------------ outcome proportions ---------------------------------------------------

both$scenario <- 'both'
male_only$scenario <- 'male only'
head(both)
head(male_only)

df <- rbind(male_only, both)

df$scenario <- factor(df$scenario, levels = c('male only', 'both'))

head(df)
table(df$scenario)

df <- df %>% mutate(outcome = case_when(
  Xd_equil == 1 ~ 'fixes',
  (Xd_equil > 1e-05 & Xd_equil < 1) ~ 'polymorphism',
  Xd_equil < 1e-05 ~ 'does not invade',
))

head(df)

101*101 # 10201 total area

df_area <- df %>% group_by(scenario) %>% count(outcome)

df_area %>% group_by(scenario) %>% summarise(sum = sum(n))

df_area$proportion <- df_area$n/10201
df_area

df_area$outcome <- factor(df_area$outcome, levels = c('does not invade', 'polymorphism', 'fixes'))
# check
df_area %>% group_by(scenario) %>% summarize(sum1 = sum(n), sum2 = sum(proportion))

head(df_area)
df_area$scenario <- factor(df_area$scenario, levels=c("both", "male only"))

outcome_proportions <- ggplot(df_area, aes(y = proportion, x = outcome, fill = scenario)) + 
  geom_bar(stat = "identity", position = 'dodge', alpha=0.85, colour="black", linewidth = 0.4, width = 0.8) +
  theme_classic(base_size = 10) +
  labs(x = "", y = "Proportion of area in A and B") + 
  scale_fill_manual(values=c('plum3', 'gray'), labels = c('both', 'male only')) +
  theme(axis.line=element_line()) + 
  theme(axis.line.x = element_line(colour = 'black', linewidth=0.3, linetype='solid')) +
  theme(axis.line.y = element_line(colour = 'black', linewidth=0.3, linetype='solid')) +
  scale_y_continuous(expand = c(0,0), limits = c(0,0.52)) + # expand brings the y axis line to the origin
  theme(legend.key.size = unit(0.35, 'cm'), #change legend key size
        #legend.key.height = unit(1, 'cm'), #change legend key height
        #legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_blank(), #change legend title font size
        legend.text = element_text(size=8)) +
  #theme(legend.position = c(0.9, 0.9)) +
  theme(legend.justification = c(0,1)) + # "top"
  guides(fill = guide_legend(reverse = TRUE)) + #flips the legend labels
  coord_flip() 
outcome_proportions


#----------------------------------------------------------------------------------------
#pdf(glue("../figures/{date}_model_numerics3.pdf"), width=12, height=6) # width=4.76, height=2)
library("cowplot")

#save <- ggdraw() +
#  draw_plot(numerics_male_only, x = 0, y = .5, width = .5, height = .5) +
#  draw_plot(numerics_both, x = .5, y = .5, width = .5, height = .5) +
#  draw_plot(outcome_proportions2, x = 0, y = 0.2, width = 0.6, height = 0.3) 
  #draw_plot(stability_overlap, x = 0, y = 0.2, width = 0.6, height = 0.3) +
#save

#http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/81-ggplot2-easy-way-to-mix-multiple-graphs-on-the-same-page/
save <- ggdraw() +
  draw_plot(numerics_male_only, x = 0, y = 0.5, width = 0.5, height = 0.5) +
  draw_plot(numerics_both, x = 0.5, y = 0.5, width = 0.5, height = 0.5) +
  draw_plot(outcome_proportions, x = 0, y = 0.1, width = 0.5, height = 0.3) + # width 0.5
  draw_plot(stability_overlap, x = 0.5, y = 0, width = 0.5, height = 0.5) +
  draw_plot_label(label = c("A", "B", "C", "D"), size = 14,
                  x = c(0, 0.5, 0, 0.5), y = c(1, 1, 0.5, 0.5))
save

# SAVED AS PDF WITH DEFAULT SIZE IN R (6.83 in x 4.83 in)


#-------------------- outdated plotting commands -------------------------------------------------------
date <- Sys.Date()

pdf(glue("../figures/{date}_save.pdf"), width=12, height=6)
#draw_plot_label(label = c("A", "B", "C"), size = 15, x = c(0, 0.5, 0), y = c(1, 1, 0.5))
#ggsave("../figures/{date}_model_numerics3.pdf", width=12, height=6)
grid.arrange(save, ncol = 1)
dev.off()



#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
pdf(glue("../figures/{date}_model_numerics.pdf"), width=12, height=6) # width=4.76, height=2)
grid.arrange(numerics_male_only, numerics_both, stability_overlap, ncol = 3)
dev.off()

pdf(glue("../figures/{date}_model_numerics2.pdf"), width=12, height=6) # width=4.76, height=2)
grid.arrange(numerics_male_only, numerics_both, outcome_proportions, stability_overlap, ncol = 2)
dev.off()

grid.arrange(numerics_male_only, numerics_both,
             outcome_proportions2, stability_overlap, nrow=2, ncol=2)
#--------------------------------------------------------------------------------------------
plot_grid(numerics_male_only, numerics_both, outcome_proportions1, nrow = 1, rel_widths = c(2/3, 2/3, 1/2))

gs <- list(numerics_male_only, numerics_both, outcome_proportions1)
grid.arrange(grobs = gs, widths = c(4,4,2), heights=c(6,6,6))
dev.off()
#----------------------------------------------------------------------------------------
pdf(glue("../figures/{date}_outcome_proportions.pdf"), width=4.76*(2/3), height=4)
grid.arrange(outcome_proportions, ncol = 1)
dev.off()


#------------------------------------- OLD save plot to pdf ------------------------
date <- Sys.Date()

# 6.26 x 3.73 (default)
pdf(glue("../figures/{date}_model_plot.pdf"), width=4.76, height=3.73)
grid.arrange(numerics, ncol = 1)
dev.off()

pdf(glue("../figures/{date}_numerics_with_estimates.pdf"), width=6.26, height=3.73)
grid.arrange(numerics_with_estimates, ncol = 1)
dev.off()

pdf(glue("../figures/{date}_outcome_proportions.pdf"), width=4.76*(2/3), height=4)
grid.arrange(area_bars, ncol = 1)
dev.off()


#----------------- combine plots -------------------

grid.arrange(all_plot, area_bars, ncol=2)

my_grobs <- list(all_plot, area_bars)
lay <- rbind(c(1,1), c(1,1), c(1,1), c(2,3), c(2,3))
lay
grid.arrange(grobs = my_grobs, layout_matrix = lay)

grid.arrange(grobs = my_grobs, widths = c(4, 1)) #layout_matrix = rbind(c(1, 2),c(1, 2))

#--------------------------------------------------------------------
df_area$outcome <- factor(df_area$outcome, levels = c('fixes', 'polymorphic', 'does not invade'))

area_bars <- ggplot(df_area, aes(y = proportion, x = outcome, fill = outcome)) + 
  geom_bar(stat = "identity", alpha=0.7, colour="black", size = 0.4, width = 0.8) +
  theme_classic(base_size = 8) +
  labs(x = "", y = "Proportion") + 
  #scale_fill_manual(values=c(myPalette(100)[1], myPalette(100)[50], myPalette(100)[100])) + 
  scale_fill_manual(values=c(myPalette(100)[100], myPalette(100)[50], myPalette(100)[1])) + 
  theme(axis.line=element_line()) + 
  theme(axis.line.x = element_line(colour = 'black', size=0.3, linetype='solid')) +
  theme(axis.line.y = element_line(colour = 'black', size=0.3, linetype='solid')) +
  scale_y_continuous(expand = c(0,0)) + # expand brings the y axis line to the origin
  coord_cartesian(xlim = c(0, (max(df_area$proportion) + 0.01))) + # prevents black outline around bars from being cut off at the top
  #geom_hline(yintercept = 0.5, linetype = "dashed", colour = "red") +
  guides(fill = guide_legend(title = "Offspring Genotype")) +
  #theme(axis.text.x = element_text(size = 8)) +
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  theme(legend.position = "none") +
  facet_wrap(~run) 

#facet_wrap(~run, ncol = 1, nrow = 3, strip.position="top") +
#coord_flip()
area_bars



stacked_bar <- ggplot(df_area, aes(y = n, x = run, fill = outcome)) + 
  geom_bar(position = "fill", stat = "identity", colour="black", size = 0.2, width = 0.7) +
  theme_bw(base_size = 10) +
  labs(x = "Scenario", y = "Proportion") + 
  scale_fill_manual(values=c(myPalette(100)[100], myPalette(100)[50], myPalette(100)[1])) + 
  theme(axis.line=element_line()) + 
  theme(axis.line.x = element_line(colour = 'black', linewidth = 0.2, linetype='solid')) +
  theme(axis.line.y = element_line(colour = 'black', linewidth = 0.2, linetype='solid')) +
  scale_y_continuous(breaks=seq(0,1.0,0.1), expand = c(0,0)) + # expand brings the y axis line to the origin
  coord_cartesian(ylim = c(0, 1.01)) + # prevents black outline around bars from being cut off at the top
  guides(fill = guide_legend(title = "Outcome")) +
  theme(legend.justification = "top")
stacked_bar

#-----------------------------------------------------------------------------------------------
male_only_plot <- ggplot(male_only, aes(x=XdXd, y=XdY, fill=Xd_equil)) + 
  geom_tile() +
  scale_fill_gradientn(colours = myPalette(100)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  #coord_equal() +
  ggtitle("male drive only (m = 0.999; n = 0)") +
  theme(plot.title = element_text(size = 10)) +
  labs(x = expression('X'^'D'*'Y'*' fitness (fY)'), y = expression('X'^'D'*'X'^'D'*' fitness (fW)')) + 
  theme_bw() +
  theme(legend.position="none")

male_only_plot

female_only_plot <- ggplot(female_only, aes(XdXd, XdY, fill=Xd_equil)) + 
  geom_tile() +
  scale_fill_gradientn(colours = myPalette(100)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  #coord_equal() +
  ggtitle("female drive only (m = 0; n = 0.24)") +
  theme(plot.title = element_text(size = 10)) +
  labs(x = expression('X'^'D'*'Y'*' fitness (fY)'), y = expression('X'^'D'*'X'^'D'*' fitness (fW)')) +
  theme_bw() +
  theme(legend.position="none")
female_only_plot


both_plot <- ggplot(both, aes(XdXd, XdY, fill=Xd_equil)) + 
  geom_tile() +
  scale_fill_gradientn(colours = myPalette(100)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  #coord_equal() +
  ggtitle("female drive only (m = 0; n = 0.24)") +
  theme(plot.title = element_text(size = 10)) +
  labs(x = expression('X'^'D'*'Y'*' fitness (fY)'), y = expression('X'^'D'*'X'^'D'*' fitness (fW)')) +
  theme_bw() +
  theme(legend.position="none")
both_plot


grid.arrange(male_only_plot, female_only_plot, both_plot, ncol=3)



#---------------------------- testing other colours -------------------------




# plotting the heat map with facets
numerics <- ggplot() + 
  #geom_tile() +
  geom_raster(data=df, mapping = aes(x=XdXd, y=XdY, fill=Xd_equil)) +
  scale_fill_gradientn(colours = myPalette(100), name=expression('X'^'D'*' frequency')) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_equal() +
  labs(x = expression('X'^'D'*'X'^'D'*' fitness'), y = expression('X'^'D'*'Y'*' fitness')) +
  theme_bw(base_size=10) +
  theme(legend.justification = "top") + 
  geom_line(data = df_PSR, aes(x = XdXd, y = XdY)) +
  geom_line(data = df_PST, aes(x = XdXd, y = XdY)) +
  geom_text(data = labels, aes(x = x, y = y, label = label), size=2, colour='black') +
  #geom_rect(data = rectangle_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), color='black', fill='transparent', linetype='dotted') +
  facet_wrap(~run, labeller = as_labeller(panel_labels)) +
  #geom_vline(data = df, aes(xintercept = 0.379), linetype='dotted') + 
  #geom_hline(data = df, aes(yintercept = 0.826), linetype='dotted') +
  #geom_vline(data = subset(df, run == "both"), aes(xintercept = 0.379), linetype='dotted') + # https://stackoverflow.com/questions/34686217/how-can-i-add-a-line-to-one-of-the-facets
  #geom_hline(data = subset(df, run == "both"), aes(yintercept = 0.826), linetype='dotted') +
  theme(panel.spacing = unit(1.5, "lines"), strip.background = element_blank()) +
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(0.35, 'cm'), #change legend key height
        legend.key.width = unit(0.35, 'cm'), #change legend key width
        legend.title = element_text(size=8), #change legend title font size
        legend.text = element_text(size=8)) + #change legend text font size
  labs(tag = "A")
numerics



ggplot(df, aes(XdXd, XdY, fill=Xd_equil)) + 
  geom_tile(alpha=0.8) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(plot.title = element_text(size = 10)) +
  labs(x = expression('X'^'D'*'X'^'D'*' fitness (fW)'), y = expression('X'^'D'*'Y'*' fitness (fY)')) +
  theme_bw(base_size=10) +
  facet_wrap(~run) + # scales='free'
  theme(panel.spacing = unit(1.5, "lines"), strip.background = element_blank())

ggplot(df, aes(XdXd, XdY, fill=Xd_equil)) + 
  geom_tile() +
  scale_fill_viridis(discrete=FALSE)

df[df$XdXd == 0.6,]
df[df$XdY == 0.6,]

ggplot(df, aes(x = XdXd, y = XdY, z = Xd_equil)) +
  geom_contour()

rectangle_data <- data.frame(
  run = factor("both", levels=c('male_only', 'female_only', 'both')),  # Facet where the rectangle will be added
  xmin = 0,   # Rectangle's x-min
  xmax = XDXD_estimate,    # Rectangle's x-max
  ymin = fY_PST_both,   # Rectangle's y-min
  ymax = XDY_estimate    # Rectangle's y-max
)
panel_labels <- c(
  'male_only'='male drive only (m = 0.99)',
  'female_only'='female drive only (n = 0.29)',
  'both'="both (m = 0.99; n = 0.29)"
)