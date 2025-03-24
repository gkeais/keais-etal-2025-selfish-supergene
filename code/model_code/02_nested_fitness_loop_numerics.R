
rm(list = ls())

setwd("~/projects/testacea_female_drive/code/model_code/")


source("01_model_functions.R")
library(glue)
date <- Sys.Date()

# generations
generations = 1e4 #10,000
generations = 10

#--- genotype fitnesses ---
fU = 1 # 'XX'
fV = 1 # 'XdX'
fW = 1  # 'XdXd' (empirical estimate = 0.379)
fX = 1 # 'XY'
fY = 1 # 'XdY' (empirical estimate = 0.826)
fitnesses <- data.frame(fU,fV,fW,fX,fY)

#--- starting frequencies ---
# Ut = XX frequency
# Vt = XdX frequency
# Wt = XdXd frequency
# Xt = XY frequency
# Yt = XdY frequency

XX = 0.5 # Ut
XdX = 0 # Vt
XdXd = 0 # Wt
XY = 0.49 # Xt
XdY = 0.01 # Yt

#--- initiate matrix with starting frequencies ---
# make an empty matrix with a column for each genotype and rows equal to the number of generations
A = matrix(data = 0, ncol = 5, nrow = generations)

#--- fill the first row of the matrix ---
# 'XX', 'XdX', 'XdXd', 'XY', 'XdY'
A[1,] = c(XX, XdX, XdXd, XY, XdY)
head(A)

# the number of fitness increments (i/step)
step <- 100
# first genotype to vary (XdY)
genotype1 <- 'fY'
# second genotype to vary (XdXd)
genotype2 <- 'fW' 

# run number, for file name
run_number <- '1' # to save unique output files if running the same scenario more than once per day

#--- set scenario ---
#scenario <- 'male_only'
#n <- 0
#m <- 0.994
scenario <- 'both'
n <- 0.317
m <- 0.994
#scenario <- 'female_only'
#n <- 0.29
#m <- 0

#--- save parameters to csv ---
# note that fW and fY are passed in with fitness of 1 but these are each re-written in the loop starting from 0
parameter <- c('generations',
               'Ut', 'Vt', 'Wt', 'Xt', 'Yt',
               'fU', 'fV', 'fW', 'fX', 'fY', 
               'n', 'm', 'step', 'genotype1', 'genotype2', 'scenario')

value <- c(generations,
           XX, XdX, XdXd, XY, XdY,
           fU, fV, fW, fX, fY, 
           n, m, step, genotype1, genotype2, scenario)

df_parameters <- data.frame(parameter, value)
df_parameters
write.csv(df_parameters, glue("../../data/model_output/{date}_run{run_number}_{scenario}_nested_fitness_loop_parameters.csv"), row.names=FALSE)


#--- run the nested_fitness_loop function ---
#double_loop <- double_fitness_loop(A, generations, fitnesses, n, m, step, genotype1, genotype2)
nested_loop <- nested_fitness_loop(A, generations, fitnesses, n, m, step, genotype1, genotype2)
df <- as.data.frame(nested_loop)
#colnames(df) <- c('fY', 'fW', 'Xd_equil', 'XX', 'XdX', 'XdXd', 'XY', 'XdY')
head(df)
tail(df)
dim(df)

#--- write output data to csv ---
write.csv(df, glue("../../data/model_output/{date}_run{run_number}_{scenario}_nested_fitness_loop_ouput.csv"), row.names=FALSE)


# fin
