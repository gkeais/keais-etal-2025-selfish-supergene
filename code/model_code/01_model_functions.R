# Functions for the numerical analysis of the female drive model 

rm(list = ls())

#---------------------------------------------------------------------------------
# This function contains the model recursion equations
#     x = a row of the A matrix (specifically, A[t-1,])
#     fitnesses = a dataframe containing the fitnesses of each genotype
#     n = the female drive strength (numeric value) 
#     m = the male drive strength (numeric value)

model = function(x, fitnesses, n, m){
  
  fU = fitnesses$fU # 'XX'
  fV = fitnesses$fV # 'XdX'
  fW = fitnesses$fW # 'XdXd'
  fX = fitnesses$fX # 'XY'
  fY = fitnesses$fY # 'XdY'
  
  y = rep(0, 5) # make a vector y of length 5 filled with 0s = [0, 0, 0, 0, 0]
  
  t = ((1/2)*fU*fX*x[1]*x[4]) + 
    (((1-n)/4)*fV*fX*x[2]*x[4]) +
    (((1+n)/4)*fV*fX*x[2]*x[4]) + 
    (((1+m)/2)*fU*fY*x[1]*x[5]) +
    ((((1+m)*(1-n))/4)*fV*fY*x[2]*x[5]) + 
    ((1/2)*fW*fX*x[3]*x[4]) +
    ((((1+n)*(1+m))/4)*fV*fY*x[2]*x[5]) +
    (((1+m)/2)*fW*fY*x[3]*x[5]) +
    ((1/2)*fU*fX*x[1]*x[4]) + 
    (((1-n)/4)*fV*fX*x[2]*x[4]) + 
    (((1-m)/2)*fU*fY*x[1]*x[5]) + 
    ((((1-m)*(1-n))/4)*fV*fY*x[2]*x[5]) +
    (((1+n)/4)*fV*fX*x[2]*x[4]) + 
    ((((1-m)*(1+n))/4)*fV*fY*x[2]*x[5]) + 
    (((1-m)/2)*fW*fY*x[3]*x[5]) +
    ((1/2)*fW*fX*x[3]*x[4])
  
  y[1] <- (((1/2)*fU*fX*x[1]*x[4]) + 
    (((1-n)/4)*fV*fX*x[2]*x[4])) / t
  
  y[2] <- ((((1+n)/4)*fV*fX*x[2]*x[4]) + 
    (((1+m)/2)*fU*fY*x[1]*x[5]) +
    ((((1+m)*(1-n))/4)*fV*fY*x[2]*x[5]) + 
    ((1/2)*fW*fX*x[3]*x[4])) / t
  
  y[3] <- (((((1+n)*(1+m))/4)*fV*fY*x[2]*x[5]) +
    (((1+m)/2)*fW*fY*x[3]*x[5])) / t 
  
  y[4] <- (((1/2)*fU*fX*x[1]*x[4]) + 
    (((1-n)/4)*fV*fX*x[2]*x[4]) + 
    (((1-m)/2)*fU*fY*x[1]*x[5]) + 
    ((((1-m)*(1-n))/4)*fV*fY*x[2]*x[5])) / t
    
  y[5] <- ((((1+n)/4)*fV*fX*x[2]*x[4]) + 
    ((((1-m)*(1+n))/4)*fV*fY*x[2]*x[5]) + 
    (((1-m)/2)*fW*fY*x[3]*x[5]) +
    ((1/2)*fW*fX*x[3]*x[4])) / t
  
  print(t)
  return(y)
}



#---------------------------------------------------------------------------------
# This function runs the model for a set number of generations
#     A = data matrix to be populated (its first row is initiated with the starting genotype frequencies)
#     generations = the number of generations the model will run for (integer value) 
#     fitnesses = a dataframe containing the fitnesses of each genotype (columns are indexed in the model function)
#     n = the female drive strength (numeric value)
#     m = the male drive strength (numeric value)

run = function(A, generations, fitnesses, n, m){
  for (t in 2:generations){
    A[t,] = model(A[t-1,], fitnesses, n, m)
  }
  return(A)
}



#---------------------------------------------------------------------------------
# This function runs the model over a range of fitness values for one pre-defined genotype 
# The fitness of the chosen genotype ranges from 0 to 1 in 0.01 steps
# Each run of the model is for x generations (usually 10,000)
# The frequency of all the genotypes (and the summed Xd frequency) at the end of each run is saved into a row of a matrix
# 
# The function takes the following parameters:
#     A = data matrix to be populated (its first row is initiated with the starting genotype frequencies)
#     generations = the number of generations the model will run for (integer value) 
#     fitnesses = a dataframe containing the fitnesses of each genotype (columns are indexed in the model function)
#     n = the female drive strength (numeric value) 
#     m = the male drive strength (numeric value)
#     step = the number of fitness increments over which to iterate (integer value)
#     genotype = the genotype whose fitness will be changed with each step (and that is defined by i/step)
#
# Returns a 101x7 matrix (if step = 100)

fitness_loop <- function(A, generations, fitnesses, n, m, step, genotype){
  # make an empty matrix with 7 columns and 101 rows (step+1)
  # the 7 columns are listed below
  # nrow = step+1 because step starts at 0 and dataframes start at 1 (0+1 = 1, the start of the dataframe)
  A2 <- matrix(data = 0, ncol = 7, nrow = step+1)
  head(A2)
  for (i in 0:step){
    # set the fitness of the chosen genotype to i/step (e.g. for i=1 and a step of 100: 1/100 = 0.01)
    fitnesses[genotype] <- i/step
    # call the model and convert to a dataframe whose columns are the five genotypes
    output <- run(A, generations, fitnesses, n, m)
    df <- as.data.frame(output)
    colnames(df) <- c('XX', 'XdX', 'XdXd', 'XY', 'XdY')
    # compute the "equilibrium" Xd frequency, the sum of the Xds in the final row of the output dataframe
    Xd_equil <- (df[generations,]$XdXd) + (df[generations,]$XdX/2) + (df[generations,]$XdY)
    # define a row of the A2 matrix column-by-column
    A2[i+1,][1] <- i/step
    A2[i+1,][2] <- Xd_equil
    A2[i+1,][3] <- df[generations,]$XX
    A2[i+1,][4] <- df[generations,]$XdX
    A2[i+1,][5] <- df[generations,]$XdXd
    A2[i+1,][6] <- df[generations,]$XY
    A2[i+1,][7] <- df[generations,]$XdY
  }
  return(A2)
}



#---------------------------------------------------------------------------------
# This function is deprecated - the nested_fitness_loop function below is a faster version
# This function runs the model over a range of fitness values for two pre-defined genotypes (nested for loop)
# This is done by calling the fitness_loop function (which calls the run function, which calls the model function)
# The fitnesses of the chosen genotypes each range from 0 to 1 in 0.01 steps
# Each run of the model is for x generations (usually 10,000)
# The frequency of all the genotypes (and the summed Xd frequency) at the end of each run is saved into a row of a matrix
# 
# The function takes the following parameters:
#     A = data matrix to be populated (its first row is initiated with the starting genotype frequencies)
#     generations = the number of generations the model will run for (integer value) 
#     fitnesses = a dataframe containing the fitnesses of each genotype (columns are indexed in the model function)
#     n = the female drive strength (numeric value) 
#     m = the male drive strength (numeric value)
#     step = the number of fitness increments over which to iterate (integer value)
#     genotype1 = the first genotype whose fitness will be changed with each step (and that is defined by i/step)
#     genotype2 = the second genotype whose fitness will be changed with each step (and that is defined by i/step)

double_fitness_loop <- function(A, generations, fitnesses, n, m, step, genotype1, genotype2){
  start.time <- Sys.time()
  for (i in 0:step){
    fitnesses[genotype1] <- i/step
    # for the first run (i == 0), output1 does not exist yet, so there is nothing to rbind
    if (i == 0) {
      # fitness_loop returns the matrix A2 
      # A2 is a 101x7 matrix with columns: genotype2, 'Xd_equil', 'XX', 'XdX', 'XdXd', 'XY', 'XdY'
      output1 <- fitness_loop(A, generations, fitnesses, n, m, step, genotype2) # 101x7 matrix
      # bind fitnesses[genotype1] as the first column of output1
      # the result is a 101x8 matrix with columns: genotype1, genotype2, 'Xd_equil', 'XX', 'XdX', 'XdXd', 'XY', 'XdY' 
      # where genotype1 is equal to i/step for all 101 rows, while genotype2 ranges from 0 to 1 in 0.01 steps (101 rows)
      output1 <- cbind(fitnesses[genotype1], output1) # 101x8 matrix
    } else {
      # step up genotype1, generate a new 101x8 matrix, and bind to output1
      output2 <- fitness_loop(A, generations, fitnesses, n, m, step, genotype2) # 101x7 matrix
      output2 <- cbind(fitnesses[genotype1], output2) # 101x8 matrix
      output1 <- rbind(output1, output2)
    }
  }
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(time.taken)
  return(output1)
}



#--------------------------------------------------------------------------------
# This is a faster implementation of the double_fitness_loop function
# This function runs the model over a range of fitness values for two pre-defined genotypes (nested for loop)
# This is done by calling the fitness_loop function (which calls the run function, which calls the model function)
# The fitnesses of the chosen genotypes each range from 0 to 1 in 0.01 steps
# Each run of the model is for x generations (usually 10,000)
# The frequency of all the genotypes (and the summed Xd frequency) at the end of each run is saved into a row of a matrix
# 
# The function takes the following parameters:
#     A = data matrix to be populated (its first row is initiated with the starting genotype frequencies)
#     generations = the number of generations the model will run for (integer value) 
#     fitnesses = a dataframe containing the fitnesses of each genotype (columns are indexed in the model function)
#     n = the female drive strength (numeric value) 
#     m = the male drive strength (numeric value)
#     step = the number of fitness increments over which to iterate (integer value)
#     genotype1 = the first genotype whose fitness will be changed with each step (and that is defined by i/step)
#     genotype2 = the second genotype whose fitness will be changed with each step (and that is defined by i/step)

nested_fitness_loop <- function(A, generations, fitnesses, n, m, step, genotype1, genotype2){
  start.time <- Sys.time()
  # Pre-allocate output matrix
  output <- matrix(NA, nrow = (step + 1) * (step + 1), ncol = 8)
  colnames(output) <- c(genotype1, genotype2, "Xd_equil", "XX", "XdX", "XdXd", "XY", "XdY")
  # Loop to fill in the output matrix
  k <- 1
  for (i in 0:step){
    fitnesses[genotype1] <- i/step
    col1 <- matrix(rep(fitnesses[genotype1][[1]], 101)) 
    A2 <- fitness_loop(A, generations, fitnesses, n, m, step, genotype2)
    output[k:(k+100), ] <- cbind(col1, A2)
    k <- k + 101
  }
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(time.taken)
  return(output)
}



#---------------------------------------------------------------------------------
# this function iteratively runs the model over a range of drive values for either m or n
drive_loop <- function(A, generations, fitnesses, n, m, step, scenario){
  A2 <- matrix(data = 0, ncol = 2, nrow = step+1)
  if (scenario == 'male') {
    for (i in 0:step){
      m <- i/step
      output <- run(A, generations, fitnesses, n, m)
      df <- as.data.frame(output)
      colnames(df) <- c('XX', 'XdX', 'XdXd', 'XY', 'XdY')
      df$Xd_freq <- (df$XdXd) + (df$XdX/2) + (df$XdY)
      fixed <- df$Xd_freq > 0.99999999
      row_number <- which(fixed)[1]
      A2[i+1,][1] <- i/step 
      A2[i+1,][2] <- row_number
    }
    return(A2)
  } else {
    for (i in 0:step){
      n <- i/step
      output <- run(A, generations, fitnesses, n, m)
      df <- as.data.frame(output)
      colnames(df) <- c('XX', 'XdX', 'XdXd', 'XY', 'XdY')
      df$Xd_freq <- (df$XdXd) + (df$XdX/2) + (df$XdY)
      fixed <- df$Xd_freq > 0.99999999
      row_number <- which(fixed)[1]
      A2[i+1,][1] <- i/step 
      A2[i+1,][2] <- row_number
    }
    return(A2)
  }
}

# fin