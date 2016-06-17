library(MASS)
library(ggplot2)
library(dplyr)

#### Reading in the data

op <- options(digits = 3)

rm(list = ls())

basepath <- 'C:/Users/evanm_000/Documents/GitHub/mixture-models'
setwd(basepath)

data <- read.table("EvanExercise.dat", header=FALSE)
colnames(data) <- 'population'

#### Getting an overview of what the data looks like

ggplot(data, aes(x = population)) + geom_density() + xlim(0, 1.1) +theme_minimal() + ggtitle('Overview of population')

#### Setting up the initial values for the Mixture Model

set.seed(737)

#### Initial guesses will just be:
#### mu = random element
#### sd = population variance
#### mixing proportion = 0.5

pop = as.matrix(data)

#### Starting points are as below. Will use K means to generate more reasonable estimates in future

alpha1 <- 5
beta1 <- 1
alpha2 <- 10
beta2 <- 5
pi <- 0.5

init <- c(alpha1, beta1, alpha2, beta2, pi)

initial <- pop[1]

comp_resp <- function(x, init_params){
  
  ### x: the observations for which responsibilties will be calculated
  ### init_params: the initial values for the beta distributions
  ###             1 = alpha parameter for first dist
  ###             2 = beta parameter for first dist
  ###             3 = alpha parameter for second dist
  ###             4 = beta parameter for second dist
  ###             5 = the class split
  
  alpha1 <- init_params[1]
  beta1 <- init_params[2]
  alpha2 <- init_params[3]
  beta2 <- init_params[4]
  pi <- init_params[5]

  resp <- NULL
  
  for (i in 1:length(x)){
    
    initial_pop1[i] <- pbeta(x[i], alpha1, beta1, ncp = 0, lower.tail = TRUE, log.p = FALSE)
    initial_pop2[i] <- pbeta(x[i], alpha2, beta2, ncp = 0, lower.tail = TRUE, log.p = FALSE)
      
    resp[i] = (pi * initial_pop2[i]) / ((1 - pi) * initial_pop1[i] + pi * initial_pop2[i])
    
    #### Now want to add bernoulli random draws to split into two classes
    #### Iterate over these draws 8000-10000 times to get median parameters for next iteration
    #### Keep going until we get convergence
  }
  
  return(resp)
}

first_results <- comp_resp(pop, init)

data.frame(resp = first_results) %>%
    ggplot(., aes(x = resp)) + geom_density()

#### Psuedo code
#### 1. Use pi to randomly determine which dist a given sample is from
#### 2. Update mu and sd based on upated sample
#### 3. Use sample to generate distribution from MASS package, save MLE
#### 3.5 This will be done twice for each adjustment of pi as there are two distributions
#### Do we take the higher of the MLE and go with that? For instance if a decrease in pi drove ML up then logic would dictate you
#### should keep doing it until you get a result
#### 4. Iterate until convergence of pi and mu/sd
#### 5. After convergence return the parameters for each of the distributions
