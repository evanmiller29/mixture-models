library(MASS)
library(ggplot2)
library(dplyr)

#### Reading in the data

op <- options(digits = 7)

rm(list = ls())

basepath <- 'C:/Users/evanm_000/Documents/GitHub/mixture-models'
setwd(basepath)

data <- read.table("EvanExercise.dat", header=FALSE)
colnames(data) <- 'population'

#### Getting an overview of what the data looks like

ggplot(data, aes(x = population)) + geom_density() + xlim(0, 1.1) +theme_minimal() + ggtitle('Overview of population')

set.seed(737)
pop = as.matrix(data)

#### Starting points are as below. Will use K means to generate more reasonable estimates in future
#### Reducing the size of the starting points to make them more feasible

alpha1 <- 0.55
beta1 <- 0.124
alpha2 <- 0.5
beta2 <- 0.25
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
  class_member <- NULL
  initial_pop1 <- NULL
  initial_pop2 <- NULL
  
  for (i in 1:length(x)){
    
    initial_pop1[i] <- pbeta(x[i], alpha1, beta1, ncp = 0, lower.tail = TRUE, log.p = FALSE)
    initial_pop2[i] <- pbeta(x[i], alpha2, beta2, ncp = 0, lower.tail = TRUE, log.p = FALSE)
      
    resp[i] = (pi * initial_pop2[i]) / ((1 - pi) * initial_pop1[i] + pi * initial_pop2[i])

    #### Now want to add bernoulli random draws to split into two classes. Split determined on pi
    
    class_member[i] <- rbinom(n = 1, size = 1, prob = pi)
    
    #### Iterate over these draws 8000-10000 times to get median parameters for next iteration
    #### Keep going until we get convergence
  }
  
  #### Do I fit the beta distributions based on the observations or the responsiveness measure for the classes?
    
  results <- data.frame(obs = x, resp = resp, class = class_member)
  dists <- fit_beta(results, init_params)
  
  return(dists)
}

fit_beta <- function(df, params){
    
  pop1 <- as.matrix(df$population[df$class == 1])
  pop2 <- as.matrix(df$population[df$class == 0])
  
  beta_dist1 <- fitdistr(pop1, dbeta, list(shape1 = params[1], shape2 = params[2]))
  beta_dist2 <- fitdistr(pop2, dbeta, list(shape1 = params[3], shape2 = params[4]))
  
  dists <- c(beta_dist1, beta_dist2)
  
  return(dists)
  
}

init_dists <- comp_resp(pop, init)
init_dists[1]

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

## avoid spurious accuracy
op <- options(digits = 3)
set.seed(123)
x <- rgamma(100, shape = 5, rate = 0.1)
fitdistr(x, "gamma")


x <- rbeta(100, shape1 = 5, shape2 = 10)
fitdistr(x, dbeta, list(shape1 = 1, shape2 = 2))

## now do this directly with more control.
fitdistr(x, dgamma, list(shape = 1, rate = 0.1), lower = 0.001)

set.seed(123)
x2 <- rt(250, df = 9)
fitdistr(x2, "t", df = 9)
## allow df to vary: not a very good idea!
fitdistr(x2, "t")
## now do fixed-df fit directly with more control.
mydt <- function(x, m, s, df) dt((x-m)/s, df)/s
fitdistr(x2, mydt, list(m = 0, s = 1), df = 9, lower = c(-Inf, 0))
set.seed(123)
x3 <- rweibull(100, shape = 4, scale = 100)
fitdistr(x3, "weibull")
set.seed(123)
x4 <- rnegbin(500, mu = 5, theta = 4)
fitdistr(x4, "Negative Binomial")
options(op)

