library(MASS)
library(ggplot2)

#### Reading in the data

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


row_sample_1 <- sample(nrow(data),size=1,replace=TRUE)
row_sample_2 <- sample(nrow(data),size=1,replace=TRUE)

m1 <- data[row_sample_1, ]
m2 <- data[row_sample_2, ]
sd1 <- var(data)
sd2 <- var(data)
pi <- 0.5

init <- c(m1, m2, sd1, sd2, pi)

## avoid spurious accuracy
op <- options(digits = 3)

#### Psuedo code
#### 1. Use pi to randomly determine which dist a given sample is from
#### 2. Update mu and sd based on upated sample
#### 3. Use sample to generate distribution from MASS package, save MLE
#### 3.5 This will be done twice for each adjustment of pi as there are two distributions
#### Do we take the higher of the MLE and go with that? For instance if a decrease in pi drove ML up then logic would dictate you
#### should keep doing it until you get a result
#### 4. Iterate until convergence of pi and mu/sd
#### 5. After convergence return the parameters for each of the distributions

mixture_model <- function(initial_params, X){
  
  ### initial_params : a list of variables holding first estimates of mu, sd and pi
  ### data : the population from which the underlying subpopulations are being investigated
  ### returns : TBC
  
  pop <- as.matrix(X)
  rows <- nrow(pop)
  
  elements_1 <- pi * rows
  elements_2 <- (1 - pi) * rows
  
  pop1 <- sample(pop, elements_1, replace=FALSE)
  pop2 <- setdiff(pop, pop1)
      
}

list[pop1, pop2] <- mixture_model(init, data)
length(pop)
### initial_params : a list of variables holding first estimates of mu, sd and pi
### data : the population from which the underlying subpopulations are being investigated
### returns : TBC

pop <- as.matrix(data)
rows <- nrow(pop)

elements_1 <- pi * rows
elements_2 <- (1 - pi) * rows

pop1 <- sample(pop, elements_1, replace=FALSE)
pop2 <- pop[sample(pop, elements_2, replace=TRUE), ]

  
head(pop)


nrow(pop)

sample(data, pi * nrow(data), replace=TRUE)