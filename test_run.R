library(MASS)
library(ggplot2)

pop <- c(-0.39, 0.12, 0.94, 1.67, 1.76, 2.44, 3.72, 4.28, 4.92, 5.53, 
         0.06, 0.48, 1.01, 1.68, 1.80, 3.25, 4.12, 4.60, 5.28, 6.22)

data.frame(sample = pop) %>%
    ggplot(., aes(x = sample)) + geom_density()

set.seed(800)

mu1 <- sample(pop, 1, replace=FALSE)
mu2 <- sample(pop, 1, replace=FALSE)
sd1 <- var(pop)
sd2 <- var(pop)
pi <- 0.5

init <- c(mu1, sd1, mu2, sd2, pi)

get_params <- function(x, params){
  
  mu1 <- params[1]
  sd1 <- params[2]
  mu2 <- params[3]
  sd2 <- params[4]
  pi <-  params[5]
  
  initial_pop1 <- sapply(1:length(x), function(y) pnorm(x[y], mu1, sd1))
  initial_pop2 <- sapply(1:length(x), function(y) pnorm(x[y], mu2, sd2))
  
  resp <- NULL
  
  for (i in 1:length(x)){
    
    resp[i] <- pi * initial_pop2[i] / ((1 - pi) * initial_pop1[i] + pi * initial_pop2[i])
    
  }
  
  mu1_1 <- sum((1 - resp) * x) / sum(1 - resp) 
  sd1_1 <- sum((1 - resp) * (x - mu1_1) ^ 2) / sum(1 - resp)
  mu2_1 <- sum(resp  * x) / sum(resp)
  sd2_1 <- sum(resp * (x - mu1_1) ^ 2) / sum(resp)
  pi_1 <- sum(resp) / length(x)
  
  init_1 <- c(mu1_1, sd1_1, mu2_1, sd2_1, pi_1)
  
   return(init_1)
  
}

result <- get_params(pop, init)

results <- list()

for (i in 1:20){
  
  if (i == 1){
    results[[i]] <- get_params(pop, init)
  }
  else{
    results[[i]] <- get_params(pop, results[[i - 1]])
  }
}

results