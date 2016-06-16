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

m1 <- data[row_sample1, ]
m2 <- data[row_sample2, ]
sd1 <- var(data)
sd2 <- var(data)
pi <- 0.5

## avoid spurious accuracy
op <- options(digits = 3)
set.seed(123)
x <- rgamma(100, shape = 5, rate = 0.1)

pop <- data.frame(pop = x)
ggplot(pop, aes(x = pop)) + geom_density() + xlim(0, 150)

fitdistr(x, "gamma")

# now do this directly with more control.
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
