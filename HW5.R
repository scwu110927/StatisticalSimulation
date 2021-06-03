#HW5
#1################################


#2################################


#3################################
# Antithetic Variate
z1 <- NULL
z2 <- NULL
k <- 0
rho <- 0
n <- 1000
for (i in 1:n){
rn <- rnorm(100, 0, 2 + rho)
z1 <- c(z1, sum(rn < k)/100)
z2 <- c(z2, sum(-rn < k)/100)
}
c(mean(c(z1, z2)), var(c(z1, z2)))

# Stratified Sampling
library('truncnorm')
n <- 1000
z1 <- rtruncnorm(n * .2, b = qnorm(0.2))
z1 <- rtruncnorm(n * .2, b = qnorm(0.2))
z1 <- rtruncnorm(n * .2, b = qnorm(0.2))
z1 <- rtruncnorm(n * .2, b = qnorm(0.2))


