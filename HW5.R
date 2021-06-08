#HW5
#1################################


#2################################


#3################################
library('truncnorm')

var.reduct <- function(k = 0, rho = 0, n1 = 1000, n2 = 100){
  z1 <- NULL
  z2 <- NULL
  z3 <- NULL
  z4 <- NULL
  for (i in 1:n1){
    # Opposite
    rn1 <- rnorm(n2, 0, sqrt(2 + rho))
    z1 <- c(z1, sum(rn1 < k)/n2)
    z2 <- c(z2, sum(-rn1 < k)/n2)
    # Stratified
    rn2 <- NULL
    for (j in 1:5){
      rn2 <- c(rn2, rtruncnorm(0.2*n2, a = qnorm((j-1)/5, sd = sqrt(2 + rho)), 
                               b = qnorm(j/5, sd = sqrt(2 + rho)), 
                               mean = 0,  sd = sqrt(2 + rho)))
    }
    z3 <- c(z3, sum(rn2 < k)/n2)
    # Control
    index <- runif(n2, 0, k)
    index2 <- index^2
    rn3 <- k*dnorm(index, 0, sqrt(2 + rho))
    b <- cov(rn3, index)/var(index)
    b[is.na(b)] <- 0
    z4 <- c(z4, 0.5 + rn3 - b*(index2-mean(index2)))
  }
  return(list("Antithetic mean and variance" = c(mean(c(z1, z2)), var(c(z1, z2))), 
              "Stratified mean and variance" = c(mean(z3), var(z3)),
              "Control mean and variance" = c(mean(z4), var(z4)),
              "Theoretical mean" = pnorm(k, 0, sqrt(2 + rho))))
}

var.reduct(k = 0, rho = -0.9)
var.reduct(k = 0, rho = -0.5)
var.reduct(k = 0, rho = 0)
var.reduct(k = 0, rho = 0.3)
var.reduct(k = 0, rho = 0.7)
var.reduct(k = 1, rho = -0.9)
var.reduct(k = 1, rho = -0.5)
var.reduct(k = 1, rho = 0)
var.reduct(k = 1, rho = 0.3)
var.reduct(k = 1, rho = 0.7)
var.reduct(k = 2, rho = -0.9)
var.reduct(k = 2, rho = -0.5)
var.reduct(k = 2, rho = 0)
var.reduct(k = 2, rho = 0.3)
var.reduct(k = 2, rho = 0.7)
var.reduct(k = 3, rho = -0.9)
var.reduct(k = 3, rho = -0.5)
var.reduct(k = 3, rho = 0)
var.reduct(k = 3, rho = 0.3)
var.reduct(k = 3, rho = 0.7)
var.reduct(k = 4, rho = -0.9)
var.reduct(k = 4, rho = -0.5)
var.reduct(k = 4, rho = 0)
var.reduct(k = 4, rho = 0.3)
var.reduct(k = 4, rho = 0.7)


#4################################




#5################################
library(zoo)
wide <- read.csv("StatisticalSimulation/maledeathrates.csv", header = T)
plot(wide$year, wide$p_2019, type = 'l', xlim = c(80, 100))
lines(wide$year[2:100], rollmean(wide$p_2019, k=3), col = 2)
lines(ksmooth(wide$year, wide$p_2019, "normal", bandwidth = 2), col = 3)
lines(smooth.spline(wide$year, wide$p_2019, df = 10), col = 4)

#6################################
library(MCMCpack)
bikes <- read.csv("StatisticalSimulation/bikes.csv", h = T)
summary(lm(riders_registered~temp_feel, data = bikes))
posterior <- MCMCregress(riders_registered~temp_feel, data = bikes)
summary(posterior)
plot(posterior)



