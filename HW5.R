#HW5
#1################################
#Impartance_sampling##########
importance_sampling <- function(k){
  t1 <- NULL
  for (i in 1:1000) {
    x <- runif(1000,0,k)
    y <- k/(pi*(1+x^2))
    z <- 0.5-y
    a <- mean(z)
    t1 <- c(t1,a)
  }
  return(cbind(mean(t1),var(t1)))
}

#conrtol variate###########
control <- function(k){
  t1 <- NULL
  for(i in 1:1000){
    x <- runif(1000,0,k)
    y <- x^2
    z <- x^4
    fx <- k/(pi*(1+x^2))
    g <- lm(fx~y+z)
    a1 <- g$coefficients[2]*(y-mean(y))
    a2 <- g$coefficients[3]*(z-mean(z))
    a <- 0.5 - fx + a1 + a2
    t1 <- c(t1, a)
  }
  return(cbind(mean(t1), var(t1)))
}

#antithetic#############
antithetic <- function(k){
  t1 <- NULL
  for (i in 1:1000) {
    a <- runif(1000)
    b <- 1-a
    x <- qcauchy(a)
    x2 <- qcauchy(b)
    y <- sum(x>k)/1000
    y2 <- sum(x2>k)/1000
    t1 <- c(t1,(y+y2)/2)
  }
  return(cbind(mean(t1),var(t1)))
}

#result#####
t2 <- NULL
t3 <- NULL
t4 <- NULL
for(i in 1:7){
  b <- c(6, 5, 4, 3.5, 3, 2.5, 2)
  k2 <- importance_sampling(b[i])
  k3 <- control(b[i])
  k4 <- antithetic(b[i])
  t2 <- rbind(t2,k2)
  t3 <- rbind(t3,k3)
  t4 <- rbind(t4,k4)
}
t2 #importance_sampling result
t3 #control result
t4 #antithetic result

#2################################
x <- matrix(rexp(10000), byrow = T, ncol = 5)
y <- x[,1]+2*x[,2]+3*x[,3]+4*x[,4]+5*x[,5]
sum(y >= 21.6)/2000
a <- y >= 21.6
var(a)
#第一個方式
#直接生成2000組5筆服從指數分配(1)的亂數
#看>=21.6 的機率

#antithetic#############
antithetic <- function(k){
  t1 <- NULL
  for (i in 1:k) {
    a <- runif(1000*5)
    b <- 1-a
    x <- qexp(a)
    x2 <- qexp(b)
    y1 <- matrix(x, byrow = T, ncol = 5)
    y2 <- matrix(x2, byrow = T, ncol = 5)
    z1 <- sum(y1[,1]+2*y1[,2]+3*y1[,3]+4*y1[,4]+5*y1[,5]>=21.6)/1000
    z2 <- sum(y2[,1]+2*y2[,2]+3*y2[,3]+4*y2[,4]+5*y2[,5]>=21.6)/1000
    t1 <- c(t1,(z1+z2)/2)
  }
  return(cbind(mean(t1),var(t1)))
}

antithetic(1000)

#a#########
x1 <- matrix(runif(10000), byrow = T, ncol = 5)
y1 <- -log(x1)/1
z1 <- y1[,1]+2*y1[,2]+3*y1[,3]+4*y1[,4]+5*y1[,5]
sum(z1 >= 21.6)/2000
a1 <- z1 >= 21.6
var(a1)
#先生成10000筆隨機均勻亂數
#透過-log轉換使之變指數分配
#然後看>=21.6的機率為何

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
var.reduct(k = 1, rho = -0.5)
var.reduct(k = 2, rho = 0)
var.reduct(k = 3, rho = 0.3)
var.reduct(k = 4, rho = 0.7)


#4################################
library(LaplacesDemon)
rn <- rnormm(100, p = c(0.5, 0.5), mu = c(-2, 2), sigma = c(1, 1))
hist(rn, freq = FALSE, main = "Hist of Mixture Normal")
x <- seq(-5, 5, 0.01)
lines(x, dnormm(x, p = c(0.5, 0.5), mu = c(-2, 2), sigma = c(1, 1)))
lines(density(rn2, kernel = "gaussian", bw = 1), lty = 2)
lines(density(rn2, kernel = "rectangular", bw = 1), lty = 3)
lines(density(rn2, kernel = "triangular", bw = 1), lty = 4)
legend('topright', legend = c("Theoretical", "gaussian", "rectangular", 
                             "triangular"), lty = 1:4)

#5################################
library(zoo)
wide <- read.csv("StatisticalSimulation/maledeathrates.csv", header = T)
plot(wide$year, wide$p_2019, type = 'l', xlim = c(90, 100), 
     main = "Age-specific Mortality Rates")
lines(wide$year[2:100], rollmean(wide$p_2019, k=3), lty = 2)
lines(ksmooth(wide$year, wide$p_2019, "normal", bandwidth = 2), lty = 3)
lines(smooth.spline(wide$year, wide$p_2019, df = 10), lty = 4)
legend('bottomright', legend = c("Data Density", "Rollmean(k=3)", "Kernel(bw=2)", 
                              "Spline(df=10)"), lty = 1:4)


#6################################
library(MCMCpack)
bikes <- read.csv("StatisticalSimulation/bikes.csv", h = T)
summary(lm(riders_registered~temp_feel, data = bikes))
posterior <- MCMCregress(riders_registered~temp_feel, data = bikes)
summary(posterior)
windows()
plot(posterior)



