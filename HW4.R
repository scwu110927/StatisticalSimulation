#HW4
#1################################
data <- sunspot.year
train.data <- data[9:279]
test.data <- data[280:289]

zeroma <- matrix(0, 10000, ncol = 10)
dr <- train.data[-1] - train.data[-271]

for(i in 1:10000){
  sa <- sample(1:27, 1)*10-9
  zeroma[i,] <- dr[c(sa:(sa+9))] 
}
md <- apply(zeroma, 2, median)
predict.value <- NULL
for(i in 1:10){
  predict.value <- c(predict.value, data[279] + sum(md[c(1:i)]))
}
b <- ts(predict.value, frequency = 1 ,start = c(1979,1))

x <- ar.ols(train.data, order = 2)
z <- predict(x, n.ahead = 10)
a <- ts(z$pred, frequency = 1 ,start = c(1979,1))

ts.plot(data, a, b, ylab ="Sunspot", lty = 1:3, 
        main = "Prediction last 10 year of Sunspot with Block Bootstrap & AR(2)")
legend(1700, 195, lty = 1:3, legend = c("Data", "AR(2)", "Block"), cex = 0.75)


#2################################
type1.test <- function(n){
  type1.vector <- NULL
  for(i in 1:n){
    x <- cbind(rep(1:3, each = 16), rnorm(48))
    fit <- lm(x[, 2] ~ factor(x[, 1]))
    type1 <- anova(fit)$"Pr(>F)"[1]
    type1.vector <- c(type1.vector, type1)
  }
  return(type1.vector)
}

type1 <- type1.test(n=1000)
hist(type1, main = "Type I Error with all zero means in simu 1000 times")
ks.test(type1, 'punif')

power.test <- function(sigma = 1, n){
  power.vector <- NULL
  mu <- seq(0, 2, .2)
  for(i in 1:length(mu)){
    power.sum <- 0
    for(j in 1:n){
    x <- cbind(rep(1:3, each = 16), c(rnorm(32), rnorm(16, mu[i], sigma)))
    fit <- lm(x[, 2] ~ factor(x[, 1]))
    F.value <- anova(fit)$'F value'[1]
    power <- pf(qf(.95, 2, 45), 2, 45, ncp=F.value*2, lower.tail = F)
    power.sum <- power.sum + power
    }
    power.vector <- c(power.vector, mean(power.sum))
  }
  return(cbind(mu, power.vector))
}

constant.var <- power.test(n = 1000)
nonconstant.var <- power.test(sigma = sqrt(2), n = 1000)
plot(constant.var, type = 'l', 
     main = "Compare Constant and Non-Constant Variance with Power in simu 1000 times")
lines(nonconstant.var, lty = 2)
legend('topleft', lty = 1:2, legend = c("Constant variance", "Non-Constant variance"))


#3################################
func.a <- function(x){
  x^3 + 2*(x^2) + 3*x - 1
}
func.b <- function(x){
  exp(1) - 1/(3.5 + x)
}
func.c <- function(x){
  exp(-x)/((1 + x^2)^2) - 0.5
}

bisection <- function(f, a, b){
  h <- abs(b - a)/10 
  i <- 0 
  j <- 0 
  a1 = b1 = 0
  while(i < 10){
    a1 = a + i * h 
    b1 = a1 + h
    if(f(a1) == 0){
      cat('Root:', a1, 'n:', j+1)
    }else if(f(b1) == 0){
      cat('Root:', b1,  'n:', j+1)
    }else if(f(a1) * f(b1) < 0){
      repeat{
        if(abs(b1 - a1) < 1e-6){break} 
          c <- (a1 + b1)/2
        if(f(a1) * f(c) < 0){
          b1 <- c
          j <- j+1
        }else{ 
          a1 <- c
          j <- j+1
        }
      }
      c <- (a1 + b1)/2
      cat('Root:', c, 'n:', j+1, 'f:', f(c), '\n')
    }
    i <- i+1 
  }
  if(j == 0){ 
    cat("there is no root between", a, "and", b)
  }
}

false.posi <- function(f, init1, init2, maxiter=1000, tol=1e-06) {
  init1[2] <- f(init1)
  init2[2] <- f(init2)
  if (init1[2] == 0.0) {return(init1)}
  if (init2[2] == 0.0) {return(init2)}
  for (i in 1:maxiter) {
    dummy <- init2[1]
    init2[1] <- init1[1]-init1[2]*(init1[1]-init2[1])/(init1[2]-init2[2])
    init2[2] <- f(init2[1])
    if (abs(init2[1]-dummy) < tol) break
  }
  list(root = init2[1], f.root = init2[2], iter = i)
}

curve(func.a, xlim = c(-3,3), lwd = 2)
abline(h=0, lty = 2)
points(0.275, 0, cex = 2, pch = 21, bg = 1)
bisection(func.a, -4, 3) 
false.posi(func.a, -4, 3)
uniroot(func.a, c(-4, 3)) 

curve(func.b, xlim = c(-5,-2), lwd = 2)
abline(h=0, lty = 2)
points(-3.132, 0, cex = 2, pch = 21, bg = 1)
bisection(func.b, -3.4, -2) 
false.posi(func.b, -3.4, -2)
uniroot(func.b, c(-3.4, -2))

curve(func.c, xlim = c(-2,2), lwd = 2)
abline(h=0, lty = 2)
points(-1.315, 0, cex = 2, pch = 21, bg = 1)
points(0.398, 0, cex = 2, pch = 21, bg = 1)
bisection(func.c, -2, 1) 
false.posi(func.c, -2, -1)
false.posi(func.c, 0, 1)
uniroot(func.c, c(-2, -1)) 
uniroot(func.c, c(0, 1)) 



#4################################
mitinom.loglike <- function(theta){
  1997/(2+theta) - 906/(1-theta) - 904/(1-theta) + 32/(theta)
}
curve(mitinom.loglike, from=-0.3, to=0.3, lwd = 2)
abline(h=0, lty = 2)
points(0.03571232, 0, cex = 2, pch = 21, bg = 1)

secant <- function(f, init1, init2, maxiter=1000, tol=1e-06){
  init1[2] <- f(init1)
  init2[2] <- f(init2)
  if (init1[2] == 0.0) {return(init1)}
  if (init2[2] == 0.0) {return(init2)}
  for (i in 1:maxiter) {
    dummy <- init2
    init2[1] <- init1[1]-init1[2]*(init1[1]-init2[1])/(init1[2]-init2[2])
    init2[2] <- f(init2[1])
    init1 <- dummy
    if (abs(init2[1]-init1[1]) < tol) break
  }
  list(root = init2[1], f.root = init2[2], iter = i)
}

ridder <- function(f, init1, init2, maxiter=1000, tol=1e-06){
  init1[2] <- f(init1)
  init2[2] <- f(init2)
  if (init1[2] == 0.0) {return(init1)}
  if (init2[2] == 0.0) {return(init2)}
  if (sign(init1[2]*init2[2]) != -1)
    stop("f() values at end points not of opposite sign")
  for (i in 1:maxiter) {
    init3 <- (init1[1] + init2[1])/2
    init3[2] <- f(init3)
    init4 <- init3[1] + (init3[1]-init1[1])*sign(init1[2])*init3[2]/
      sqrt(init3[2]^2-init1[2]*init2[2])
    init4[2] <- f(init4)
    if (sign(init3[2]*init4[2]) == -1){init1 <- init3}
    else if (sign(init1[2]*init4[2]) == -1){init1 <- init1}
    else if (sign(init2[2]*init4[2]) == -1){init1 <- init2}
    else {stop("f() values at iteration points not of opposite sign")}
    init2 <- init4
    if (abs(init2[1]-init1[1]) < tol) break
  }
  list(root = init2[1], f.root = init2[2], iter = i)
}

newton.raphson <- function(f, init, maxiter=1000, tol=1e-06){
  require(numDeriv) 
  if (f(init) == 0.0) {return(init)}
  for (i in 1:maxiter) {
    dx <- genD(func = f, x = init)$D[1] 
    init.next <- init - f(init)/dx
    if (abs(init.next - init) < tol) break
    init <- init.next
  }
  list(root = init.next, f.root = f(init.next), iter = i)
}

uniroot(mitinom.loglike, c(0.01, 0.99), tol = 1e-06)
secant(mitinom.loglike, .01, .99)
ridder(mitinom.loglike, .01, .99)
newton.raphson(mitinom.loglike, .01)


#5################################

wide <- read.csv("StatisticalSimulation/maledeathrates.csv", header = T)
long_p <- reshape(wide, direction = "long", 
                       varying = list(names(wide)[2:4]),
                       v.names = "p", idvar = "year",
                       timevar = "t", times = 1:3)
long_n <- reshape(wide, direction = "long", 
                  varying = list(names(wide)[5:7]),
                  v.names = "n", idvar = "year",
                  timevar = "t", times = 1:3)
long_d <- reshape(wide, direction = "long", 
                  varying = list(names(wide)[8:10]),
                  v.names = "d", idvar = "year",
                  timevar = "t", times = 1:3)
long <- cbind(long_p, long_n, long_d)[, c(1, 8, 9, 18, 27)]
long <- long[order(long$year),]
write.table(long, file = "StatisticalSimulation/maledeathrates2.CSV", 
            sep=",", row.names = F, na = "NA")

long <- read.csv("StatisticalSimulation/maledeathrates2.csv", header = T)
y <- log(-log(1-k[-86, 2]))
y2 <- 1-k[-86, 2]
d <- k[-86, 4]
x <- k[-86, 1]
w <- k[-86, 3]


#weights nx
g <- lm(y ~ x, weights = k[-86, 3])
summary(g)
beta <- g$coefficients[2]
alpha <- g$coefficients[1]
c <- exp(beta)
b <- exp(alpha + log(log(c)) - log(c - 1))
b*c

#weights sqrt(nx)
g1 <- lm(y~x, weights = sqrt(as.numeric(k[2:86, 3])))
summary(g1)
beta_1 <- g1$coefficients[2]
alpha_1 <- g1$coefficients[1]
c_1 <- exp(beta_1)
b_1 <- exp(alpha_1 + log(log(c_1)) - log(c_1 - 1))
b_1*c_1


wls <- function(par) { 
  x1 <- par[1]
  x2 <- par[2] 
  sum(w*(y-x1-x2*x)^2) }
nlminb(start=c(-10, 1), obj=f)

nlm <- function(par) { 
  x1 <- par[1]
  x2 <- par[2] 
  sum(w*(y2-exp(-x1*(x2^x)*(x2-1)/log(x2)))^2) }
nlminb(start=c(0.5, 2), obj=f)

mle <- function(par) { 
  x1 <- par[1]
  x2 <- par[2] 
  sum((w-d)*x1*x2^x*(x2-1)/log(x2)-d*log(1-exp(-x1*(x2^x)*(x2-1)/log(x2)))) }
nlminb(start=c(0.5, 2), obj=f)


