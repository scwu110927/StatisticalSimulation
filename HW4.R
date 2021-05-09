#HW4
#1################################
data <- sunspot.year
train.data <- data[1:279]
test.data <- data[280:289]


zeroma <- matrix(0, 10000, ncol = 10)
dr <- train.data[-1] - train.data[-279]

set.seed(10)
for(i in 1:10000){
  sa <- sample(1:269, 1)
  zeroma[i,] <- dr[c(sa:(sa+9))] 
}
md <- apply(zeroma, 2, median)
predict.value <- NULL
for(i in 1:10){
  predict.value <- c(predict.value, 
                     train.data[279] + sum(md[c(1:i)]))
}
b <- ts(predict.value, frequency = 1 ,start = c(1979,1))

x <- ar.ols(train.data, order = 2)
z <- predict(x, n.ahead = 10)
a <- ts(z$pred, frequency = 1 ,start = c(1979,1))

ts.plot(data, a, b, ylab ="Sunspot", col=c("gray47", "red", "blue"))
legend(1700, 195,                                
       lty=1:1,                                   
       col = c("black", "red", "blue"),        
       legend = c("Data", "AR(2)", "Block"), cex = 0.75)


#2################################
ano.pvalue <- function(mu=0, sigma=1, n){
  p.vector <- NULL
  power.v <- NULL
  #k <- 0
  for(i in 1:n){
    x <- matrix(rnorm(48), nrow = 3, byrow = T)
    x[2,] <- rnorm(16, mu, sigma)
    mui <- apply(x,1,mean)
    mu <- mean(x)
    mstr <- sum(16*(mu-mui)^2)/2
    mse <- (sum((x-mu)^2)-mstr*2)/45
    f <- mstr/mse
    p <- pf(f, 2, 45, lower.tail = F)
    power <- ?power.anova.test(groups = 3, n = 16, between.var = mstr, 
                              within.var = mse, sig.level = p)$power
    p.vector <- c(p.vector, p)
    power.v <- c(power.v, power)
  }
  return(cbind(p.vector,power.v))
}
a <- ano.pvalue(n=1000)
b <- ano.pvalue(sigma=sqrt(2),n=1000)

hist(a[,1])
hist(a[,2])
mean(a[,2])
hist(b[,1])
ks.test(a[,1], "punif")
ks.test(b[,1], "punif")


ano.pvalue <- function(mu=0, sigma=1, n){
  p.vector <- NULL
  power.v <- NULL
  #k <- 0
  for(i in 1:n){
    x <- cbind(rep(1:3, each = 16),
               c(rnorm(32), rnorm(16, mu, sigma)))
    fit <- anova(lm(x[, 2] ~ factor(x[, 1])))
    ms <- fit$"Mean Sq"
    p <- fit$"Pr(>F)"[1]
    power <- power.anova.test(groups = 3, n = 16, between.var = ms[1], 
                              within.var = ms[2])$power
    p.vector <- c(p.vector, p)
    power.v <- c(power.v, power)
  }
  return(cbind(p.vector,power.v))
}


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
        if(abs(b1 - a1) < 1e-5){break} 
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

bisection(func.a, -4, 3) 
curve(func.a, xlim = c(-3,3), lwd = 1.5, lty = 2)
abline(h=0)
polyroot(c(-1, 3, 2, 1)) 
uniroot(func.a, c(-4, 3)) 
false.posi(func.a, -4, 3)

bisection(func.b, -4, 3) 
bisection(func.b, -3.4, 3) 
curve(func.b, xlim = c(-5,2), lwd = 1.5, lty = 2)
abline(h=0)
uniroot(func.b, c(-3.4, -2)) 
false.posi(func.b, -3.4, -2)

bisection(func.c, -4, 3) 
curve(func.c, xlim = c(-2,2), lwd = 1.5, lty = 2)
abline(h=0)
uniroot(func.c, c(-2, -1)) 
uniroot(func.c, c(1, 0)) 
false.posi(func.c, -2, -1)
false.posi(func.c, 1, 0)


#4################################
mitinom.loglike <- function(theta){
  1997/(2+theta) - 906/(1-theta) - 904/(1-theta) + 32/(theta)
}
curve(mitinom.loglike, from=-0.3, to=0.3)
abline(h=0)

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
#secant(mitinom.loglike, .01, .1)
ridder(mitinom.loglike, .01, .99)
newton.raphson(mitinom.loglike, .01)
#newton.raphson(mitinom.loglike, .1)

#5################################

f=function(x) { x+1/x }
nlminb(start=2,obj=f,lower=0)

