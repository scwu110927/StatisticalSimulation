#HW4
#1################################
data <- sunspot.year
train.data <- data[1:279]
test.data <- data[280:289]

#block bootstrap
zeroma <- matrix(0, 10000, ncol = 10)
dr <- train.data[-1] - train.data[-279]
#dr 為前後差異
set.seed(10)
for(i in 1:10000){
  sa <- sample(1:269, 1)
  zeroma[i,] <- dr[c(sa:(sa+9))] 
  #將每十個數存入矩陣row
}
  
md <- apply(zeroma, 2, median)
# 變動量的中位數

predict.value <- NULL
for(i in 1:10){
  predict.value <- c(predict.value, 
                     train.data[279] + sum(md[c(1:i)]))
  #最後一個數加上中位數
}
b <- ts(predict.value, frequency = 1 ,start = c(1979,1))

#AR(2)
x <- ar.ols(train.data, order = 2)
z <- predict(x, n.ahead = 10)
a <- ts(z$pred, frequency = 1 ,start = c(1979,1))

ts.plot(data, a, b, ylab ="Sunspot", col=c("gray47", "red", "blue"))
legend(1700, 195,                                
       lty=1:1,                                   
       col = c("black", "red", "blue"),        
       legend = c("Data", "AR(2)", "Block"), cex = 0.75)
#AR(2) 預測較佳

#2################################
ano.pvalue <- function(mu, n){
  p.vector <- NULL
  power.v <- NULL
  #k <- 0
  for(i in 1:n){
    x <- matrix(rnorm(48), nrow = 3, byrow = T)
    x[2,] <- rnorm(16, mu, 4)
    mui <- apply(x,1,mean)
    mu <- mean(x)
    sstr <- sum(16*(mu-mui)^2)
    ssto <- sum((x-mu)^2)
    f <- (sstr/2)/((ssto-sstr)/45)
    p <- pf(f, 2, 45, lower.tail = F)
    power <- pf(qf(p, 2, 45, lower.tail = F),2,45, lower.tail=F)
    #if(p < 0.05){k <- k+1}
    p.vector <- c(p.vector, p)
    power.v <- c(power.v, power)
  }
  return(cbind(p.vector,power.v))
}
a <- ano.pvalue(1,10000)
a
hist(a[,1])
#p.value不是均勻分配
ks.test(a[,1], "punif")
#?mu要不要算機率
#power?

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

#bisection function
bisection <- function(f, a, b){
  h <- abs(b - a)/10 #切十等分
  i <- 0 #找兩個函數相乘小於零的次數
  j <- 0 #找兩者中間數的次數
  a1 = b1 = 0
  #options(digits = 5)
  while(i < 10){
    a1 = a + i * h #從左邊開始加一等分
    b1 = a1 + h
    if(f(a1) == 0){
      cat('Root:', a1, 'n:', j+1)
    }else if(f(b1) == 0){
      cat('Root:', b1,  'n:', j+1)
    }else if(f(a1) * f(b1) < 0){
      repeat{
        if(abs(b1 - a1) < 1e-5){break} #當兩個值得距離小於1e-10則停止
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
    i <- i+1 #找兩個相乘<0
  }
  if(j == 0){ #表示十等分裡都沒有根
    cat("there is no root between", a, "and", b)
  }
}

bisection(func.a, -4, 3) #0.275684
func.a(0.2756839752)
curve(func.a, xlim = c(-3,3), col='blue', lwd = 1.5, lty = 2)
abline(h=0)
abline(v=0)
polyroot(c(-1, 3, 2, 1)) 
#[1] 0.2756822+0.000000i 
#[2] -1.1378411+1.527312i
#[3] -1.1378411-1.527312i
#bisection 找不到非實數解?
uniroot(func.a, c(-4, 3)) #check r 0.2756814

bisection(func.b, -4, 3) #-3.500085、-3.132141
curve(func.b, xlim = c(-5,2), col='blue', lwd = 1.5, lty = 2)
abline(h=0)
abline(v=0)
uniroot(func.b, c(-4.5, -3.3)) #-3.500085
uniroot(func.b, c(-3.4, -2)) #-3.132141

bisection(func.c, -4, 3) #-1.315119、0.3984478
curve(func.c, xlim = c(-2,2), col='blue', lwd = 1.5, lty = 2)
abline(h=0)
abline(v=0)
uniroot(func.c, c(-2, -1)) #-1.31514
uniroot(func.c, c(1, 0)) #0.3984464

#library(NLRoot)
#BFfzero(func.a, -4, 3)
#BFfzero(func.b, -4, 3)
#BFfzero(func.c, -4, 3)

#convergence criterion?

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
#0.0357123
secant(mitinom.loglike, .01, .99)
#secant(mitinom.loglike, .01, .1)
false.posi(mitinom.loglike, .01, .99)
ridder(mitinom.loglike, .01, .99)
newton.raphson(mitinom.loglike, .01)
#newton.raphson(mitinom.loglike, .1)

#5################################

f=function(x) { x+1/x }
nlminb(start=2,obj=f,lower=0)

