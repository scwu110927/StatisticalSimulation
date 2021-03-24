##-HW1-################
#1(a)
rep(seq(0,4), each = 5)

#1(b)
seq(1:5) + rep(0:4, each = 5)

#1(C) 
x <- c("red","yellow","blue","green","magenta","cyan")
i <- seq(1:3) + rep(0:3, each = 3)
x[i]


#2(a)-################
library(spatstat)

mindist <- function(x, y){
 min(nndist(x, y))
}

x <- runif(20)
y <- runif(20)
mindist(x, y)

#2(b)
grid <- function(x){
  if(x[1] > 0.5 && x[2] > 0.5){index <- 1
  }else if(x[1] < 0.5 && x[2] > 0.5){index <- 2
  }else if(x[1] < 0.5 && x[2] < 0.5){index <- 3
  }else {index <- 4}
}
xy <- rbind(x,y)
index <- apply(xy, 2, grid)

plot(x, y, pch=index)
abline(h = 0.5, v = 0.5, lty=2)

#2(c)
symbols(x, y, circles = index, add = T)


#3-################
gcd = function(a,b){
  if (b==0) a else gcd(b, a %% b) 
}

lcm = function(c,d){
  return(c * d / gcd(c,d)) 
}


#4(a)-################
midsqur <- function(seed,times){
  numvector <- NULL
  for(i in 1:times){
    num <- seed * seed
    seed <- (num%/%1000) %% 1000000
    numvector <- c(numvector, seed)
  }
  numvector <- (numvector / 10^6)
  return(numvector)
}

x <- ceiling(runif(1, 0, 999999)) 
x1 <- midsqur(x, 10000)
hist(x1)
ks.test(x1, y = "punif")

for(i in 1:length(x1)){
  y <- x1[i] - x1
  yc <- x1[which(y == 0)]
}
table(yc)


#4(b)
x1 <- 0.6
x2 <- 0
for (i in 1:10000){
  x1 <- (69069*x1) %% 2^32
  x2 <- c(x2, x1)
}
x2 <- x2[-1] / 2^32
hist(x2)

v <- floor(x2 * 10) 
chisq.test(table(v))
ks.test(x2, y = "punif")


#4(c)
xi <- rnorm(1, 0, 1)
yi <- rnorm(1, 0, 1)
zi <- rnorm(1, 0, 1)
vector <- NULL

for (j in 1:10000) {
  xi <- (171*xi) %% 30269
  yi <- (172*yi) %% 30307
  zi <- (170*zi) %% 30323
  ui <- ((xi/30269) + (yi/30307) + (zi/30323)) %% 1
  vector <- c(vector, ui)
}

v <- floor(x2 * 10) 
chisq.test(table(v))
ks.test(vector, y = "punif")


#5(a)-################
t1 <- NULL
be <- 10000/15
for (i in 1:1000) {
  a1 <- sample(c(1:15), 10000, T)
  a2 <- ceiling(15 * runif(10000))
  b1 <- table(a1)
  b2 <- table(a2)
  c1 <- sum((b1-be)^2 / be)
  c2 <- sum((b2-be)^2 / be)
  d1 <- pchisq(c1, 14)
  d2 <- pchisq(c2, 14)
  t1 <- cbind(t1,c(d1,d2))
}
par(mfrow = c(1, 2))
hist(t1[1,])
hist(t1[2,])

v1 <- table(floor(t1[1,] * 10))
chisq.test(v1)
v2 <- table(floor(t1[2,] * 10))
chisq.test(v2)

#5(b)
casio <- function(seed, times){
  phi <- (1+sqrt(5)) / 2
  uvector <- NULL
  for(i in 1:times){
    seed <- ((phi+seed)^5) %% 1
    uvector <- c(uvector,seed)
  }
  return(uvector)
}

casio2 <- function(seed,times){
  uvector <- NULL
  for(i in 1:times){
    seed <- ((sqrt(2)+seed)^5) %% 1
    uvector <- c(uvector,seed)
  }
  return(uvector)
}

l <- runif(1, 0, 1)
k <- casio(l, 10000)
k2 <- casio2(l, 10000)
par(mfrow = c(1, 2))
hist(k)
hist(k2)

v1 <- floor(k * 10)
chisq.test(table(v1))
v2 <- floor(k2 * 10)
chisq.test(table(v2))

#6
fibonacci <- function(seed, n){ 
  m <- length(seed) - 1 
  for (j in 1:n) { 
    x <- (seed[j] + seed[j+m]) %% 1
    seed <- c(seed,x)
  }
  return(seed[-c(1:(m+1))])
}

 
p <- NULL
pid <- NULL
for(i in 1:10){
  k <- runif(10*i, 0, 1)
  k <- fibonacci(k, 10000)
  pvector <- chisq.test(table(ceiling(k*10)/10))$p.value
  p <- c(p, pvector)
  mat <- matrix(k[-1], ncol=3333, byrow = F)
  mat2 <- apply(mat, 2, rank)
  mat3 <- mat2[1,]*100 + mat2[2,]*10 + mat2[3,]
  pidvector <- chisq.test(table(mat3))$p.value
  pid <- c(pid, pidvector)
}

plot(seq(10, 100, 10), p, type = "b", main = "difference of p-values")
lines(seq(10, 100, 10), pid, type = "b", col = 2)
legend(70, .2, c("GOF", "Indept"), col = c(1, 2), lty = 1, pch = 1)
abline(h = 0.05, col = 3)
