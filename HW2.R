##-HW2-################
#1################################
#(a)
gap.test <- function(data, a, b){
  data <- data/max(data)
  n <- length(data)
  x <- c(1:n) * (a < data & data < b)
  x1 <- x[x > 0]
  y <- x1[-1] - x1[-length(x1)]-1
  t.y <- table(y)
  e <- (b-a) * (1-(b-a))^(c(1:dim(t.y))-1) * sum(t.y)
  e.p <- e / sum(e)
  c.t <- chisq.test(t.y, p = e.p)
  return(list("gaps" = t.y, "expected count" = round(e), "chisq.test" = c.t))
}

permutation.test <- function(data){
  cond <- length(data) %% 3
  if(cond == 0){
    x1 <- matrix(data, ncol = length(data) %/% 3, byrow = F)
  }else{
    x1 <- matrix(data[-c(1:cond)], ncol = length(data) %/% 3, byrow = F)
  }
  y1 <- apply(x1, 2, rank)
  y2 <- y1[1,]*100 + y1[2,]*10 + y1[3,]
  c.t <- chisq.test(table(y2))
  return(c.t)
}

rn.excel <- as.matrix(read.csv("StatisticalSimulation/rn10000.csv", header = F))
set.seed(2)
rn.r <- runif(10000)
gap.test(rn.excel, 0.2, 0.7)
permutation.test(rn.excel)

gap.test(rn.r, 0.2, 0.7)
permutation.test(rn.r)


#(b)
updown.test <- function(num,runs){
  k.vector <- NULL
  n.r <- 0
  r <- 0
  for(i in 1:runs){
    x <- runif(num)
    x1 <- (x[-1]>x[-num])
    x2 <- sum((x1[-1] != x1[-(num-1)]))
    z <- (x2 - (2*num-1)/3) / sqrt((16*num-29)/90)
    k <- pnorm(z)
    k.vector <- c(k.vector,k)
    if(k > 0.025 & k < 0.975){
      n.r <- n.r + 1
    }else{
      r <- r + 1
    }
  }
  cat(paste("Numbers of not reject (alpha=0.05): ", n.r),
      paste("Numbers of reject (alpha=0.05): ", r),
      paste("Proportion of not reject: ", n.r/runs),
      sep = "\n")
  return(k.vector)
}

set.seed(2)
k <- updown.test(10000, 1000)


#2################################
#(a)
pidata <- readLines("StatisticalSimulation/pi.txt")
pi.digit <- as.numeric(strsplit(as.character(pidata), "")[[1]][-c(1:2)])
barplot(table(pi.digit))
chisq.test(table(pi.digit))
gap.test(pi.digit, .2, .8)


#3################################
q3 <- read.table("StatisticalSimulation/hw2_3.txt")
numbers <- as.vector(as.matrix(q3[, -c(1, 9)]))
hist(numbers, breaks = seq(0, 42, 4.2))

v <- floor(numbers * .25)
chisq.test(table(v))
gap.test(numbers, .2, .8)


#4################################
#(a)
#Box-Muller
boxmuller.method <- function(runs){
  p.value <- NULL
  p1 <- 0
  p2 <- 0
  p3 <- 0
  for(i in 1:runs){
    u1 <- runif(10000)
    u2 <- runif(10000)
    theta <- 2 * pi * u1
    k <- -log(u2)
    r <- sqrt(2 * k)
    x <- r * cos(theta)
    y <- r * sin(theta)
    p <- ks.test(c(x,y), "pnorm")$p.value
    if(p <= 0.01){
      p1 <- p1 + 1
    }else if(p <= 0.05){
      p2 <- p2 + 1
    }else if(p <= 0.1){
      p3 <- p3 + 1
    }
    p.value <- c(p.value, p)
  }
  cat("p-value under 0.01 :",p1,"\n",
      "p-value under 0.05 :",p1 + p2,"\n",
      "p-value under 0.1 :",p1 + p2 + p3,"\n")
  return(p.value)
}

#Polar Method
polar.method <- function(runs){
  p.value <- NULL
  p1 <- 0
  p2 <- 0
  p3 <- 0
  for(i in 1:runs){
    v1 <- runif(10000, min = -1, max = 1)
    v2 <- runif(10000, min = -1, max = 1)
    w <- v1^2 + v2^2
    w1 <- which(w < 1)
    w2 <- cbind(v1, v2, w)
    w2 <- w2[c(w1),]
    c <- sqrt(-2 * log(w2[,3]) / w2[,3])
    x <- c * w2[,1]
    y <- c * w2[,2]
    p <- ks.test(c(x,y), "pnorm")$p.value
    if(p <= 0.01){
      p1 <- p1 + 1
    }else if(p <= 0.05){
      p2 <- p2 + 1
    }else if(p <= 0.1){
      p3 <- p3 + 1
    }
    p.value <- c(p.value, p)
  }
  cat("p-value under 0.01 :",p1,"\n",
      "p-value under 0.05 :",p1 + p2,"\n",
      "p-value under 0.1 :",p1 + p2 + p3,"\n")
  return(p.value)
}

#Ratio of uniforms
ratio.of.uniforms <- function(runs){
  p.value <- NULL
  p1 <- 0
  p2 <- 0
  p3 <- 0
  for(i in 1:runs){
    u1 <- runif(10000)
    u2 <- runif(10000)
    v <- sqrt(2/exp(1)) * (2*u2-1)
    x <- v / u1
    z <- x^2 / 4
    z1 <- which(z <= (0.259/u1) + 0.35 & z <= -log(u1))
    z2 <- cbind(x, z)
    z2 <- z2[c(z1),]
    p <- ks.test(z2[,1], "pnorm")$p.value
    if(p <= 0.01){
      p1 <- p1 + 1
    }else if(p <= 0.05){
      p2 <- p2 + 1
    }else if(p <= 0.1){
      p3 <- p3 + 1
    }
    p.value <- c(p.value, p)
  }
  cat("p-value under 0.01 :",p1,"\n",
      "p-value under 0.05 :",p1 + p2,"\n",
      "p-value under 0.1 :",p1 + p2 + p3,"\n")
  return(p.value)
}

#rnorm()
r.norm <- function(runs){
  p.value <- NULL
  p1 <- 0
  p2 <- 0
  p3 <- 0
  for(i in 1:runs){
    u <- rnorm(10000)
    p <- ks.test(u, "pnorm")$p.value
    if(p <= 0.01){
      p1 <- p1 + 1
    }else if(p <= 0.05){
      p2 <- p2 + 1
    }else if(p <= 0.1){
      p3 <- p3 + 1
    }
    p.value <- c(p.value, p)
  }
  cat("p-value under 0.01 :",p1,"\n",
      "p-value under 0.05 :",p1 + p2,"\n",
      "p-value under 0.1 :",p1 + p2 + p3,"\n")
  return(p.value)
}

set.seed(2)
a <- boxmuller.method(1000)
b <- polar.method(1000)
c <- ratio.of.uniforms(1000)
d <- r.norm(1000)
par(mfrow = c(2, 2))
hist(a, main = 'p-value of ks.test(boxmuller.method)')
hist(b, main = 'p-value of ks.test(polar.method)')
hist(c, main = 'p-value of ks.test(ratio.of.uniforms)')
hist(d, main = 'p-value of ks.test(r.norm)')
ks.test(a, "punif")
ks.test(b, "punif")
ks.test(c, "punif")
ks.test(d, "punif")

#4(b)
u <- runif(2)
theta <- 2 * pi * u[1]
k <- -log(u[2])
r <- sqrt(2 * k)
x <- r * cos(theta)
y <- r * sin(theta)
xy.vector <- NULL
for(i in 1:100){
  x <- (131*(x)) %% (2^35)
  y <- (131*(y)) %% (2^35)
  xy.vector <- c(xy.vector,x,y)
}
sum((-3.3 < xy.vector & xy.vector < 3.6))/length(xy.vector)
 

#5################################
#(a)
x <- seq(-100, 100, 0.1)
c <- max(dcauchy(x) / dt(x, 0.5))

x.p <- NULL
t <- 0
set.seed(2)
repeat{
  u1 <- runif(1)
  u2 <- runif(1)
  x <- tan(pi*(u1 - 1/2))
  k <- dcauchy(x) / dt(x, 0.5) / c
  if(u2 <= k){
    x.p <- c(x.p, x)
  }
  t <- t + 1
  if (length(x.p) == 1000) break
}

(acc.rate <- 1000 / t)
x <- seq(-10, 10, 0.1)
y <- dcauchy(x)
plot(x, y, type="l", lty = 2, main = "Density and Simulated of Cauchy(0, 1)")
lines(density(x.p), lty = 1)
legend("topright", inset=.05,
       c("Simulated", "Theoretic"), lty = c(1, 2))

ks.test(x.p, "pcauchy")
permutation.test(x.p)

#(b)
leng1 <- NULL
ks1 <- NULL
ks2 <- NULL
set.seed(2)
for(i in 1:1000){
  x1<- NULL
  for(j in 1:1000){
    u1 <- runif(1)
    u2 <- runif(1)
    v <- 2*u2 - 1
    if (u1^2 + v^2 < 1){
      x <- v/u1
      x1 <- c(x1, x)
      }else{
        next
      }
    }

  x2 <- NULL
  for(j in 1:1000){
    u1 <- rnorm(1)
    u2 <- rnorm(1)
    x <- u1/u2
    x2 <- c(x2, x)
    }

  leng1 <- c(leng1, length(x1)/10000)
  ks1 <- c(ks1, ks.test(x1, "pcauchy")$p.value)
  ks2 <- c(ks2, ks.test(x2, "pcauchy")$p.value)
}

par(mfrow = c(1, 2))
hist(ks1, main = 'p-value of ks.test(r.o.unif)')
hist(ks2, main = 'p-value of ks.test(r.o.norm)')
ks.test(ks1, "punif")
ks.test(ks2, "punif")
mean(leng1)



#6################################
#Table method
table.method <- function(runs){
  x.vector <- NULL
  for(i in 1:runs){
    u <- runif(1)
    x1 <- c(rep(0,2),rep(1,4),rep(2,2),rep(3,0))
    x2 <- c(rep(0,9),rep(1,4),rep(2,2),rep(3,3))
    x3 <- c(rep(0,6),rep(1,4),rep(2,2),rep(3,7))
    if(u < 0.8){
      j <- floor(10*u) + 1        
      x <- x1[j]
    }else if(u < 0.98){
      j <- floor(100*u) - 80 + 1
      x <- x2[j]
    }else{
      j <- floor(1000*u) - 980 + 1
      x <- x3[j]
    }
    x.vector <- c(x.vector, x)
  }
  return(table(x.vector))
}

#Alias method
alias.run <- function(n) {
  temp <- NULL
  for (i in 1:n) {
    x <- floor(4*runif(1))
    x1 <- (runif(1) < (27-4)/27)
    x2 <- (runif(1) < (27-24)/27)
    x3 <- (runif(1) < (27-25)/27)
    xx <- c(0,1,2,2)*c(0,x3,x2,x1)
    y <- x - xx[c(x+1)]
    temp <- c(temp,y)
  }
  cat(x1,x2,x3,"\n")
  return(temp)
} 

set.seed(2)
p <- table.method(10000)
proc.time()
k <- alias.run(10000)
proc.time()
sum(abs(p-c(8/27, 12/27, 6/27, 1/27)*10000))
sum(abs(table(k)-c(8/27, 12/27, 6/27, 1/27)*10000))

