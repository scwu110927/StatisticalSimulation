##-HW2-################
#1################################
#(a)
gap.test <- function(data, a, b){
  data <- data/max(data)
  n <- length(data)
  x <- c(1:n) * (a < data & data < b)
  x1 <- x[x>0]
  y <- x1[-1]-x1[-length(x1)]-1
  t.y <- table(y)
  e <- (b-a) * (1-(b-a))^(c(1:dim(t.y))-1) * sum(t.y)
  e.p <- e / sum(e)
  c.t <- chisq.test(t.y, p = e.p)
  return(list("gaps" = t.y, "expected count" = round(e), "chisq.test" = c.t))
}

rn10000 <- read.csv("StatisticalSimulation/rn10000.csv", header = F)
rn10000 <- as.matrix(rn10000)
gap.test(rn10000, 0.2, 0.7)


permutation.test <- function(data){
  x1 <- matrix(data[-1], ncol = length(data) %/% 3, byrow = F)
  y1 <- apply(x1, 2, rank)
  y2 <- y1[1,]*100 + y1[2,]*10 + y1[3,]
  return(table(y2))
}

y <- permutation.test(rn10000)
chisq.test(y)

#(b)
updown.test <- function(num,runs){
  k.vector <- NULL
  p <- 0
  q <- 0
  g <- 0
  for(i in 1:runs){
    x <- runif(num)
    x1 <- (x[-1]>x[-num])
    x2 <- sum((x1[-1] != x1[-(num-1)]))
    z <- (x2-(2*num-1)/3)/sqrt((16*num-29)/90)
    k <- pnorm(z)
    k.vector <- c(k.vector,k)
    if(k <= 0.01){
      p <- p+1
    }else if(k <= 0.05){
      q <- q+1
    }else if(k <= 0.1){
      g <- g+1
    }
  }
  cat("p-value under 0.01", p,
      "p-value between 0.01 and 0.05", q,
      "p-value between 0.05 and 0.1", g,
      "p-value better than 0.1", runs-p-q-g,
      sep = "\n")
  return(k.vector)
}

k <- updown.test(10000, 1000)
hist(k, xlab = 'p-value', main = 'z.table.pvalue')
ks.test(k, y = "punif") #reject H0 not uniform
chisq.test(table(ceiling(k*10))) #reject H0 not uniform


#2################################
#(a)
pidata <- read.table("StatisticalSimulation/pi.txt", colClasses="character")
pi.digit <- as.numeric(strsplit(as.character(pidata), "")[[1]][-c(1:2)])

hist(pi.digit)
chisq.test(table(pi.digit))
gap.test(pi.digit, .2, .8)


#3################################
q3 <- read.table("StatisticalSimulation/hw2_3.txt")
numbers <- as.vector(as.matrix(q3[, -c(1, 9)]))
hist(numbers)

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
    theta <- 2*pi*u1
    k <- -log(u2)
    r <- sqrt(2*k)
    x <- r*cos(theta)
    y <- r*sin(theta)
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
      "p-value under 0.05 :",p1+p2,"\n",
      "p-value under 0.1 :",p1+p2+p3,"\n")
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
    w <- v1^2+v2^2
    w1 <- which(w < 1)
    w2 <- cbind(v1, v2, w)
    w2 <- w2[c(w1),]
    c <- sqrt(-2*log(w2[,3])/w2[,3])
    x <- c*w2[,1]
    y <- c*w2[,2]
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
      "p-value under 0.05 :",p1+p2,"\n",
      "p-value under 0.1 :",p1+p2+p3,"\n")
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
    v <- sqrt(2/exp(1))*(2*u2-1)
    x <- v/u1
    z <- x^2/4
    z1 <- which(z <= (0.259/u1)+0.35 & z <= -log(u1))
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
      "p-value under 0.05 :",p1+p2,"\n",
      "p-value under 0.1 :",p1+p2+p3,"\n")
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
      "p-value under 0.05 :",p1+p2,"\n",
      "p-value under 0.1 :",p1+p2+p3,"\n")
  return(p.value)
}

a <- boxmuller.method(1000)
b <- polar.method(1000)
c <- ratio.of.uniforms(1000)
d <- r.norm(1000)

#best is rnorm 
  
#4(b)
u <- runif(2)
theta <- 2*pi*u[1]
k <- -log(u[2])
r <- sqrt(2*k)
x <- r*cos(theta)
y <- r*sin(theta)
xy.vector <- c(x,y)
for(i in 1:100){
  x <- (131*(x)) %% (2^35)
  y <- (131*(y)) %% (2^35)
  xy.vector <- c(xy.vector,x,y)
}
xy.vector
sum(-3.3 < xy.vector & xy.vector < 3.6)/length(xy.vector)
 


#5
#(a)
x <- seq(-10, 10, 0.1)
c <- max(dcauchy(x)/dt(x, 0.5))
x.p <- NULL
t <- 0

for(i in 1:1000){
  t1 <- rt(1, 0.5)
  t2 <- rt(1, 0.5)
  k <- (dcauchy(t1)/dt(t1, 0.5))/c
  if(t2 <= k){
    x.p <- c(x.p, t1)
  }else{
    next
  }
}

length(x.p)/1000
hist(x.p)
ks.test(x.p, "pcauchy")
  
#(b)
leng1 <- NULL
ks1 <- NULL
ks2 <- NULL
for(i in 1:10){
  x1<- NULL
  for(j in 1:10000){
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
  for(j in 1:10000){
    u1 <- rnorm(1)
    u2 <- rnorm(1)
    x <- u1/u2
    x2 <- c(x2, x)
    }

  leng1 <- c(leng1, length(x1)/10000)
  ks1 <- c(ks1, ks.test(x1, "pcauchy")$p.value)
  ks2 <- c(ks2, ks.test(x2, "pcauchy")$p.value)
}

mean(leng1)
mean(ks1)
mean(ks2)


#6
#table.method
x0 <- dbinom(0,3,1/3) #0.2962
x1 <- dbinom(1,3,1/3) #0.4444
x2 <- dbinom(2,3,1/3) #0.2222
x3 <- dbinom(3,3,1/3) #0.0370

table.method <- function(runs){
  x.vector <- NULL
  for(i in 1:runs){
    u <- runif(1)
    x1 <- c(rep(0,2),rep(1,4),rep(2,2),rep(3,0))
    x2 <- c(rep(0,9),rep(1,4),rep(2,2),rep(3,3))
    x3 <- c(rep(0,6),rep(1,4),rep(2,2),rep(3,7))
    if(u < 0.8){
      j <- floor(10*u)+1        
      x <- x1[j]
    }else if(u < 0.98){
      j <- floor(100*u)-80+1
      x <- x2[j]
    }else{
      j <- floor(1000*u)-980+1
      x <- x3[j]
    }
    x.vector <- c(x.vector, x)
  }
  return(table(x.vector))
}

p <- table.method(10000)
c(p[1]/10000, p[2]/10000, p[3]/10000, p[4]/10000)
c(x0, x1, x2,x3)

#the Alias method
alias.run=function(n) {
 temp=NULL
for (i in 1:n) {
x=floor(4*runif(1))
x2=(runif(1) < 2/27)*1
x3=(runif(1) < 3/27)*1
 x4=(runif(1) < 23/27)*1
 xx=c(0,1,2,3)*c(0,x2,x3,x4)
 y=x-xx[c(x+1)]
temp=c(temp,y)
}
return(temp)
 } 
alias.run(100)


