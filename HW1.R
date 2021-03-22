#1(a)
rep(seq(0,4),each=5)

#1(b)
seq(1:5) + rep(0:4,each=5)

#1(C) 123234345456
x=c("red","yellow","blue","green","magenta","cyan")
i=seq(1:3)+rep(0:3,each=3)
x[i]


#2(a)
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
index <- apply(xy, 2, grid)

plot(x, y, pch=index)
abline(h = 0.5, v = 0.5, lty=2)

#(c)rectangles
symbols(x, y, circles=index, add=T)


#3
gcd = function(a,b)
{ if (b==0) a else gcd(b,a%%b) }

lcm = function(c,d)
{ return(c*d/gcd(c,d)) }



#4(a)
midsqur <- function(seed,times){
  numvector<- NULL
  for(i in 1:times){
    num <- seed*seed
    seed <- (num %/% 1000) %% 1000000
    #num (1)先除以1000取整數部分，等同於去掉後三位數
    #(2)將剩下的數值除以1000000取餘數部分，
    #若剩下的數大於六位數，則取六位數
    #若少於六位數，則直接取該數
    numvector<-c(numvector,seed)
  }
  return(numvector)
}

(x <- ceiling(runif(1,0,999999))) #seed
(x1 <- midsqur(x,10000))
x1 <- (x1/10^6)
hist(x1)
ks.test(x1, y="punif")

for(i in 1:length(x1)){
  y <- x1[i]-x1
  ys <- sum(y == 0)
  yc <- x1[which(y == 0)]
}
ys
yc


#4(b)
x1 <- 0.6
x2 <- 0
for (i in 1:10000){
  x1 <- (69069*x1) %% 2^32
  x2 <- c(x2, x1)
}
(x2 <- x2[-1]/2^32)
hist(x2)

v1<-floor(x2*10) #分組以0.1為級距
table(v1)
ks.test(x2, y = "punif")
chisq.test(table(v1))


#4(c)
xi<-rnorm(1,0,1)
yi<-rnorm(1,0,1)
zi<-rnorm(1,0,1)
vector <- NULL

for (j in 1:10000) {
  xi <- (171*xi) %% 30269
  yi <- (172*yi) %% 30307
  zi <- (170*zi) %% 30323
  ui <- ((xi/30269)+(yi/30307)+(zi/30323)) %% 1
  vector <- c(vector,ui)
}

v1<-floor(vector*10) #分組以0.1為級距
table(v1)
ks.test(vector,y="punif")
chisq.test(table(v1))

#5(b)
#fi
casio <- function(seed,times){
  fi <- (1+sqrt(5))/2
  uvector <- NULL
  for(i in 1:times){
    seed <- ((fi+seed)^5) %% 1
    uvector<-c(uvector,seed)
  }
  return(uvector)
}

(l <- runif(1,0,1))
k <- casio(l,10000)
hist(k)

v1<-floor(k*10) #分組以0.1為級距
#hist(v1)
table(v1)


#sqrt(2)
casio2 <- function(seed,times){
  uvector <- NULL
  for(i in 1:times){
    seed <- ((sqrt(2)+seed)^5) %% 1
    uvector<-c(uvector,seed)
  }
  return(uvector)
}

k2<-casio2(l,10000)
hist(k2)

v2<-floor(k2*10) #分組以0.1為級距
#hist(v1)
table(v2)

chisq.test(table(v1))
chisq.test(table(v2))
