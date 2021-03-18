#1(a)
rep(seq(0,4),each=5)

#(b)
seq(1:5) + rep(0:4,each=5)

#(C) 123234345456
x=c("red","yellow","blue","green","magenta","cyan")
i=seq(1:3)+rep(0:3,each=3)
x[i]


#2.(a)
library(spatstat)

mindist <- function(x, y){
 min(nndist(x, y))
}

x <- runif(20)
y <- runif(20)
mindist(x, y)

#(b)
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


#3.
gcd = function(a,b)
{ if (b==0) a else gcd(b,a%%b) }

lcm = function(c,d)
{ return(c*d/gcd(c,d)) }



#4.
#(b)
x1 <- 0.6
x2 <- 0
for (i in 1:10){
#  x1 <- (69069*x1)-2^32*floor(69069*x1/2^32)
  x1 <- 69069*x1 %% 2^32
  x2 <- c(x2, x1)
}
x2 <- x2[-1]/2^32
hist(x2)
x2[1:10]
x1 <- (69069*x1)-2^32*floor(69069*x1/2^32)
69069*2862316056.6 %% 2^32
0.77777%%1

ks.test(x2, 1)

0.7%%1
