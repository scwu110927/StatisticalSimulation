#HW3
#1################################
library(ISR3)
x1 <- c(rep(1300,6),rep(1200,6), rep(1100,4))
x2 <- c(7.5, 9.0, 11.0, 13.5, 17.0,
        23.0, 5.3, 7.5, 11.0, 13.5,
        17.0, 23.0, 5.3, 7.5, 11.0,
        17.0)
x3 <- c(0.012, 0.012, 0.0115, 0.013, 0.0135,
        0.012, 0.040, 0.038, 0.032, 0.026,
        0.034, 0.041, 0.084, 0.098, 0.092,
        0.086)
y <- c(49.0, 50.2, 50.5, 48.5, 47.5,
       44.5, 28.0, 31.5, 34.5, 35.0,
       38.0, 38.5, 15.0, 17.0, 20.5,
       29.5)

x4 <- x1 * x1
x5 <- x2 * x2
x6 <- x3 * x3
x7 <- x1 * x2
x8 <- x1 * x3
x9 <- x2 * x3

X <- matrix(c(rep(1,16), x1, x2, x3, x4, x5, x6, x7, x8, x9), ncol = 10)

#Cholesky
A <- t(X) %*% X
L <- t(chol(A))
B <- solve(t(L), solve(L, t(X) %*% y)) 
round(B, 2)

#SWEEP
M <- matrix(rbind(cbind(A, t(X) %*% y), cbind(t(y) %*% X, t(y) %*% y)), nrow = 11)
M2 <- SWP(M, 1:10)
M2

XX <- -M2[1:10,1:10]
round(sqrt(diag(XX * (M2[11,11]/6))) , 2)

g <- lm(y~x1+x2+x3+x4+x5+x6+x7+x8+x9)
summary(g)

#2################################
plot(lynx, type= 'l', main = 'Time series plot of Lynx')
ar.ols(lynx, order = 1)
ar.ols(lynx, order = 2)


#3################################
library(forecast)
mortality <- read.csv("StatisticalSimulation/malemortality.csv", h = T)
x <- scale(log(mortality[mortality$year <= 99, -1]), scale = FALSE)
svd.x <- svd(x, 1, 1)
b <- svd.x$v
k <- svd.x$u * svd.x$d[1]
fit.k <- ar.ols(k, order = 1)
pred.k <- NULL
p.k <- k[19]
for(i in 1:5){
  p.k <- fit.k$x.intercept + fit.k$ar * p.k
  pred.k <- c(pred.k, p.k)
}
xnew.appr <- pred.k %*% t(b)
round(xnew.appr, 2)
xnew <- scale(log(mortality[mortality$year > 99, -1]), scale = FALSE)
sum((xnew - xnew.appr)^2)


pca.x <- princomp(x)

#4################################
#(a)
library(gtools)
DDT <- c(65, 98, 117, 122, 130)
egg <- c(0.52, 0.53, 0.49, 0.49, 0.37)
p <- permutations(5,5)
DDT2 <- matrix(DDT[t(p)], byrow = T ,ncol = 5)
length(which(DDT2 %*% sort(egg) <= sum(DDT * egg)))/nrow(DDT2) 
cor.test(DDT, egg, method = "pearson") 
cor.test(DDT, egg, method ="spearman") 


x <- c(585, 1002, 472, 493, 408, 690, 291)
y <- c(0.1, 0.2, 0.5, 1, 1.5, 2, 3)
p2 <- permutations(7,7)
x2 <- matrix(x[t(p2)], byrow = T, ncol = 7)
length(which(x2 %*% sort(y) <= sum(x * y)))/nrow(x2)
cor.test(x, y, method="pearson") 
cor.test(x, y, method="spearman") 

#(b)
choln2u <- function(rho = 0.2){
  A <- matrix(c(1, rho, rho, 1), ncol=2)
  B <- t(chol(A))
  t1 <- NULL
  t2 <- NULL
  t3 <- NULL
  for (i in 1:100) { 
    x1 <- rnorm(10)
    y1 <- rnorm(10)
    xy1 <- rbind(x1, y1)
    xy <- B %*% xy1
    xy.unif <- floor(10*pnorm(xy)) 
    z1 <- cor.test(xy.unif[1, ], xy.unif[2, ], method = "pearson")$p.value
    z2 <- cor.test(xy.unif[1, ], xy.unif[2, ], method = "spearman")$p.value
    t1 <- c(t1, z1)
    t2 <- c(t2, z2)
    t <- NULL
    z0 <- sum(x1*y1)
    for (j in 1:1000) {
      x0 <- sample(x1, 10, F)
      y0 <- sample(y1, 10, F)
      t <- c(t, sum(x0*y0))
    }
    z3 <- sum(t >= z0)/1000
    t3 <- c(t3, z3)
  }
  print(paste('pearson p<0.05個數:', length(t1[t1<0.05])))
  print(paste('spearman p<0.05個數:', length(t2[t2<0.05])))
  print(paste('permutation test p<0.05個數:', length(t2[t2<0.05])))
}

choln2u(0.2)
choln2u(0.8)

#5################################
MWWTest = function(seq1, seq2){
  seq = c(seq1, seq2)
  n1 = length(seq1)
  n2 = length(seq2)
  
  r1 = rank(seq)[1:n1]
  r2 = rank(seq)[-c(1:n1)]
  
  W1 = sum(r1)
  W2 = sum(r2)
  
  U1 = n1 * n2 + (n1 * (n1 + 1) / 2) - W1
  U2 = n1 * n2 + (n2 * (n2 + 1) / 2) - W2
  
  U = min(U1, U2)
  Z = (U - n1 * n2 / 2) / sqrt(n1 * n2 * (n1 + n2 +1) / 12)
  df = data.frame(U1, U2)
  return (U1)
}

data = NULL
for (i in 2:10){
  for (j in 2:10){
    x = NULL
    for (k in 1:10000){
      y = MWWTest(runif(i), runif(j))
      x = c(x, y)
    }
    data = cbind(data, sort(as.numeric(x)))
  }
}

df = data.frame(data, row.names = c(1:10000))
critical = df[500, ]
cri_mat = matrix(critical, ncol = 9, byrow = TRUE)
cri_df = data.frame(cri_mat, row.names = c(2:10))
colnames(cri_df) = c(2:10)
print(cri_df)


#6################################
library(bootstrap)

law.bs <- function(x,y){
  t1 <- NULL
  t2 <- NULL
  
  for (i in y) { 
    x1 <- sample(1:15, x, T)
    x2 <- law[x1,]
    x3 <- array(apply(x2, 2, var))
    t1 <- c(t1, x3[1]) 
    t2 <- c(t2, x3[2])
    
  }
  my_list <- list("Original Variance (LSAT & GPA)" = c(var(law$LSAT), var(law$GPA)),
                  "Bootstrap Variance (LSAT & GPA) of 10 r.v" =  ,
                  "Bootstrap Variance (LSAT & GPA) of 15 r.v" =  ,
                  "Bootstrap Variance (LSAT & GPA) of 20 r.v" =  ,
                  "Bootstrap Variance (LSAT & GPA) of 25 r.v" =  )
  return(my_list)
  return(c(mean(t1), mean(t2)))
}


law.bs <- function(x,y){
  t1 <- NULL
  t2 <- NULL
  t3 <- NULL
  t4 <- NULL
  for(j in 1:200){
    k <- 50*j
    for (i in k) { 
      x1 <- sample(1:15, x, T)
      x2 <- law[x1,]
      x3 <- array(apply(x2, 2, var))
      t1 <- c(t1, x3[1]) 
      t2 <- c(t2, x3[2])
    }
    t3 <- c(t3, mean(t1))
    t4 <- c(t4, mean(t2))
  }
  par(mfrow=c(1,2))
  p1 <- rbind(plot(t3, type = 'l'), abline(h = var(law$LSAT), lty = 2))
  p2 <- rbind(plot(t4, type = 'l'), abline(h = var(law$GPA), lty = 2))
  return(list(p1, p2))
}
law.bs(10)
law.bs(15)
law.bs(20)
law.bs(25)



