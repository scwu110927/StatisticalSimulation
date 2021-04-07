##-HW1-################
#1(a)



#2(a)
pi <- read.table("StatisticalSimulation/pi.txt", colClasses="character")
pi.digit <- as.numeric(strsplit(as.character(pi), "")[[1]][-c(1:2)])

hist(pi.digit)
chisq.test(table(pi.digit))





#3
q3 <- read.table("StatisticalSimulation/hw2_3.txt")
numbers <- as.vector(as.matrix(q3[, -c(1, 9)]))
hist(numbers)

v <- floor(numbers * .25)
chisq.test(table(v))
hist(data)

gap.test <- function(data, a, b){
  data <- data/42
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

gap.test(numbers, .2, .8)

