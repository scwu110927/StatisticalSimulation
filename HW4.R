#HW3
#1################################




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

