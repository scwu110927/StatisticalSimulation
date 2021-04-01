##-HW1-################
#1(a)



#2(a)
pi <- read.table("StatisticalSimulation/pi.txt", colClasses="character")
pi.digit <- as.numeric(strsplit(as.character(pi), "")[[1]][-c(1:2)])

hist(pi.digit)
chisq.test(table(pi.digit))


