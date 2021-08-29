# Nomogram

library(rms)
library(Hmisc)
library(lattice)
library(survival)
library(Formula)
library(ggplot2)
library(foreign)
options(stringsAsFactors=FALSE)

bc <- read.csv("Nomogram1.csv")
bc[1:5,1:5]
length(bc[1, ])
length(bc[ ,1])
dim(bc)

bc[,"futime"]=bc[,"futime"]*365
dim(bc)
bc[1:5,]
bc
max(bc$futime)

bc <- na.omit(bc)
dim(bc)
bc[1:5,]
length(bc[1, ])
length(bc[ ,1])
dim(bc)

dd <- datadist(bc)
options(datadist="dd")


f <- cph(Surv(futime, fustat) ~ riskScore+M_Stage+IFNg_signature_6+CTL, x=T, y=T, surv=T, data=bc, time.inc=365*3)
surv <- Survival(f)
surv(365)
surv(730)

med <- Quantile(f)
surv <- Survival(f)

nom <- nomogram(f, fun=function(x) med(lp=x),
	 funlabel="Median Survival Time")
plot(nom)

nom <- nomogram(f, fun=list(function(x) surv(365, x), function(x) surv(730, x)), lp=F, funlabel=c("Probability of 1-year survival", "Probability of 2-year survival"), maxscale=100, fun.at=c(0.99, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3,0.2,0.1,0.05,0.03)) 
nom <- nomogram(f, fun=list(function(x) surv(365, x), function(x) surv(730, x)), lp=F, funlabel=c("Probability of 1-year survival", "Probability of 2-year survival"), maxscale=100, fun.at=c(0.987, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3,0.2,0.1,0.05,0.03)) 
plot(nom)
print(nom)
validate(f, method="boot", B=1000, dxy=T)
rcorrcens(Surv(futime, fustat) ~ predict(f), data = bc)



# Calibration curve
f1 <- cph(Surv(futime, fustat) ~ riskScore+M_Stage+IFNg_signature_6+CTL, x=T, y=T, surv=T, data=bc, time.inc=365)
cal1 <- calibrate(f1, cmethod="KM", method="boot", u=365, m=14,B=1000)
plot(cal1)

f3 <- cph(Surv(futime, fustat) ~ riskScore+M_Stage+IFNg_signature_6+CTL, x=T, y=T, surv=T, data=bc, time.inc=730)
cal3 <- calibrate(f3, cmethod="KM", method="boot", u=730, m=14, B=1000)
plot(cal3)