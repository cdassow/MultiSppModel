#place to build the species portion of the multispecies fishery model
#1.17.2020

rm(list = ls())
setwd("C:/Users/jones/BoxSync/NDstuff/Dissertation/4/MultiSppModel/") #working directory for Colin's computer

#start with a stage-structured model with 2 stages for a single species

#places to store model output
years=1000
a.a=numeric(years) #adults spec. A
j.a=numeric(years) #juveniles spec. A
a.b=numeric(years) #adults spec. B
j.b=numeric(years) #juveniles spec. B
p.a=numeric(years) #other prey for spec. A
p.b=numeric(years) #other prey for spec. B
j.r=numeric(years) #shared resource for juveniles

#natural mortality of each species
mortA=0.15
mortB=0.15
#maturation of each species
matA=0.2
matB=0.2
#beverton-holt recruitment parameters for each species
xa=150 #stock size at which recruitment = half the max
xb=150 #stock size at which recruitment = half the max
ya=200 #maximum recruits
yb=200 #maximum recruits
#predation parameters for each species
ha=100 #handle time 
hb=100 #handle time
ea=001 #encounter rate for a juveniles with adults of b
eb=0.1 #encounter rate for b juveniles with adults of a
j.bRate=0.2 #Percentage of resources consumed by adults of A that is juveniles of B
j.aRate=0.2 #Percentage of resources consumed by adults of B that is juveniles of A
#alternative prey & juvenile resource parameters
K.p.a=400 #carrying capacity for spec. A other prey
K.p.b=400 #carrying capacity for spec. B other prey
r.p.a=0.8 #growth coef for spec. A other prey
r.p.b=0.8 #growth coef for spec. B other prey
ha.p=100 # A adult handling time for other prey
hb.p=100 # B adult handling time for other prey
ea.p=0.01 # A adult encounter rate for other prey
eb.p=0.01 # B adult encounter rate for other prey
#trying a linear effect of each species juveniles on the other
efc.a=0.3 #seq(0.05, 0.1, length.out = years)
efc.b=0.3 #seq(0.1, 0.05, length.out = years)
#harvest parameters
qa=0.1 #catchability of A
qb=0.1 #catchability of B
Ea=1 #effort for A
Eb=1 #effort for B
harvA=c(rep(0, years/2), rep(0, years/2)) #harvest rate for A
harvB=c(rep(0, years/2), rep(0, years/2)) #harvest rate for B

#initializing
start=c(500,100,100,100, 400, 400, 400)
a.a[1]=start[1]
j.a[1]=start[2]
a.b[1]=start[3]
j.b[1]=start[4]
p.a[1]=start[5]
p.b[1]=start[6]
j.r[1]=start[7]

for(i in 1:(years-1)){
  #predation calculations
  Apred=(ea.p*p.a[i]/(1+ea.p*ha.p*p.a[i]))*a.a[i]
  Bpred=(eb.p*p.b[i]/(1+eb.p*hb.p*p.b[i]))*a.b[i]
  #catch & harvest calculations
  Acatch=(qa*Ea*a.a[i]*harvA[i])
  bcatch=(qb*Eb*a.b[i]*harvB[i])
  #recruitment calculations
  Arecruit=((ya*a.a[i])/(xa+a.a[i]))
  Brecruit=((yb*a.b[i])/(xb+a.b[i]))
  #juveniles
  j.a[i+1]=j.a[i]+Arecruit-(mortA*j.a[i])-matA*j.a[i]-(Bpred*j.aRate)-efc.a*j.b[i]
  j.b[i+1]=j.b[i]+Brecruit-(mortB*j.b[i])-matB*j.b[i]-(Apred*j.bRate)-efc.b*j.a[i]
  #alternative prey & shared resources
  p.a[i+1]=p.a[i]+r.p.a*(1-(p.a[i]/K.p.a))*p.a[i]-Apred
  p.b[i+1]=p.b[i]+r.p.b*(1-(p.b[i]/K.p.b))*p.b[i]-Bpred
  #adults
  a.a[i+1]=a.a[i]+matA*j.a[i]-mortA*a.a[i]-Acatch+Apred+(Apred*j.bRate)
  a.b[i+1]=a.b[i]+matB*j.b[i]-mortB*a.b[i]-bcatch+Bpred+(Bpred*j.aRate)
  
}
#collecting output
abund=data.frame(j.a, j.b, a.a, a.b, p.a, p.b, j.r)
#looking at output
plot(1:years, abund$j.a, type = "l", lty = 3, ylim = c(0,max(abund, na.rm = T)), ylab = "abundance", xlab = "time")
lines(1:years, abund$j.b, lty = 3, col = "red")
lines(1:years, abund$a.a, lty = 1, col = "black")
lines(1:years, abund$a.b, lty = 1, col = "red")
lines(1:years, abund$p.a, lty = 5, col = "grey")
lines(1:years, abund$p.b, lty = 5, col = "darkred")
#lines(1:years, abund$j.r, lty = 1, col = "green")
#abline(v=years/2, col = "grey")
legend("bottomright", legend = c("j.a", "j.b", "a.a", "a.b"), lty = c(3,3,1,1), col = c("black", "red", "black", "red"), bty = "n")
#legend("bottomright", legend = c("j.a", "j.b", "a.a", "a.b", "p.a", "p.b"), lty = c(3,3,1,1,5,5), col = c("black", "red", "black", "red", "grey", "darkred"), bty = "n")


#storing juvenile resource predation calculations here while I try something else

# j.K=500 #carrying capacity for juvenile shared resource
# j.R=0.8 #growth coef for juvenile shared resource
# jPred.a=0.5 #encounter rate of A juveniles on shared resource
# jPred.b=0.5 #encounter rate of B juveniles on shared resource
# j.h=30 #handling time for juveniles shared resource
# juv A +((jPred.a*j.r[i])/(1+jPred.a*j.h*j.r[i])*j.a[i])
# juv B +((jPred.b*j.r[i])/(1+jPred.b*j.h*j.r[i])*j.b[i])
# resource pop dynamics j.r[i+1]=j.r[i]+j.R*(1-(j.r[i]/j.K))*j.r[i]-((jPred.a*j.r[i])/(1+jPred.a*j.h*j.r[i])*j.a[i])-((jPred.b*j.r[i])/(1+jPred.b*j.h*j.r[i])*j.b[i])


#playing around with some recruitment relationships

#beverton-holt recruitment parameters for each species
xa=400 #stock size at which recruitment = half the max
xb=500 #stock size at which recruitment = half the max
ya=100 #maximum recruits
yb=100 #maximum recruits
sizes=seq(0, 500, length.out = 10)

plot(sizes, ((ya*sizes)/(xa+sizes)), type = "l", ylab = "numRecruits", xlab = "adult abund")
lines(sizes, ((yb*sizes)/(xb+sizes)), col="red")
