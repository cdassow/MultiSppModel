#place to build the species portion of the multispecies fishery model
#1.17.2020

rm(list = ls())
setwd("C:/Users/jones/BoxSync/NDstuff/Dissertation/4/MultiSppModel/") #working directory for Colin's computer

#start with a stage-structured model with 2 stages for a single species
years=100
A.A=numeric(years)
A.J=numeric(years)
B.A=numeric(years)
B.J=numeric(years)

start=c(500,3000,400,3000)

Amort=.01 #natural mortality
Bmort=.01 #natural mortality
Arep=.5 #reproduction
Brep=.5 #reproduction
Amat=.2 #maturation
Bmat=.2 #maturation
Apred=3 #number of juveniles of B consumbed by A adults each year
Bpred=3 #number of juveniles of A consumbed by B adults each year
K.aa=1000 #carrying capacity of A adults
K.aj=5000 #carrying capacity of A juv
K.ba=1000 #carrying capacity of B adults
K.bj=5000 #carrying capacity of B juv
q.a=.1 #catchability of A
e.a=1 #effort for A
q.b=.1 #catchability of B
e.b=1 #effort for B
Aharv=0.9 #proportion of the catch harvested for A
Bharv=0.01 #proportion of the catch harvested for B

A.A[1]=start[1]
A.J[1]=start[2]
B.A[1]=start[3]
B.J[1]=start[4]
for(i in 1:years){
  A.A[i+1]=ifelse(A.A[i]+(A.J[i]*Amat-A.A[i]*Amort)*(1-(A.A[i]/K.aa))-(q.a*e.a*A.A[i]*Aharv)<0,
                  0,
                  A.A[i]+(A.J[i]*Amat-A.A[i]*Amort)*(1-(A.A[i]/K.aa))-(q.a*e.a*A.A[i]*Aharv))
  A.J[i+1]=ifelse(A.J[i]+(A.A[i]*Arep-B.A[i]*Bpred)*(1-(A.J[i]/K.aj))<0,
                  0,
                  A.J[i]+(A.A[i]*Arep-B.A[i]*Bpred)*(1-(A.J[i]/K.aj)))
  B.A[i+1]=ifelse(B.A[i]+(B.J[i]*Bmat-B.A[i]*Bmort)*(1-(B.A[i]/K.ba))-(q.b*e.b*B.A[i]*Bharv)<0,
                  0,
                  B.A[i]+(B.J[i]*Bmat-B.A[i]*Bmort)*(1-(B.A[i]/K.ba))-(q.b*e.b*B.A[i]*Bharv))
  B.J[i+1]=ifelse(B.J[i]+(B.A[i]*Brep-A.A[i]*Apred)*(1-(B.J[i]/K.bj))<0,
                  0,
                  B.J[i]+(B.A[i]*Brep-A.A[i]*Apred)*(1-(B.J[i]/K.bj)))
  
    
}
abund=data.frame(A.A, A.J, B.A, B.J) #storing the abundances of adults an juveniles at each time step


#visualize output

plot(1:nrow(abund), abund$A.A, type = "l", ylab = "Abundance", xlab = "Time", ylim = range(abund))
lines(1:nrow(abund), abund$A.J, lty=2)
lines(1:nrow(abund), abund$B.A, lty=1, col="red")
lines(1:nrow(abund), abund$B.J, lty=2, col="red")
legend("right", legend = c("A.adult", "A.juv", "B.adult", "B.juv"), lty = c(1,2,1,2), col = c("black", "black", "red", "red"))
