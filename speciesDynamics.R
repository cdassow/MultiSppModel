# 4.28.20 CD
# adapted version of the Biggs et al. 2009 model & some code do show how the model produces alternative stable states
library(deSolve)
rm(list=ls())
#parameter definitions
# qE - harvest, species specific
# s - survial, species specific
# cJA - effect of adults of a given species on juveniles of a given species (cover cannibalism or interspecific predation)
# cJJ - effect of juveniles of one species on juveniles of the other (can be predation or competition)
# h - rate at which juveniles leave foraging arena for refuge, species specific
# v - rate at which juveniles enter foraging arena from refuge, species specific
# f - fecundity, species specific


simBiggs<-function(t,y,params){
  A1<-y[1]
  A2<-y[2]
  J1<-y[3]
  J2<-y[4]
  with(as.list(params),{
    dA1dt=-qE1*A1+(s1-1)*A1+s1*J1
    dA2dt=-qE2*A2+(s2-1)*A2+s2*J2
    dJ1dt=-cJ1A1*J1*A1-cJ1J2*J2*J1-(cJ1A2*v1*A2*J1)/(h1+v1+cJ1A2*A2)+f1*A1
    dJ2dt=-cJ2A2*J2*A2-cJ2J1*J1*J2-(cJ2A1*v2*A1*J2)/(h2+v2+cJ2A1*A1)+f2*A2
    return(list(c(dA1dt,dA2dt,dJ1dt,dJ2dt)))
  })
}
times=1:500
p=c(qE1=1.8,s1=0.5,cJ1A1=0.008,cJ1A2=0.003,cJ1J2=0.001,v1=1,h1=6,f1=2,
    qE2=1.8,s2=0.5,cJ2A2=0.008,cJ2A1=0.003,cJ2J1=0.001,v2=1,h2=6,f2=2)
y0=c(50,20,10,10)
sim=ode(y=y0,times=times,func=simBiggs,parms=p)
sim[nrow(sim),]
plot(sim[,1],sim[,2],type='l',ylim=c(0,max(sim[,2:5])))
lines(sim[,1],sim[,3],col="blue")
lines(sim[,1],sim[,4],lty=2)
lines(sim[,1],sim[,5],lty=2,col="blue")

### RUNNING OVER A RANGE OF HARVESTS - TEST FOR ALTERNATE STABLE STATES

store=data.frame(qEs=seq(.05,5,length.out=30),A1=0,A2=0,J1=0,J2=0)
y0=c(10,200,40,40)
for(i in 1:nrow(store)){
  p=c(c(qE1=store$qEs[i],s1=0.5,cJ1A1=0.007,cJ1A2=0.003,cJ1J2=0.001,v1=1,h1=6,f1=1),
      c(qE2=1.8,s2=0.5,cJ2A2=0.007,cJ2A1=0.003,cJ2J1=0.001,v2=1,h2=6,f2=1))
  sim=ode(y=y0,times=times,func=simBiggs,parms=p)
  store$A1[i]=sim[nrow(sim),2]
  store$A2[i]=sim[nrow(sim),3]
  store$J1[i]=sim[nrow(sim),4]
  store$J2[i]=sim[nrow(sim),5]
}
store2=data.frame(qEs=seq(.05,5,length.out=30),A1=0,A2=0,J1=0,J2=0)
y0=c(200,10,40,40)
for(i in 1:nrow(store)){
  p=c(c(qE1=store2$qEs[i],s1=0.5,cJ1A1=0.007,cJ1A2=0.003,cJ1J2=0.001,v1=1,h1=6,f1=1),
      c(qE2=1.8,s2=0.5,cJ2A2=0.007,cJ2A1=0.003,cJ2J1=0.001,v2=1,h2=6,f2=1))
  sim=ode(y=y0,times=times,func=simBiggs,parms=p)
  store2$A1[i]=sim[nrow(sim),2]
  store2$A2[i]=sim[nrow(sim),3]
  store2$J1[i]=sim[nrow(sim),4]
  store2$J2[i]=sim[nrow(sim),5]
}
plot(store$qEs,store$A1,lwd=3,type='l',ylim=c(0,max(store[,2:3])),ylab = "Abundance",xlab = "Harvest (q*E)")
lines(store$qEs,store$A2,lwd=3,col='grey')
lines(store2$qEs,store2$A1,lwd=3,lty=3)
lines(store2$qEs,store2$A2,lwd=3,col='grey',lty=3)
legend("topright",legend = c("sp 1", "sp 2","run1","run2"), col = c("black","grey","black","black"),lwd=2,lty = c(1,1,1,3),bty="n")

