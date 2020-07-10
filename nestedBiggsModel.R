#recreating Biggs et al. 2009 results using our version of her model
#trying to show how we can nest Bigg's original model in our version

setwd("C:/Users/jones/BoxSync/NDstuff/Dissertation/4/MultiSppModel/")
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
tstep=1:300

qE1Fun=approxfun(x=tstep,y=c(seq(0,by=0.1111,length.out = 113),rep(0,187)))
# qE2Fun=approxfun(x=tstep,y=seq(0,5,length.out = length(tstep)))
# qE1Fun=approxfun(x=tstep,y=rep(1.8, length(tstep)))
qE2Fun=approxfun(x=tstep,y=rep(0, length(tstep)))

# h1Fun=approxfun(x=tstep,y=seq(0,10,length.out = length(tstep)))
# h2Fun=approxfun(x=tstep,y=seq(0,10,length.out = length(tstep)))
h1Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
h2Fun=approxfun(x=tstep,y=rep(8,length(tstep)))

simBiggs<-function(t,y,params){
  A1<-y[1]
  A2<-y[2]
  J1<-y[3]
  J2<-y[4]
  with(as.list(params),{
    dA1dt=-qE1Fun(t)*A1+(s1-1)*A1+s1*J1
    dA2dt=-qE2Fun(t)*A2+(s2-1)*A2+s2*J2
    dJ1dt=-cJ1A1*J1*A1-cJ1J2*J2*J1-(cJ1A2*v1*A2*J1)/(h1Fun(t)+v1+cJ1A2*A2)+f1*A1
    dJ2dt=-cJ2A2*J2*A2-cJ2J1*J1*J2-(cJ2A1*v2*A1*J2)/(h2Fun(t)+v2+cJ2A1*A1)+f2*A2
    return(list(c(dA1dt,dA2dt,dJ1dt,dJ2dt)))
  })
}
p=c(s1=0.1,cJ1A1=0.001,cJ1A2=0.5,cJ1J2=0.000,v1=1,f1=2,
    s2=0.1,cJ2A2=0.001,cJ2A1=0.3,cJ2J1=0.000,v2=1,f2=2)
y0=c(100,10,0,0)
sim=ode(y=y0,times=tstep,func=simBiggs,parms=p)
#sim[nrow(sim),]
plot(sim[,1],sim[,2],type='l',ylim=c(0,max(sim[,2:3],na.rm = T)),col='grey',lwd=2,xlab = "Year",ylab = "pop. size")
lines(sim[,1],sim[,3],col='black',lwd=2)
lines(sim[,1],qE1Fun(sim[,1]),lty=2)
legend("topleft", legend = c("piscivore","planktivore"), col = c('grey','black'),lty = 1,lwd=2,bty = "n")
#lines(sim[,1],sim[,4],lty=2,col='grey')
#lines(sim[,1],sim[,5],lty=2,col='black')


### RECREATING BIGGS HARVEST RESULTS ####

store=data.frame(qEs=seq(0,8,length.out=30),A1=0,A2=0,J1=0,J2=0)
y0=c(100,10,0,0)
for(i in 1:nrow(store)){
  p=c(c(qE1=store$qEs[i],s1=0.5,cJ1A1=0.007,cJ1A2=0.003,cJ1J2=0.001,v1=1,h1=6,f1=1),
      c(qE2=0,s2=0.5,cJ2A2=0.007,cJ2A1=0.003,cJ2J1=0.001,v2=1,h2=6,f2=1))
  sim=ode(y=y0,times=times,func=simBiggs,parms=p)
  store$A1[i]=sim[nrow(sim),2]
  store$A2[i]=sim[nrow(sim),3]
  store$J1[i]=sim[nrow(sim),4]
  store$J2[i]=sim[nrow(sim),5]
}

plot(store$qEs,store$A1,lwd=3,type='l',ylim=c(0,max(store[,2:3])),ylab = "Abundance",xlab = "Harvest (q*E)",col='grey')
lines(store$qEs,store$A2,lwd=3,col='black')
legend("topright",legend = c("sp 1", "sp 2"), col = c("grey","black"),lwd=2,lty = 1,bty="n")


### RECREATING BIGGS RESULTS FOR SHORELINE DEVELOPMENT ####

store=data.frame(hs=seq(0,8,length.out=30),A1=0,A2=0,J1=0,J2=0)
y0=c(100,10,40,40)
for(i in 1:nrow(store)){
  p=c(c(qE1=0,s1=0.5,cJ1A1=0.007,cJ1A2=0.003,cJ1J2=0.001,v1=1,h1=store$hs[i],f1=1),
      c(qE2=0,s2=0.5,cJ2A2=0.007,cJ2A1=0.003,cJ2J1=0.001,v2=1,h2=6,f2=1))
  sim=ode(y=y0,times=times,func=simBiggs,parms=p)
  store$A1[i]=sim[nrow(sim),2]
  store$A2[i]=sim[nrow(sim),3]
  store$J1[i]=sim[nrow(sim),4]
  store$J2[i]=sim[nrow(sim),5]
}

plot(store$hs,store$A1,lwd=3,type='l',ylim=c(0,max(store[,2:3])),ylab = "Abundance",xlab = "Shoreline restoration", col='grey')
lines(store$hs,store$A2,lwd=3,col='black')
legend("topright",legend = c("sp 1", "sp 2"), col = c("grey","black"),lwd=2,lty = 1,bty="n")


