# some code to explore my version of the Biggs model and see how robust it is to changes in parameter values.
#plan is to run model over a range of different values for a few different parameters using the same two loop structure I used earlier to see if hyperstability occurs

#parms to loop through that change their hyperstability effect depending on magnitute of difference (f,s,cannibalism(cJJ))
library(deSolve)
rm(list = ls())
simBiggs3<-function(t,y,params){
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
#### FECUNDITY LOOP ####
times=1:500

fs=data.frame(fs=seq(0.5,10,length.out = 30),hyst=0) # 1=yes, 0=no
for(f in 1:nrow(fs)){
  
  store=data.frame(qEs=seq(.05,8,length.out=30),A1=0,A2=0,J1=0,J2=0)
  y0=c(10,200,40,40)
  for(i in 1:nrow(store)){
    p=c(c(qE1=store$qEs[i],s1=0.5,cJ1A1=0.002,cJ1A2=0.002,cJ1J2=0.001,v1=1,h1=1,f1=fs$fs[f]),
        c(qE2=1.8,s2=0.5,cJ2A2=0.002,cJ2A1=0.002,cJ2J1=0.001,v2=1,h2=20,f2=1))
    sim=ode(y=y0,times=times,func=simBiggs3,parms=p)
    store$A1[i]=sim[nrow(sim),2]
    store$A2[i]=sim[nrow(sim),3]
    #store$J1[i]=sim[nrow(sim),4]
    #store$J2[i]=sim[nrow(sim),5]
  }
  store2=data.frame(qEs=seq(.05,8,length.out=30),A1=0,A2=0,J1=0,J2=0)
  y0=c(200,10,40,40)
  for(i in 1:nrow(store)){
    p=c(c(qE1=store2$qEs[i],s1=0.5,cJ1A1=0.002,cJ1A2=0.002,cJ1J2=0.001,v1=1,h1=6,f1=fs$fs[f]),
        c(qE2=1.8,s2=0.5,cJ2A2=0.002,cJ2A1=0.002,cJ2J1=0.001,v2=1,h2=6,f2=1))
    sim=ode(y=y0,times=times,func=simBiggs3,parms=p)
    store2$A1[i]=sim[nrow(sim),2]
    store2$A2[i]=sim[nrow(sim),3]
    #store2$J1[i]=sim[nrow(sim),4]
    #store2$J2[i]=sim[nrow(sim),5]
  }
  if(any(abs(store$A2-store2$A2)>1) & any(abs(store$A1-store2$A1)>1)){fs$hyst[f]=1}
  
}
plot(fs$fs,fs$hyst,pch=16,xlab = "fecundity",ylab = "Hysteresis Present?")
plot(fs$fs-1,fs$hyst,pch=16,xlab = "Diff fecundity",ylab = "Hysteresis Present?")

#### SURVIVAL LOOP ####
times=1:500

ss=data.frame(ss=seq(0.05,.9,length.out = 30),hyst=0) # 1=yes, 0=no
for(f in 1:nrow(ss)){
  
  store=data.frame(qEs=seq(.05,8,length.out=30),A1=0,A2=0,J1=0,J2=0)
  y0=c(10,200,40,40)
  for(i in 1:nrow(store)){
    p=c(c(qE1=store$qEs[i],s1=ss$ss[f],cJ1A1=0.002,cJ1A2=0.002,cJ1J2=0.001,v1=1,h1=1,f1=1),
        c(qE2=1.8,s2=0.5,cJ2A2=0.002,cJ2A1=0.002,cJ2J1=0.001,v2=1,h2=20,f2=1))
    sim=ode(y=y0,times=times,func=simBiggs3,parms=p)
    store$A1[i]=sim[nrow(sim),2]
    store$A2[i]=sim[nrow(sim),3]
    #store$J1[i]=sim[nrow(sim),4]
    #store$J2[i]=sim[nrow(sim),5]
  }
  store2=data.frame(qEs=seq(.05,8,length.out=30),A1=0,A2=0,J1=0,J2=0)
  y0=c(200,10,40,40)
  for(i in 1:nrow(store)){
    p=c(c(qE1=store2$qEs[i],s1=ss$ss[f],cJ1A1=0.002,cJ1A2=0.002,cJ1J2=0.001,v1=1,h1=6,f1=1),
        c(qE2=1.8,s2=0.5,cJ2A2=0.002,cJ2A1=0.002,cJ2J1=0.001,v2=1,h2=6,f2=1))
    sim=ode(y=y0,times=times,func=simBiggs3,parms=p)
    store2$A1[i]=sim[nrow(sim),2]
    store2$A2[i]=sim[nrow(sim),3]
    #store2$J1[i]=sim[nrow(sim),4]
    #store2$J2[i]=sim[nrow(sim),5]
  }
  if(any(abs(store$A2-store2$A2)>1) & any(abs(store$A1-store2$A1)>1)){ss$hyst[f]=1}
  
}
plot(ss$ss,ss$hyst,pch=16,xlab = "Survival",ylab = "Hysteresis Present?")
plot(ss$ss-0.5,ss$hyst,pch=16,xlab = "Diff Survival",ylab = "Hysteresis Present?")

#### CANNIBALISM LOOP ####
times=1:500

can=data.frame(can=seq(0.05,.9,length.out = 30),hyst=0) # 1=yes, 0=no
for(f in 1:nrow(can)){
  
  store=data.frame(qEs=seq(.05,8,length.out=30),A1=0,A2=0,J1=0,J2=0)
  y0=c(10,200,40,40)
  for(i in 1:nrow(store)){
    p=c(c(qE1=store$qEs[i],s1=0.5,cJ1A1=can$can[f],cJ1A2=0.002,cJ1J2=0.001,v1=1,h1=1,f1=1),
        c(qE2=1.8,s2=0.5,cJ2A2=0.002,cJ2A1=0.002,cJ2J1=0.001,v2=1,h2=20,f2=1))
    sim=ode(y=y0,times=times,func=simBiggs3,parms=p)
    store$A1[i]=sim[nrow(sim),2]
    store$A2[i]=sim[nrow(sim),3]
    #store$J1[i]=sim[nrow(sim),4]
    #store$J2[i]=sim[nrow(sim),5]
  }
  store2=data.frame(qEs=seq(.05,8,length.out=30),A1=0,A2=0,J1=0,J2=0)
  y0=c(200,10,40,40)
  for(i in 1:nrow(store)){
    p=c(c(qE1=store2$qEs[i],s1=0.5,cJ1A1=can$can[f],cJ1A2=0.002,cJ1J2=0.001,v1=1,h1=6,f1=1),
        c(qE2=1.8,s2=0.5,cJ2A2=0.002,cJ2A1=0.002,cJ2J1=0.001,v2=1,h2=6,f2=1))
    sim=ode(y=y0,times=times,func=simBiggs3,parms=p)
    store2$A1[i]=sim[nrow(sim),2]
    store2$A2[i]=sim[nrow(sim),3]
    #store2$J1[i]=sim[nrow(sim),4]
    #store2$J2[i]=sim[nrow(sim),5]
  }
  if(any(abs(store$A2-store2$A2)>1) & any(abs(store$A1-store2$A1)>1)){can$hyst[f]=1}
  
}
plot(can$can,can$hyst,pch=16,xlab = "Cannibalism Sp1",ylab = "Hysteresis Present?")
plot(can$can-0.002,can$hyst,pch=16,xlab = "Diff in Cannibalism",ylab = "Hysteresis Present?")
