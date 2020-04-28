library(deSolve)
simBiggs<-function(t,y,params){
  A<-y[1]
  F<-y[2]
  J<-y[3]
  with(as.list(params),{
    dAdt=-qE*A+(s-1)*A+s*J
    dFdt=Df*(Fr-F)-cFA*F*A
    dJdt=-cJA*J*A-(cJF*v*F*J)/(h+v+cJF*F)+f*A
    return(list(c(dAdt,dFdt,dJdt)))
  })
}
times=1:500
p=c(qE=1.8,s=0.5,Df=0.1,Fr=100,cFA=0.3,cJA=0.001,cJF=0.5,v=1,h=8,f=2)
y0=c(200,50,10)
sim=ode(y=y0,times=times,func=simBiggs,parms=p)
sim[nrow(sim),]
plot(sim[,1],sim[,2],type='l',ylim=c(0,max(sim[,2:3])))
lines(sim[,1],sim[,3],lty=2)
store=data.frame(qEs=seq(1.5,8,length.out=30),A=0,P=0,J=0,ref=0)
y0=c(20,100,40)
for(i in 1:nrow(store)){
  p=c(qE=store$qEs[i],s=0.5,Df=0.1,Fr=100,cFA=0.3,cJA=0.001,cJF=0.5,v=1,h=8,f=2)
  sim=ode(y=y0,times=times,func=simBiggs,parms=p)
  store$A[i]=sim[nrow(sim),2]
  store$P[i]=sim[nrow(sim),3]
  store$J[i]=sim[nrow(sim),4]
  store$ref[i]=(0.5*1*sim[nrow(sim),3]*sim[nrow(sim),4])/(8+1+0.5*sim[nrow(sim),3])
}
store2=data.frame(qEs=seq(1.5,8,length.out=30),A=0,P=0,J=0,ref=0)
y0=c(200,10,40)
for(i in 1:nrow(store)){
  p=c(qE=store2$qEs[i],s=0.5,Df=0.1,Fr=100,cFA=0.3,cJA=0.001,cJF=0.5,v=1,h=8,f=2)
  sim=ode(y=y0,times=times,func=simBiggs,parms=p)
  store2$A[i]=sim[nrow(sim),2]
  store2$P[i]=sim[nrow(sim),3]
  store2$J[i]=sim[nrow(sim),4]
  store2$ref[i]=(0.5*1*sim[nrow(sim),3]*sim[nrow(sim),4])/(8+1+0.5*sim[nrow(sim),3])
}
par(mar=c(5,5,2,5))
plot(store$qEs,store$A,lwd=3,type='l', ylab = "Abundance", xlab = "Harvest (q*E)")
lines(store$qEs,store$P,lwd=3,col='grey')
lines(store2$qEs,store2$A,lwd=3,lty=3)
lines(store2$qEs,store2$P,lwd=3,col='grey',lty=3)
par(new=T)
plot(store$qEs,store$ref,axes = F, type = 'l',xlab = NA,ylab = NA,col="blue",ylim=c(0,25),lwd=3)
lines(store2$qEs,store2$ref,lty=2,lwd=3,col="blue")
axis(side=4)
mtext(side=4,line=3,'cJF*v*F*J / h+v+cJF*F')

#### ADDING SECOND SPECIES IN ####

simBiggs2<-function(t,y,params){
  A1<-y[1]
  A2<-y[2]
  J1<-y[3]
  J2<-y[4]
  with(as.list(params),{
    dA1dt=-qE1*A1+(s1-1)*A1+s1*J1
    dA2dt=-qE2*A2+(s2-1)*A2+s2*J2
    dJ1dt=-cJ1A1*J1*A1-(cJ1A2*v1*A2*J1)/(h1+v1+cJ1A2*A2)+f1*A1
    dJ2dt=-cJ2A2*J2*A2-(cJ2A1*v2*A1*J2)/(h2+v2+cJ2A1*A1)+f2*A2
    return(list(c(dA1dt,dA2dt,dJ1dt,dJ2dt)))
  })
}
times=1:500
p=c(c(qE1=1.8,s1=0.5,cJ1A1=0.001,cJ1A2=0.001,v1=1,h1=6,f1=2),
    c(qE2=1.8,s2=0.5,cJ2A2=0.001,cJ2A1=0.003,v2=1,h2=10,f2=2))
y0=c(20,20,10,10)
sim=ode(y=y0,times=times,func=simBiggs2,parms=p)
sim[nrow(sim),]
plot(sim[,1],sim[,2],type='l',ylim=c(0,max(sim[,2:5])))
lines(sim[,1],sim[,3],col="blue")
lines(sim[,1],sim[,4],lty=2)
lines(sim[,1],sim[,5],lty=2,col="blue")


### RUNNING OVER A RANGE OF HARVESTS

store=data.frame(qEs=seq(.05,10,length.out=30),A1=0,A2=0,J1=0,J2=0)
y0=c(10,200,40,40)
for(i in 1:nrow(store)){
  p=c(c(qE1=store$qEs[i],s1=0.5,cJ1A1=0.001,cJ1A2=0.001,v1=1,h1=6,f1=2),
      c(qE2=1.8,s2=0.5,cJ2A2=0.001,cJ2A1=0.005,v2=1,h2=10,f2=2))
  sim=ode(y=y0,times=times,func=simBiggs2,parms=p)
  store$A1[i]=sim[nrow(sim),2]
  store$A2[i]=sim[nrow(sim),3]
  store$J1[i]=sim[nrow(sim),4]
  store$J2[i]=sim[nrow(sim),5]
}
plot(store$qEs,store$A1,lwd=2,type='l',ylim=c(0,max(store[,2:3])),ylab = "Abundance",xlab = "q*E")
lines(store$qEs,store$A2,lwd=2,col='grey')
legend("topright",legend = c("species 1", "species 2"), col = c("black","grey"),lwd=2)
store2=data.frame(qEs=seq(.05,10,length.out=30),A1=0,A2=0,J1=0,J2=0)
y0=c(200,10,40,40)
for(i in 1:nrow(store)){
  p=c(c(qE1=store$qEs[i],s1=0.5,cJ1A1=0.001,cJ1A2=0.001,v1=1,h1=6,f1=2),
      c(qE2=1.8,s2=0.5,cJ2A2=0.001,cJ2A1=0.005,v2=1,h2=10,f2=2))
  sim=ode(y=y0,times=times,func=simBiggs2,parms=p)
  store$A1[i]=sim[nrow(sim),2]
  store$A2[i]=sim[nrow(sim),3]
  store$J1[i]=sim[nrow(sim),4]
  store$J2[i]=sim[nrow(sim),5]
}
lines(store$qEs,store$A1,lwd=4,lty=3)
lines(store$qEs,store$A2,lwd=4,col='grey',lty=3)

plot((store$qEs-1.8),store$A1,lwd=2,type='l',ylim=c(0,max(store[,2:3])),ylab = "Abundance",xlab = "Mortality Difference")
lines((store$qEs-1.8),store$A2,lwd=2,col='grey')
legend("topright",legend = c("species 1", "species 2"), col = c("black","grey"),lwd=2)

#### ADDING JUVENILE COMP. TO 2SPP MODEL ####

simBiggs3<-function(t,y,params){
  A1<-y[1]
  A2<-y[2]
  J1<-y[3]
  J2<-y[4]
  with(as.list(params),{
    dA1dt=-qE1*A1+(s1-1)*A1+s1*J1
    dA2dt=-qE2*A2+(s2-1)*A2+s2*J2
    dJ1dt=-cJ1A1*J1*A1-cJ1A2*A2*J1-(cJ1J2*v1*J2*J1)/(h1+v1+cJ1J2*J2)+f1*A1
    dJ2dt=-cJ2A2*J2*A2-cJ2A1*A1*J2-(cJ2J1*v2*J1*J2)/(h2+v2+cJ2J1*J1)+f2*A2
    return(list(c(dA1dt,dA2dt,dJ1dt,dJ2dt)))
  })
}
times=1:500
p=c(c(qE1=1.8,s1=0.5,cJ1A1=0.003,cJ1A2=0.003,cJ1J2=0.001,v1=1,h1=6,f1=2),
    c(qE2=1.8,s2=0.5,cJ2A2=0.003,cJ2A1=0.003,cJ2J1=0.001,v2=1,h2=6,f2=2))
y0=c(21,20,10,10)
sim=ode(y=y0,times=times,func=simBiggs3,parms=p)
sim[nrow(sim),]
plot(sim[,1],sim[,2],type='l',ylim=c(0,max(sim[,2:5])))
lines(sim[,1],sim[,3],col="blue")
lines(sim[,1],sim[,4],lty=2)
lines(sim[,1],sim[,5],lty=2,col="blue")

### RUNNING OVER A RANGE OF HARVESTS

store=data.frame(qEs=seq(.05,8,length.out=30),A1=0,A2=0,J1=0,J2=0)
y0=c(10,200,40,40)
for(i in 1:nrow(store)){
  p=c(c(qE1=store$qEs[i],s1=0.5,cJ1A1=0.002,cJ1A2=0.002,cJ1J2=0.001,v1=1,h1=1,f1=1),
      c(qE2=1.8,s2=0.5,cJ2A2=0.002,cJ2A1=0.002,cJ2J1=0.001,v2=1,h2=20,f2=1))
  sim=ode(y=y0,times=times,func=simBiggs3,parms=p)
  store$A1[i]=sim[nrow(sim),2]
  store$A2[i]=sim[nrow(sim),3]
  store$J1[i]=sim[nrow(sim),4]
  store$J2[i]=sim[nrow(sim),5]
}
plot(store$qEs,store$A1,lwd=3,type='l',ylim=c(0,max(store[,2:3])),ylab = "Abundance",xlab = "q*E")
lines(store$qEs,store$A2,lwd=3,col='grey')
store2=data.frame(qEs=seq(.05,8,length.out=30),A1=0,A2=0,J1=0,J2=0)
y0=c(200,10,40,40)
for(i in 1:nrow(store)){
  p=c(c(qE1=store2$qEs[i],s1=0.5,cJ1A1=0.002,cJ1A2=0.002,cJ1J2=0.001,v1=1,h1=6,f1=1),
      c(qE2=1.8,s2=0.5,cJ2A2=0.002,cJ2A1=0.002,cJ2J1=0.001,v2=1,h2=6,f2=1))
  sim=ode(y=y0,times=times,func=simBiggs3,parms=p)
  store2$A1[i]=sim[nrow(sim),2]
  store2$A2[i]=sim[nrow(sim),3]
  store2$J1[i]=sim[nrow(sim),4]
  store2$J2[i]=sim[nrow(sim),5]
}
lines(store2$qEs,store2$A1,lwd=4,lty=3)
lines(store2$qEs,store2$A2,lwd=4,col='grey',lty=3)
legend("topright",legend = c("sp 1", "sp 2"), col = c("black","grey"),lwd=2)

#flipping location of juvenile competition and interspecific competition in the model
simBiggs4<-function(t,y,params){
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
p=c(c(qE1=1.8,s1=0.5,cJ1A1=0.003,cJ1A2=0.003,cJ1J2=0.001,v1=1,h1=6,f1=2),
    c(qE2=1.8,s2=0.5,cJ2A2=0.003,cJ2A1=0.003,cJ2J1=0.001,v2=1,h2=6,f2=2))
y0=c(21,20,10,10)
sim=ode(y=y0,times=times,func=simBiggs4,parms=p)
sim[nrow(sim),]
plot(sim[,1],sim[,2],type='l',ylim=c(0,max(sim[,2:5])))
lines(sim[,1],sim[,3],col="blue")
lines(sim[,1],sim[,4],lty=2)
lines(sim[,1],sim[,5],lty=2,col="blue")

### RUNNING OVER A RANGE OF HARVESTS

store=data.frame(qEs=seq(.05,8,length.out=30),A1=0,A2=0,J1=0,J2=0)
y0=c(10,200,40,40)
for(i in 1:nrow(store)){
  p=c(c(qE1=store$qEs[i],s1=0.5,cJ1A1=0.002,cJ1A2=0.002,cJ1J2=0.001,v1=1,h1=6,f1=1),
      c(qE2=1.8,s2=0.5,cJ2A2=0.002,cJ2A1=0.002,cJ2J1=0.001,v2=1,h2=6,f2=1))
  sim=ode(y=y0,times=times,func=simBiggs4,parms=p)
  store$A1[i]=sim[nrow(sim),2]
  store$A2[i]=sim[nrow(sim),3]
  store$J1[i]=sim[nrow(sim),4]
  store$J2[i]=sim[nrow(sim),5]
}
plot(store$qEs,store$A1,lwd=3,type='l',ylim=c(0,max(store[,2:3])),ylab = "Abundance",xlab = "Harvest (q*E)")
lines(store$qEs,store$A2,lwd=3,col='grey')
store2=data.frame(qEs=seq(.05,8,length.out=30),A1=0,A2=0,J1=0,J2=0)
y0=c(200,10,40,40)
for(i in 1:nrow(store)){
  p=c(c(qE1=store2$qEs[i],s1=0.5,cJ1A1=0.002,cJ1A2=0.002,cJ1J2=0.001,v1=1,h1=6,f1=1),
      c(qE2=1.8,s2=0.5,cJ2A2=0.002,cJ2A1=0.002,cJ2J1=0.001,v2=1,h2=6,f2=1))
  sim=ode(y=y0,times=times,func=simBiggs4,parms=p)
  store2$A1[i]=sim[nrow(sim),2]
  store2$A2[i]=sim[nrow(sim),3]
  store2$J1[i]=sim[nrow(sim),4]
  store2$J2[i]=sim[nrow(sim),5]
}
lines(store2$qEs,store2$A1,lwd=3,lty=3)
lines(store2$qEs,store2$A2,lwd=3,col='grey',lty=3)
legend("topright",legend = c("sp 1", "sp 2"), col = c("black","grey"),lwd=2)

#taking cannibalism out of the model
simBiggs5<-function(t,y,params){
  A1<-y[1]
  A2<-y[2]
  J1<-y[3]
  J2<-y[4]
  with(as.list(params),{
    dA1dt=-qE1*A1+(s1-1)*A1+s1*J1
    dA2dt=-qE2*A2+(s2-1)*A2+s2*J2
    dJ1dt=-cJ1J2*J2*J1-(cJ1A2*v1*A2*J1)/(h1+v1+cJ1A2*A2)+f1*A1
    dJ2dt=-cJ2J1*J1*J2-(cJ2A1*v2*A1*J2)/(h2+v2+cJ2A1*A1)+f2*A2
    return(list(c(dA1dt,dA2dt,dJ1dt,dJ2dt)))
  })
}
times=1:1000
p=c(c(qE1=1.8,s1=0.2,cJ1A1=0.003,cJ1A2=0.06,cJ1J2=0.001,v1=1,h1=8,f1=2),
    c(qE2=1.8,s2=0.2,cJ2A2=0.003,cJ2A1=0.06,cJ2J1=0.001,v2=1,h2=8,f2=2))
y0=c(30,20,10,10)
sim=ode(y=y0,times=times,func=simBiggs5,parms=p)
sim[nrow(sim),]
plot(sim[,1],sim[,2],type='l',ylim=c(0,max(sim[,2:5])))
lines(sim[,1],sim[,3],col="blue")
lines(sim[,1],sim[,4],lty=2)
lines(sim[,1],sim[,5],lty=2,col="blue")

### RUNNING OVER A RANGE OF HARVESTS
# 4.24.20- haven't tried this yet the model is behaving weridly without cannibalism and I'm not sure why
store=data.frame(qEs=seq(.05,8,length.out=30),A1=0,A2=0,J1=0,J2=0)
y0=c(10,200,40,40)
for(i in 1:nrow(store)){
  p=c(c(qE1=store$qEs[i],s1=0.5,cJ1A1=0.002,cJ1A2=0.002,cJ1J2=0.001,v1=1,h1=8,f1=1),
      c(qE2=1.8,s2=0.5,cJ2A2=0.002,cJ2A1=0.002,cJ2J1=0.001,v2=1,h2=8,f2=1))
  sim=ode(y=y0,times=times,func=simBiggs4,parms=p)
  store$A1[i]=sim[nrow(sim),2]
  store$A2[i]=sim[nrow(sim),3]
  store$J1[i]=sim[nrow(sim),4]
  store$J2[i]=sim[nrow(sim),5]
}
plot(store$qEs,store$A1,lwd=3,type='l',ylim=c(0,max(store[,2:3])),ylab = "Abundance",xlab = "q*E")
lines(store$qEs,store$A2,lwd=3,col='grey')
store2=data.frame(qEs=seq(.05,8,length.out=30),A1=0,A2=0,J1=0,J2=0)
y0=c(200,10,40,40)
for(i in 1:nrow(store)){
  p=c(c(qE1=store2$qEs[i],s1=0.5,cJ1A1=0.002,cJ1A2=0.002,cJ1J2=0.001,v1=1,h1=6,f1=1),
      c(qE2=1.8,s2=0.5,cJ2A2=0.002,cJ2A1=0.002,cJ2J1=0.001,v2=1,h2=6,f2=1))
  sim=ode(y=y0,times=times,func=simBiggs4,parms=p)
  store2$A1[i]=sim[nrow(sim),2]
  store2$A2[i]=sim[nrow(sim),3]
  store2$J1[i]=sim[nrow(sim),4]
  store2$J2[i]=sim[nrow(sim),5]
}
lines(store2$qEs,store2$A1,lwd=4,lty=3)
lines(store2$qEs,store2$A2,lwd=4,col='grey',lty=3)
legend("topright",legend = c("sp 1", "sp 2"), col = c("black","grey"),lwd=2)
