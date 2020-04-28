## trying one of the models presented in Holt and Polis 1997, pg 755.
library(deSolve)
rm(list=ls())

hpMod=function(t,y,p){
  #unpack states & parms
  A=y[1]
  B=y[2]
  with(as.list(p),{
    dA.dt=A*((ba*ea*I)/(ea*A+eb*B)+beta*alpha*B-ma)
    dB.dt=B*((bb*eb*I)/(ea*A+eb*B)-alpha*A-mb)
    return(list(c(dA.dt, dB.dt)))
  })
}

y0=c(100,100)
t=1:1000
p=c(I=10,bb=1,eb=.75,ba=.5,ea=.5,beta=2,alpha=.01,mb=.3,ma=.3)

out=ode(y = y0, times = t, func = hpMod, parms = p)
plot(out[,1], out[,2], type="l", ylim = c(-5,max(out[,2:3])));lines(out[,1], out[,3], lty=2);legend("topright", legend = c("A","B"),lty = 1:2)

#parm definitions
# b1 & b2 - convert resource consumption to reproduction
# e1 & e2 - efficiency at finding the resource (or their ability to compete for it)
# I - rate of input for the shared resource
# beta - benefit of consuming prey for the predator
# alpha - mortality inflicted on prey by predator
# m1 & m2 - natural mortality rate


#run model over a range of mortalities
store=data.frame(mort=seq(.1,.9, length.out = 15), A=0, B=0)
y0=c(100,10)
for(i in 1:nrow(store)){
  p=c(I=10,bb=1,eb=.75,ba=.5,ea=.5,beta=2,alpha=.01,mb=.25,ma=store$mort[i])
  out=ode(y = y0, times = t, func = hpMod, parms = p)
  store$A[i]=out[nrow(out),1]
  store$B[i]=out[nrow(out),2]
}
store2=data.frame(mort=seq(.1,.9, length.out = 15), A=0, B=0)
y0=c(10,100)
for(i in 1:nrow(store2)){
  p=c(I=10,bb=1,eb=.75,ba=.5,ea=.5,beta=2,alpha=.01,mb=.25,ma=store2$mort[i])
  out=ode(y = y0, times = t, func = hpMod, parms = p)
  store2$A[i]=out[nrow(out),1]
  store2$B[i]=out[nrow(out),2]
}
plot(store$mort,store$A,type = "l",lwd=2,ylim = c(0,max(store[,2:3])))
lines(store$mort,store$B,col="grey",lwd=2)
lines(store2$mort,store2$A,lwd=3,lty=2)
lines(store2$mort,store2$B,lwd=3,lty=2)
