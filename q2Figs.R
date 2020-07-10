# trying to make some heat maps to explore how stocking can maintain sp1 dominance
library(deSolve)
library(ggplot2)
rm(list=ls())

#approx funs for changing parms through time
qE1Fun=approxfun(x=tstep,y=c(seq(0,by=0.1111,length.out = 113),rep(0,187)))
# qE2Fun=approxfun(x=tstep,y=seq(0,5,length.out = length(tstep)))
# qE1Fun=approxfun(x=tstep,y=rep(1.8, length(tstep)))
qE2Fun=approxfun(x=tstep,y=rep(0, length(tstep)))

# h1Fun=approxfun(x=tstep,y=seq(0,10,length.out = length(tstep)))
# h2Fun=approxfun(x=tstep,y=seq(0,10,length.out = length(tstep)))
h1Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
h2Fun=approxfun(x=tstep,y=rep(8,length(tstep)))

st1Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
st2Fun=approxfun(x=tstep,y=rep(0,length(tstep)))

#matrix to hold output, starting with different harvest levels on each species

qEs=seq(0,8,length.out = 30)
sto=seq(0,200, length.out = 30)
df=expand.grid(X=qEs,Y=sto)
df$A1=numeric(900)
df$A2=numeric(900)
df$J1=numeric(900)
df$J2=numeric(900)

source('q2Func.R')

for(i in 1:nrow(df)){
    tstep=1:300
    qE1Fun=approxfun(x=tstep,y=rep(df$X[i], length(tstep)))
    qE2Fun=approxfun(x=tstep,y=rep(0, length(tstep)))
    st1Fun=approxfun(x=tstep,y=rep(df$Y[i],length(tstep)))
    st2Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
    p=c(s1=0.1,cJ1A1=0.002,cJ1A2=0.5,cJ1J2=0.003,v1=1,f1=2,
        s2=0.1,cJ2A2=0.002,cJ2A1=0.3,cJ2J1=0.003,v2=1,f2=2)
    y0=c(10,10,0,0)
    sim=ode(y=y0,times=tstep,func=simBiggsQ2,parms=p)
    df$A1[i]=sim[nrow(sim)-1,2]
    df$A2[i]=sim[nrow(sim)-1,3]
    df$J1[i]=sim[nrow(sim)-1,4]
    df$J2[i]=sim[nrow(sim)-1,5]
}

vz=ggplot(data = df)+theme_classic()
vz+geom_tile(aes(x=df$X, y=df$Y, fill=df$A1))+scale_fill_gradient(low ="blue", high = "yellow")
vz+geom_tile(aes(x=df$X, y=df$Y, fill=df$A2))+scale_fill_gradient(low ="blue", high = "yellow")


#### OPTIM FUNC ####
source('q2Func.R')

optA1=function(parms){
  tstep=1:300
  qE1Fun=approxfun(x=tstep,y=rep(parms[1], length(tstep)))
  qE2Fun=approxfun(x=tstep,y=rep(parms[2], length(tstep)))
  st1Fun=approxfun(x=tstep,y=rep(parms[3],length(tstep)))
  st2Fun=approxfun(x=tstep,y=rep(parms[4],length(tstep)))
  p=c(s1=0.1,cJ1A1=0.002,cJ1A2=0.5,cJ1J2=0.003,v1=1,f1=2,
      s2=0.1,cJ2A2=0.002,cJ2A1=0.3,cJ2J1=0.003,v2=1,f2=2)
  y0=c(10,10,0,0)
  sim=ode(y=y0,times=tstep,func=simBiggsQ2,parms=p)
  return(sim[nrow(sim)-1,2]) #A1 abundance, what I want to maximize
}

parms=c(1.8,0,200,0)
out=optim(par=parms, optA1,control = list(fnscale=-1))
out
