#trying out a new way of coding up the multispp model 
#simplifying by using a lotka volterra model to start with
#modified version of the model on pg 245-246 in A Primer of Ecology with R book

rm(list=ls())
msm=function(years=100,
              alpha=matrix(c(0.01, 0.02, 0.02, 0.01), ncol = 2, byrow = T),
              rd=c(1,1), 
              q=.1, 
              e1=rep(0,times=years),
              e2=rep(0,times=years),
              start=c(50,25)){
  #storing output
  n1=numeric(years)
  n2=numeric(years)
  n1[1]=start[1]
  n2[1]=start[2]
  for(i in 1:years){
    n1[i+1]=n1[i]+rd[1]*n1[i]*(1-alpha[1,1]*n1[i]-alpha[1,2]*n2[i])-q*e1[i]*n1[i]
    n2[i+1]=n2[i]+rd[2]*n2[i]*(1-alpha[2,1]*n1[i]-alpha[2,2]*n2[i])-q*e2[i]*n2[i]
  }
  return(cbind(n1,n2))
}

#matrix of alphas to be used in the model 
alphs=matrix(c(0.01, 0.01, 0.01, 0.01), ncol = 2, byrow = T)
#basic run of the model without harvesting - here dominant species is controlled by priority effect
base=msm(start=c(50,40), alpha = alphs)
matplot(0:100, base, type = "l", col=1, ylim = c(0, 110))
abline(h=1/0.01, lty=3) #carrying capacity

#exploring the effects of harvesting on flipping the system. I modify harvest by increasing or decreasing effort for each species
test=msm(years=200, start = c(50, 10),
          alpha = alphs,
          e1=c(rep(0,100), rep(2, 25), rep(0,75)),
          e2=c(rep(0,200)))

matplot(0:200,test, type = "l", col = 1, ylim = c(0,110))
abline(h=1/0.01, lty=3, col="grey") # carrying capacity
abline(v=100, lty=3, col="grey") # start of harvesting
abline(v=130, lty=3, col="grey") # end of harvesting

