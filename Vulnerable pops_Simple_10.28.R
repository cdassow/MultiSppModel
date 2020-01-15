

options(scipen=999)
Years <- 150 # eventually want years and days

##### WALLEYE ######
Nw <- 30 #initial population size

rw <-  vector()
rw <- c(rw, 0.3:150)
for (j in 1:Years){
  rw[j+1]=rw[j]-0.001
}

# pop growth rate, decreases every year beacuse CLIMATE CHANGE or something
kw <- 400 #carrying capacity
qw <- 0.09 # catchability of walleye
Ew <- 1 #effort for walleye
Mw <- 0.02 # natural mortality

#this is just to look at how pop model is working
for (i in 1:Years){
  Nw[i+1] <- Nw[i] + (rw[i]*Nw[i]*(1-(Nw[i]/kw)))-(qw*Ew*Nw[i])-(Nw[i]*Mw)
}
matplot(Nw,type="l", lty=1)  


##### BASS #######
Nb <- 70 #initial population size
rb <- 0.3 # pop growth rate
kb <- 500 
qb <- 0.12 # catchability 
Eb <- 1 #effort (CPUE)
Mb <- 0.02 # natural mortality

#this is just to look at how pop model is working
for (i in 1:Years){
  Nb[i+1] <- Nb[i] + (rb*Nb[i]*(1-(Nb[i]/kb)))-(qb*Eb*Nb[i])-(Nb[i]*Mb)
}
matplot(Nb,type="l", lty=1)  


######### Angler fish/fishing choice #######
## Do we need to do this for each individual angler? 
# First we set up a expected catch vector
expCw <-  vector()
expCw <- c(expCw, 1:Years) 
expCb <- vector()
expCb <- c(expCb, 1:Years)
# Now we need to develop how this updates. - we need actual catch to do this. 
actCw <-  vector()
actCw <- c(actCw, 1:Years) 
actCb <- vector()
actCb <- c(actCb, 1:Years)
# as of right now, expected catch rate = based on previous actual with a random mixed in
#I think we need to set up a dummy for the first round
Uw <- vector()
Uw <- c(Uw, 0.7:Years) 
Ub <-  vector()
Ub <- c(Ub, 0.3:Years) 
Usit = 5
ProAngw <-  vector() #starting distribution of anglers
ProAngw <- c(ProAngw, 0:Years) 
ProAngb <-  vector()
ProAngb <- c(ProAngb, 0:Years) 
ProAngsit <-  vector()
ProAngsit <- c(ProAngsit, 0:Years) 
# preference determination:
pi <- .7 #preference for walleye (.5 = equal preference, 1 = only walleye)

####### model and figs ##############
for (i in 1:Years){
  actCw[i] = qw * Ew * Nw[i] #actual catch (fishing mortality) for time 1
  actCb[i] = qb * Eb * Nb[i] 
  expCw[i] = actCw[i]  + rnorm(1,mean=0,sd=4)# expected catch for time 1
  expCb[i] = actCb[i] + rnorm(1,mean=0,sd=4)
  Uw[i] = pi*expCw[i] #weight for walleye = 7
  Ub[i] = (1-pi)*expCb[i] #weight for LMB
  ProAngw[i] = Uw[i]/(Uw[i]+Ub[i]+Usit)
  ProAngb[i] = Ub[i]/(Uw[i]+Ub[i]+Usit)
  ProAngsit[i] = Usit/(Uw[i]+Ub[i]+Usit)
  Nb[i+1] = Nb[i] + (rb*Nb[i]*(1-(Nb[i]/kb)))-(actCb[i]*ProAngb[i])-(Nb[i]*Mb) #biomass at time step 2
  Nw[i+1] = Nw[i] + (rw[i]*Nw[i]*(1-(Nw[i]/kw)))-(actCw[i]*ProAngw[i])-(Nw[i]*Mw)
}

matplot(ProAngb,type='l',lty=1, ylim=c(0,1))
matplot(ProAngw, type='l',lty=1, ylim=c(0,1))
matplot(ProAngsit, type='l',lty=1, ylim=c(0,1))
matplot(Nb,type='l',lty=1)
matplot(Nw,type='l',lty=1)

###### TO DO #########
# make sure the utility function is reasonable
# size structured pop model w/ realistic parameters
# figure out how to make the expected catch based on previous exp and not just random







