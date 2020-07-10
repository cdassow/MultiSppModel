#script to hold the q2 function along so we can just source it in Rmarkdown 
#parameter definitions
# qE - harvest, species specific, this is what is controlled by 'regulations'
# s - survial, species specific
# cJA - effect of adults of a given species on juveniles of a given species (cover cannibalism or interspecific predation)
# cJJ - effect of juveniles of one species on juveniles of the other (can be predation or competition)
# h - rate at which juveniles leave foraging arena for refuge, species specific
# v - rate at which juveniles enter foraging arena from refuge, species specific
# f - fecundity, species specific
# stock1 - annual stocked num spp1
# stock2 - annual stocked num spp2

simBiggsQ2<-function(t,y,params){
  A1<-y[1]
  A2<-y[2]
  J1<-y[3]
  J2<-y[4]
  with(as.list(params),{
    dA1dt=-qE1Fun(t)*A1+(s1-1)*A1+s1*J1
    dA2dt=-qE2Fun(t)*A2+(s2-1)*A2+s2*J2
    dJ1dt=-cJ1A1*J1*A1-cJ1J2*J2*J1-(cJ1A2*v1*A2*J1)/(h1Fun(t)+v1+cJ1A2*A2)+f1*A1+st1Fun(t)
    dJ2dt=-cJ2A2*J2*A2-cJ2J1*J1*J2-(cJ2A1*v2*A1*J2)/(h2Fun(t)+v2+cJ2A1*A1)+f2*A2+st2Fun(t)
    return(list(c(dA1dt,dA2dt,dJ1dt,dJ2dt)))
  })
}