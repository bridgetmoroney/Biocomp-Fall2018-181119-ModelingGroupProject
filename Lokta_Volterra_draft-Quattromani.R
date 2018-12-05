rm(list=ls())
setwd("/Users/ta4ha/Documents/Biocomputing")

# Loading packages
library(ggplot2)
library(deSolve)

# Lokta-Voltera Model
LVsim<-function(t,y,p){  #assigning values to function variables
  H=y[1] # prey initial value
  P=y[2] # predator initial value
  
  b=p[1] # prey birth rate
  a=p[2] # predator attack rate
  e=p[3] # conversion efficiency of prey to predators
  s=p[4] # predator death rate
  
  #defining consumer-resource differential equations
  dHdt=b*H-(a*P*H) # prey population
  dPdt=(e*a*P*H)-s*P # predator population
  
  return(list(c(dHdt,dPdt)))
}

# Initial Simulation
times= seq(0, 100, by=0.01) # Run simulation from time unit 0 to 100, over a 0.1 time step
y0=c(25,5) # Ho=25  Po=5
params1=c(0.5,0.02,0.1,0.2) # b=0.5  a=0.02  e=0.1  s=0.2
LVsim1=ode(y=y0,times=times,func=LVsim,parms=params1) # Use differential equations to simulate defined functions
LVout1=data.frame(time=LVsim1[,1],population=LVsim1[,2],animals=LVsim1[,3]) #create dataframe
# Plot data
ggplot(LVout1,aes(x=time,y=population))+geom_line()+geom_line(data=LVout1,mapping=aes(x=time,y=animals),col='red')+theme_classic()+ggtitle("Predator and Prey Populations Over Time")

# Increase b (prey birth rate)
times= seq(0, 100, by=0.01)
y0=c(25,5)
params2=c(1,0.02,0.1,0.2)
LVsim2=ode(y=y0,times=times,func=LVsim,parms=params2)
LVout2=data.frame(time=LVsim2[,1],population=LVsim2[,2],animals=LVsim2[,3])
ggplot(LVout2,aes(x=time,y=population))+geom_line()+geom_line(data=LVout2,mapping=aes(x=time,y=animals),col='red')+theme_classic()+ggtitle("Predator and Prey Populations Over Time")

# Decrease b
times= seq(0, 100, by=0.01)
y0=c(25,5)
params3=c(0.25,0.02,0.1,0.2)
LVsim3=ode(y=y0,times=times,func=LVsim,parms=params3)
LVout3=data.frame(time=LVsim3[,1],population=LVsim3[,2],animals=LVsim3[,3])
ggplot(LVout3,aes(x=time,y=population))+geom_line()+geom_line(data=LVout3,mapping=aes(x=time,y=animals),col='red')+theme_classic()+ggtitle("Predator and Prey Populations Over Time")

# Increase a (predator attack rate)
times= seq(0, 100, by=0.01)
y0=c(25,5)
params4=c(0.5,0.04,0.1,0.2)
LVsim4=ode(y=y0,times=times,func=LVsim,parms=params4)
LVout4=data.frame(time=LVsim4[,1],population=LVsim4[,2],animals=LVsim4[,3])
ggplot(LVout4,aes(x=time,y=population))+geom_line()+geom_line(data=LVout4,mapping=aes(x=time,y=animals),col='red')+theme_classic()+ggtitle("Predator and Prey Populations Over Time")

# Decrease a
times= seq(0, 100, by=0.01)
y0=c(25,5)
params5=c(0.5,0.01,0.1,0.2)
LVsim5=ode(y=y0,times=times,func=LVsim,parms=params5)
LVout5=data.frame(time=LVsim5[,1],population=LVsim5[,2],animals=LVsim5[,3])
ggplot(LVout5,aes(x=time,y=population))+geom_line()+geom_line(data=LVout5,mapping=aes(x=time,y=animals),col='red')+theme_classic()+ggtitle("Predator and Prey Populations Over Time")

# Increase e (conversion efficiency of prey to predators)
times= seq(0, 100, by=0.01)
y0=c(25,5)
params6=c(0.5,0.02,0.2,0.2)
LVsim6=ode(y=y0,times=times,func=LVsim,parms=params6)
LVout6=data.frame(time=LVsim6[,1],population=LVsim6[,2],animals=LVsim6[,3])
ggplot(LVout6,aes(x=time,y=population))+geom_line()+geom_line(data=LVout6,mapping=aes(x=time,y=animals),col='red')+theme_classic()+ggtitle("Predator and Prey Populations Over Time")

# Decrease e
times= seq(0, 100, by=0.01)
y0=c(25,5)
params7=c(0.5,0.02,0.05,0.2)
LVsim7=ode(y=y0,times=times,func=LVsim,parms=params7)
LVout7=data.frame(time=LVsim7[,1],population=LVsim7[,2],animals=LVsim7[,3])
ggplot(LVout7,aes(x=time,y=population))+geom_line()+geom_line(data=LVout7,mapping=aes(x=time,y=animals),col='red')+theme_classic()+ggtitle("Predator and Prey Populations Over Time")

# Increase s (predator death rate)
times= seq(0, 100, by=0.01)
y0=c(25,5)
params8=c(0.5,0.02,0.1,0.4)
LVsim8=ode(y=y0,times=times,func=LVsim,parms=params8)
LVout8=data.frame(time=LVsim8[,1],population=LVsim8[,2],animals=LVsim8[,3])
ggplot(LVout8,aes(x=time,y=population))+geom_line()+geom_line(data=LVout8,mapping=aes(x=time,y=animals),col='red')+theme_classic()+ggtitle("Predator and Prey Populations Over Time")

# Decrease s
times= seq(0, 100, by=0.01)
y0=c(25,5)
params9=c(0.5,0.02,0.1,0.1)
LVsim9=ode(y=y0,times=times,func=LVsim,parms=params9)
LVout9=data.frame(time=LVsim9[,1],population=LVsim9[,2],animals=LVsim9[,3])
ggplot(LVout9,aes(x=time,y=population))+geom_line()+geom_line(data=LVout9,mapping=aes(x=time,y=animals),col='red')+theme_classic()+ggtitle("Predator and Prey Populations Over Time")