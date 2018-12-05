## Intro to Biocomputing Final Project: Dynamic Modeling
##
## Jacquie Whalen, Bridget Moroney, Andrew Quattromani
##
## 12/7/18


##Part 1: Lotka Volterra model simulations

library(deSolve) #loading packages
library(ggplot2) #loading packages
library(reshape2) #loading packages
rm(list=ls()) #clearing global environment



lotka_sim<-function(t,y,p){ #creating function for LV model
  H=y[1] #unpacking state variables and parameters, Herbivore H population state variable
  P=y[2] #Predator P population state variable
  
  b=p[1] #prey birth rate parameter
  a=p[2] #predator attack rate parameter
  e=p[3] #conversion efficiency parameter
  s=p[4] #predator death rate parameter
  
  dHdt=(b*H)-(a*P*H) #differential equation for change in prey population
  dPdt=(e*a*P*H)-(s*P) #differential equation for change in predator population
  return(list(c(dHdt, dPdt))) #returning a list of the population changes for predator and prey
}

i_params=c(0.5,0.02,0.1,0.2) #initial parameters
times=seq(from=1, to=100, by=0.1) #time period is from 1 to 100 by increment of 0.1
y0=c(25,5) #initial populations of prey (25) and predator (5)

i_sim=ode(y=y0, times=times, func = lotka_sim, parms = i_params) #simulates model using ODE function based on the parameters, LV function, initial populations, and time period
i_out=data.frame(time=i_sim[,1], H=i_sim[,2], P=i_sim[,3]) #creates a dataframe of time step and the populations of H and P
i_output=melt(i_out, id.vars="time", variable.name="Group", value.name="Density") #melts the dataframe to put it into long format to be used for graphing
ggplot(i_output, aes(x=time, y=Density))+geom_line(aes(color=Group))+theme_classic()+ylab("Population Density")+ggtitle("Initial parameters") #plots the time vs. population density of both predator H and prey P with red being H and blue being P


#changing of initial parameters done individually with same steps as the initial model simulation


#lower birth rate
params1=c(0.25,0.02,0.1,0.2)
sim1=ode(y=y0, times=times, func = lotka_sim, parms = params1)
out1=data.frame(time=sim1[,1], H=sim1[,2], P=sim1[,3])
output1=melt(out1, id.vars="time", variable.name="Group", value.name="Density")
ggplot(output1, aes(x=time, y=Density))+geom_line(aes(color=Group))+theme_classic()+ylab("Population Density")+ggtitle("Lower Birth Rate")

#higher birth rate
params2=c(1.5,0.02,0.1,0.2)
sim2=ode(y=y0, times=times, func = lotka_sim, parms = params2)
out2=data.frame(time=sim2[,1], H=sim2[,2], P=sim2[,3])
output2=melt(out2, id.vars="time", variable.name="Group", value.name="Density")
ggplot(output2, aes(x=time, y=Density))+geom_line(aes(color=Group))+theme_classic()+ylab("Population Density")+ggtitle("Higher Birth Rate")

#lower predator attack rate
params3=c(0.5,0.01,0.1,0.2)
sim3=ode(y=y0, times=times, func = lotka_sim, parms = params3)
out3=data.frame(time=sim3[,1], H=sim3[,2], P=sim3[,3])
output3=melt(out3, id.vars="time", variable.name="Group", value.name="Density")
ggplot(output3, aes(x=time, y=Density))+geom_line(aes(color=Group))+theme_classic()+ylab("Population Density")+ggtitle("Lower Predator Attack Rate")

#higher predator attack rate
params4=c(0.5,0.06,0.1,0.2)
sim4=ode(y=y0, times=times, func = lotka_sim, parms = params4)
out4=data.frame(time=sim4[,1], H=sim4[,2], P=sim4[,3])
output4=melt(out4, id.vars="time", variable.name="Group", value.name="Density")
ggplot(output4, aes(x=time, y=Density))+geom_line(aes(color=Group))+theme_classic()+ylab("Population Density")+ggtitle("Higher Predator Attack Rate")

#lower conversion efficiency
params5=c(0.5,0.02,0.01,0.2)
sim5=ode(y=y0, times=times, func = lotka_sim, parms = params5)
out5=data.frame(time=sim5[,1], H=sim5[,2], P=sim5[,3])
output5=melt(out5, id.vars="time", variable.name="Group", value.name="Density")
ggplot(output5, aes(x=time, y=Density))+geom_line(aes(color=Group))+theme_classic()+ylab("Population Density")+ggtitle("Lower Conversion Efficiency")

#higher conversion efficiency
params6=c(0.5,0.02,0.3,0.2)
sim6=ode(y=y0, times=times, func = lotka_sim, parms = params6)
out6=data.frame(time=sim6[,1], H=sim6[,2], P=sim6[,3])
output6=melt(out6, id.vars="time", variable.name="Group", value.name="Density")
ggplot(output6, aes(x=time, y=Density))+geom_line(aes(color=Group))+theme_classic()+ylab("Population Density")+ggtitle("Higher Conversion Efficiency")

#lower predator death rate
params7=c(0.5,0.02,0.1,0.05)
sim7=ode(y=y0, times=times, func = lotka_sim, parms = params7)
out7=data.frame(time=sim7[,1], H=sim7[,2], P=sim7[,3])
output7=melt(out7, id.vars="time", variable.name="Group", value.name="Density")
ggplot(output7, aes(x=time, y=Density))+geom_line(aes(color=Group))+theme_classic()+ylab("Population Density")+ggtitle("Lower Predator Death Rate")

#higher predator death rate
params8=c(0.5,0.02,0.1,0.6)
sim8=ode(y=y0, times=times, func = lotka_sim, parms = params8)
out8=data.frame(time=sim8[,1], H=sim8[,2], P=sim8[,3])
output8=melt(out8, id.vars="time", variable.name="Group", value.name="Density")
ggplot(output8, aes(x=time, y=Density))+geom_line(aes(color=Group))+theme_classic()+ylab("Population Density")+ggtitle("Higher Predator Death Rate")



## Part 2: R-M Model simulations

## Part 3: Paradox of Enrichment
