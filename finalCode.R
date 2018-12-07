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





#Simulating Predator Prey Dynamics: Part II

#Modeling the Rosenweig MacArthur prey-predator dynamics:

install.packages("deSolve")
library("deSolve") #calling desolve library
library(ggplot2) #calling ggplot2 library
library(reshape2) #calling reshape2 library

#function will have six parameters and two state variables
ddSimRM=function(t,y,p){ #defining function variables where y are state variables and p are parameters
  H=y[1] #H is herbivore state variable
  P=y[2] #P is predator state variable
  b=p[1] #Prey birth rate parameter
  a=p[2] #Competition coefficient parameter
  w=p[3] #Predator attack rate parameter
  d=p[4] #Saturation term parameter
  e=p[5] #Conversion efficiency of prey to predators parameter
  s=p[6] #Predator death rate parameter
  
  modelH=b*H*(1-a*H)-w*(H/(d+H))*P #defining two prey-predator functions
  modelP=e*w*(H/(d+H))*P-(s*P)
  
  return(list(c(modelH,modelP))) #setting function to return two models above
}

params1=c(0.8,0.001,5,400,0.07,0.2) #case 1 with b=0.8, e=0.07, s=0.2, w=5, d=400, a=0.001
initial1=c(500,120) #initial conditions Ho=500 and Po=120
times1=seq(1,100, by=0.1) #running simulation from time 1 to 100 over 0.1 time step

modelSim1=ode(y=initial1,times=times1,func=ddSimRM,parms=params1) #using ordinary differential equation to simulate functions above
modelOutput1=data.frame(time=modelSim1[,1],H=modelSim1[,2],P=modelSim1[,3]) #creating dataframe to store simulation values
modelOutput1=melt(modelOutput1,id.vars = "time",variable.name = "Group",value.name = "Density") #melting density data in order to create two line plot
ggplot(modelOutput1,aes(x=time,y=Density))+geom_line(aes(color=Group))+theme_classic()+xlab("Time")+ggtitle("Standard Parameters")+theme(plot.title = element_text(hjust = 0.5)) #creating ggplot with simulation data

#Will now consider 6 other cases where a single parameter is both increased and decreased
#The same code was used to increase and decrease each parameter, and a new parameter value was
#put in when wanting to change the parameter and all other values in global environment redefined

#Case 2 changing b parameter
#Increased b to 2.4 and decreased b to 0.26 (multiple of 3)
params2=c(0.26,0.001,5,400,0.07,0.2)
initial2=c(500,120)
times2=seq(1,100, by=0.1)

modelSim2=ode(y=initial2,times=times2,func=ddSimRM,parms=params2) #comments on new simulation model is same as standard parameters
modelOutput2=data.frame(time=modelSim2[,1],H=modelSim2[,2],P=modelSim2[,3])
modelOutput2=melt(modelOutput2,id.vars = "time",variable.name = "Group",value.name = "Density")
ggplot(modelOutput2,aes(x=time,y=Density))+geom_line(aes(color=Group))+theme_classic()+xlab("Time")+ggtitle("Decreased b Parameter")+theme(plot.title = element_text(hjust = 0.5))

#Case 3 changing e parameter
#Increased e to 0.14 and decreased e to 0.035 (multiple of 2)
params3=c(0.8,0.001,5,400,0.14,0.2)
initial3=c(500,120)
times3=seq(1,100, by=0.1)

modelSim3=ode(y=initial3,times=times3,func=ddSimRM,parms=params3) #comments on new simulation model is same as standard parameters
modelOutput3=data.frame(time=modelSim3[,1],H=modelSim3[,2],P=modelSim3[,3])
modelOutput3=melt(modelOutput3,id.vars = "time",variable.name = "Group",value.name = "Density")
ggplot(modelOutput3,aes(x=time,y=Density))+geom_line(aes(color=Group))+theme_classic()+xlab("Time")+ggtitle("Increased e Parameter")+theme(plot.title = element_text(hjust = 0.5))


#Case 4 changing s parameter
#Increased s to 0.4 and decreased s to 0.1 (multiple of 2)
params4=c(0.8,0.001,5,400,0.07,0.1)
initial4=c(500,120)
times4=seq(1,100, by=0.1)

modelSim4=ode(y=initial4,times=times4,func=ddSimRM,parms=params4) #comments on new simulation model is same as standard parameters
modelOutput4=data.frame(time=modelSim4[,1],H=modelSim4[,2],P=modelSim4[,3])
modelOutput4=melt(modelOutput4,id.vars = "time",variable.name = "Group",value.name = "Density")
ggplot(modelOutput4,aes(x=time,y=Density))+geom_line(aes(color=Group))+theme_classic()+xlab("Time")+ggtitle("Decreased S Parameter")+theme(plot.title = element_text(hjust = 0.5))

#Case 5 changing w
#Increased w to 15 and decreased w to 1.67 (multiple of 3)
params5=c(0.8,0.001,15,400,0.07,0.2)
initial5=c(500,120)
times5=seq(1,100, by=0.1)

modelSim5=ode(y=initial5,times=times5,func=ddSimRM,parms=params5) #comments on new simulation model is same as standard parameters
modelOutput5=data.frame(time=modelSim5[,1],H=modelSim5[,2],P=modelSim5[,3])
modelOutput5=melt(modelOutput5,id.vars = "time",variable.name = "Group",value.name = "Density")
ggplot(modelOutput5,aes(x=time,y=Density))+geom_line(aes(color=Group))+theme_classic()+xlab("Time")+ggtitle("Increased w Parameter")+theme(plot.title = element_text(hjust = 0.5))

#Case 6 changing d
#Increased d to 800 and decreased d to 200 (multiple of 2)
params6=c(0.8,0.001,5,200,0.07,0.2)
initial6=c(500,120)
times6=seq(1,100, by=0.1)

modelSim6=ode(y=initial6,times=times6,func=ddSimRM,parms=params6) #comments on new simulation model is same as standard parameters
modelOutput6=data.frame(time=modelSim6[,1],H=modelSim6[,2],P=modelSim6[,3])
modelOutput6=melt(modelOutput6,id.vars = "time",variable.name = "Group",value.name = "Density")
ggplot(modelOutput6,aes(x=time,y=Density))+geom_line(aes(color=Group))+theme_classic()+xlab("Time")+ggtitle("Decreased d Parameter")+theme(plot.title = element_text(hjust = 0.5))

#Case 7 changing a
#Increased a to 0.003 and decreased a to 0.00033 (multiple of 3)
params7=c(0.8,0.003,5,400,0.07,0.2)
initial7=c(500,120)
times7=seq(1,100, by=0.1)

modelSim7=ode(y=initial7,times=times7,func=ddSimRM,parms=params7) #comments on new simulation model is same as standard parameters
modelOutput7=data.frame(time=modelSim7[,1],H=modelSim7[,2],P=modelSim7[,3])
modelOutput7=melt(modelOutput7,id.vars = "time",variable.name = "Group",value.name = "Density")
ggplot(modelOutput7,aes(x=time,y=Density))+geom_line(aes(color=Group))+theme_classic()+xlab("Time")+ggtitle("Increased alpha Parameter")+theme(plot.title = element_text(hjust = 0.5))

#Simulating Predator Prey Dynamics: Part III

#Testing carrying capacity and Paradox of enrichment:
#Carrying capacity will be varied by changing alpha to affect carrying capacity in 200 increments
#Each alpha value should be enetered into the second parameter value and run to obtain 7 plots

#alphas: 0.00125, 0.001, 0.000833, 0.00071, 0.000625, 0.00056, and 0.0005
paramscc=c(0.8,0.0005,5,400,0.07,0.2)
initialcc=c(500,120)
timescc=seq(1,100, by=0.1)

modelSimcc=ode(y=initialcc,times=timescc,func=ddSimRM,parms=paramscc)

modelOutputcc=data.frame(time=modelSimcc[,1],H=modelSimcc[,2],P=modelSimcc[,3]) #comments on new simulation model is same as standard parameters
modelOutputcc=melt(modelOutputcc,id.vars = "time",variable.name = "Group",value.name = "Density")
ggplot(modelOutputcc,aes(x=time,y=Density))+geom_line(aes(color=Group))+theme_classic()+xlab("Time")+ggtitle("K=2,000")+theme(plot.title = element_text(hjust = 0.5))

