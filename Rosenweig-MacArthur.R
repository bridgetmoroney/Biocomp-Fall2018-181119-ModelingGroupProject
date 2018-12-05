install.packages("deSolve")
library("deSolve") #calling desolve library
library(ggplot2) #calling ggplot2 library
library(reshape2) #calling reshaoe2 library

#Simulating Predator Prey Dynamics: Part II

#Modeling the Rosenweig MacArthur prey-predator dynamics:

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
