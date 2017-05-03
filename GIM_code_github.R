## CODE FOR FITTING THE GIM MODEL AND FOR SIMULATING DATA FROM THE GIM MODEL
# first load "GIM_r_functions.R";
# needs package "stats";
# needs vectors of pairwise differences (x1,x2,x3) 
# and of relative mutation rates (r1,r2,r3)

#example of simulation
# (Mc denotes contemporary migration rate)
data<-sim.GIM(a=2,theta=2,b=1.5,c1=2,c2=2,tau1=0,tau0=1.5,M1=0.2,M2=0.2,
        M1c=0,M2c=0,N=c(2000,2000,4000))


r1<-rep(1,data$N1)
r2<-rep(1,data$N2)
r3<-rep(1,data$N3)
x1<-data$x1
x2<-data$x2
x3<-data$x3



# to fit a model just run, for example,
gim.list<-gim("IIM.3")

# or, if you need to input a vector of initial values,
gim.list<-gim("IIM.3",iv=rep(2.5,8))

# to calculate, for example, the 95% profile likelihood confidence interval for 
# the parameter "M2" under the IIM.3 model run:

CI.pl(gim.list,"M2","IIM.3",quantile=7.677432, ub=20) 

# CI.pl searches for the profile likelihood CI in the region (0,ub).
# Currentl, it can be applied to models "iso.1", "IM.4", "IIM.3"and "IIM.5".


#### MODELS AVAILABLE TO FIT:

### GIM MODELS (M_prime denotes contemporary migration rate; there are 
### no constraints on population sizes) 

# GIM.1: full GIM model (11 parameters)

# GIM.2: secondary contact model (ancestral migration rates fixed at zero) (9 pars)

# GIM.3: GIM model with symmetric migration rates in both periods (9 pars)

# GIM.4: secondary contact model with symmetric migration
# (ancestral migration rates fixed at zero) (8 pars)


### MODELS IN FIGURE 7 of Costa & Wilkinson-Herbots (Genetics,2017):

# iso.1 (4 parameters)
# IM.1 (6 pars)
# IIM.1 (7 pars)
# IIM.2 (9 pars)
# IIM.3 (8 pars)


### OTHER MODELS

# IM.2: "IM1" from figure 4, IIM paper, model with M1=0 (5 parameters)

# IM.3: "IM1" model with M2=0 (5 pars)

# IM.3: "IM1" model with M2=0 (5 pars)

# IM.4: "IM1" model with M1=M2 (5 pars)

# IIM.4: "IIM2" model with M2=0 (8 pars)

# IIM.5: "IIM2" model with M1=M2 and equal pop sizes during migration (7 pars) 

# IIM.6: "IIM2" model with M1=M2 and allowing for unequal population sizes (8 pars)

# iso.2: "IIM2" model with M1=M2=0 (7 pars)





