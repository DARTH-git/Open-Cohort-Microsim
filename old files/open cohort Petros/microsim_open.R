
###############  simple 3-state microsimulation model              #############
#########             for the CE16 NIHES course            #####################

# Developed by:
# the Decision Analysis in R for Technologies in Health (DARTH) group
# Fernando Alarid-Escudero, PhD (1) 
# Eva A. Enns, MS, PhD (1)	
# M.G. Myriam Hunink, MD, PhD (2,3)
# Hawre J. Jalal, MD, PhD (4) 
# Eline M. Krijkamp, MSc (2)	
# Petros Pechlivanoglou, PhD (5) 

# In collaboration of: 		
# 1 University of Minnesota School of Public Health, Minneapolis, MN, USA
# 2 Erasmus MC, Rotterdam, The Netherlands
# 3 Harvard T.H. Chan School of Public Health, Boston, USA
# 4 University of Pittsburgh Graduate School of Public Health, Pittsburgh, PA, USA
# 5 The Hospital for Sick Children, Toronto and University of Toronto, Toronto ON, Canada

################################################################################
# Please cite our publications when using this code
# darthworkgroup.com 
## Jalal H, et al. An Overview of R in Health Decision Sciences. 
# Med. Decis. Making. 2017; 37(3): 735-746. 
## Krijkamp EM, et al. Microsimulation modeling for health decision sciences 
# using R: a tutorial. Med. Decis. Making. 2018; 38(3): 400-422.

################################################################################
# Copyright 2017, 
# THE HOSPITAL FOR SICK CHILDREN AND THE COLLABORATING INSTITUTIONS. 
# All rights reserved in Canada, the United States and worldwide.  
# Copyright, trademarks, trade names and any and all associated intellectual 
# property are exclusively owned by THE HOSPITAL FOR SICK CHILDREN and the 
# collaborating institutions and may not be used, reproduced, modified, 
# distributed or adapted in any way without written permission.

################################################################################

rm(list = ls())  # Delete everything that is in R's memory

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))  #set working directory as the folder where the course material is stored

#### 01 Load packages ####

#### 02 Load Functions ####

#### 03 Input Model Parameters ####
set.seed(1987)                            # set the seed 

# Model structure
v.n <- c("healthy", "pregnant", "healthy with baby", "dead")       # vector with state names
n.s <- length(v.n)                        # number of states
c.l <- 1/12                                 # cycle length
max.age <- 80                             # maximum age
n.t <- 80 / c.l                           # number of cycles
n.i <- 10000                              # number of individuals
d.r <- 1.5                                # annual discount rate
v.M_Init <- rep("healthy", times = n.i)   # initial state for all individual at the start of the model 

# load age distribution in the population
p.male <- 0.5 
v.ageDist <- read.csv("agedist.csv")
v.ages <- as.numeric(rownames(v.ageDist))

v.ageMale <- sample(x    =v.ages  , replace = T,size = n.i * p.male, 
                    prob = c(v.ageDist[1,2],diff(v.ageDist[,2])))


v.ageFemale <- sample(x  =v.ages  , replace = T,size = n.i * (1 - p.male), 
                    prob = c(v.ageDist[1,2],diff(v.ageDist[,2])))
# generate baseline characteristics

plot(ecdf(v.ageBase))
plot(v.ageDist[,2])

#### stop here ####

# Transition probabilities
p.HS <- 0.05         # probability healthy -> sick
p.HD <- 0.02         # probability healthy -> dead 
p.SD <- 0.1          # probability    sick -> dead

# Costs inputs
c.H <- 1500           # cost of one cycle in healthy state
c.S <- 5000           # cost of one cycle in sick state
c.D <- 0

# utility inputs
u.H     <- 1               # utility when healthy 
u.S     <- 0.85            # utility when sick 
u.D     <- 0               # utility when dead 

#### 03.1 Construct functions ####
#### 03.1.1 Probability functions ####
# Probs function outputs transition probabilities for next cycle
Probs <- function(M_it) { 
  # M_it:    current health state
  
  p.it <- vector("numeric", length(v.n))     # intiatilize vector of state transition probabilities
  
  # update p.it with the appropriate probabilities   
  p.it[M_it == "healthy"] <- c(1 - p.HD - p.HS, p.HS, p.HD)     # transition probabilities when healthy 
  p.it[M_it == "sick"]    <- c(0, 1 - p.SD, p.SD)               # transition probabilities when sick 
  p.it[M_it == "dead"]    <- c(0, 0, 1)                         # transition probabilities when dead      
  
  return(p.it)  # return transition probability
}  


#### 03.1.2 Cost functions ####
# Costs function calculates the cost accrued by an individual this cycle
Costs <- function (M_it) {
  # M_it: current health state
  c.it <- c()
  c.it[M_it == "dead"]    <- c.D     # costs at dead state
  c.it[M_it == "healthy"] <- c.H     # costs accrued by being healthy this cycle
  c.it[M_it == "sick"]    <- c.S     # costs accrued by being sick this cycle
  
  return(c.it)  # return costs accrued this cycle
}

#### 03.1.2 Effectiveness functions ####
# Effs function outputs QALYs accrued by an individual for this cycle
Effs <- function (M_it) {
  # M_it: current health state
  q.it <- c() 
  q.it[M_it == "dead"]    <- u.D        # QALYs at dead state
  q.it[M_it == "healthy"] <- u.H        # QALYs accrued by being healthy this cycle
  q.it[M_it == "sick"]    <- u.S        # QALYs accrued by being sick this cycle
  
  return(q.it)  # return the QALYs accrued this cycle
}


#### 03.1.3 Initiate the matrices ####
# m.M: health state for each patient at each cycle
# m.E: outcomes (e.g. QALYs) accrued by each patient at each cycle
# m.C: costs accrued by each individual at each cycle
m.M  <- m.E <- m.C <- 
  matrix(nrow = n.i, ncol = n.t + 1, 
         dimnames = list(paste("ind", 1:n.i, sep = " "),     # name the rows ind1, ind2, ind3, etc.
                         paste("cycle", 0:n.t, sep = " ")))  # name the columns cycle0, cycle1, cycle2, cycle3, etc.


#### 04 Model proces ####
p <- Sys.time() # start simulation 

# individuals 1 through n.i
for (i in 1:n.i) { # open individuals loop
  
  m.M[i, 1] <- v.M_Init[i]       # initial health state for individual i
  m.C[i, 1] <- Costs(m.M[i, 1])  # costs accrued individual i during cycle 0
  m.E[i, 1] <- Effs(m.M[i, 1])   # QALYs accrued individual i during cycle 0
  
  # cycles 1 thru n.t
  for (t in 1:n.t) { # open time loop
    
    # get transition probabilities based on previous health state
    v.p <- Probs(m.M[i, t])
    
    # sample the current health state based on transition probabilities v.p 
    m.M[i, t + 1] <- sample(x = v.n, prob = v.p, size = 1)
    
    # calculate costs and effects accrued by individual i for cycle t
    m.C[i, t + 1] <- Costs(m.M[i, t + 1])   # costs
    m.E[i, t + 1] <- Effs(m.M[i, t + 1])    # QALYs
    
  } # close time loop
  # Display simulation progress
  if(i/(n.i/10) == round(i/(n.i/10), 0)) { # display progress every 10%
    cat('\r', paste(i/n.i * 100, "% done", sep = " "))
  }
  
} # close individuals loop
comp.time = Sys.time() - p
comp.time

#### 05 Output ####
v.dw <- 1 / (1 + d.r) ^ (0:n.t) # calculate discount weight for each cycle based on discount rate d.r

tc <- m.C %*% v.dw    # total discounted cost per individual
te <- m.E %*% v.dw    # total discounted QALYs per individual 

tc_avg <- mean(tc)     # average discounted cost 
te_avg <- mean(te)     # average discounted QALYs

results <- data.frame("Total Cost" = tc_avg, "Total QALYs" = te_avg)
results



# VISUALIZE RESULTS
source("../functions/Microsim_visualize.R")

