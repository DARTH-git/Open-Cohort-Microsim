
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
v.n <- c("healthy", "healthy with baby", "sick", "dead")       # vector with state names
n.s <- length(v.n)                        # number of states
n.t <- 60                                 # number of cycles
n.i <- 10000                              # number of individuals
d.r <- 3                                  # annual discount rate
v.M_Init <- rep("healthy", times = n.i)   # initial state for all individual at the start of the model 

# Transition probabilities
p.HS <- 0.05         # probability healthy -> sick
p.HHB <- 0.005        # probability healthy -> healthy with baby
p.HBS <- 0.05        # probability healthy with baby -> sick
p.HBD <- 0.02         # probability healthy with baby -> dead
p.HD <- 0.02         # probability healthy -> dead 
p.SD <- 0.1          # probability    sick -> dead

# Costs inputs
c.H <- 1500           # cost of one cycle in healthy state
c.HB <- 8000          # cost of one cycle in healthy with baby state
c.S <- 5000           # cost of one cycle in sick state
c.D <- 0

# utility inputs
u.H     <- 1               # utility when healthy 
u.HB    <- 0.98            # utility when healthy with baby
u.S     <- 0.85            # utility when sick 
u.D     <- 0               # utility when dead 

#### 03.1 Construct functions ####
#### 03.1.1 Probability functions ####
# Probs function outputs transition probabilities for next cycle
Probs <- function(M_it) { 
  # M_it:    current health state
  
  p.it <- vector("numeric", length(v.n))     # intiatilize vector of state transition probabilities
  
  # update p.it with the appropriate probabilities   
  p.it[M_it == "healthy"] <- c(1 - p.HD - p.HS - p.HHB, p.HHB, p.HS, p.HD)     # transition probabilities when healthy
  p.it[M_it == "healthy with baby"] <- c(1 - p.HBS - p.HBD, 0, p.HBS, p.HBD)   # transition probabilities when healthy with baby 
  p.it[M_it == "sick"]    <- c(0, 0, 1 - p.SD, p.SD)               # transition probabilities when sick 
  p.it[M_it == "dead"]    <- c(0, 0, 0, 1)                         # transition probabilities when dead      
  
  return(p.it)  # return transition probability
}  


#### 03.1.2 Cost functions ####
# Costs function calculates the cost accrued by an individual this cycle
Costs <- function (M_it) {
  # M_it: current health state
  c.it <- c()
  c.it[M_it == "dead"]    <- c.D     # costs at dead state
  c.it[M_it == "healthy"] <- c.H     # costs accrued by being healthy this cycle
  c.it[M_it == "healthy with baby"] <- c.HB     # costs accrued by being healthy with baby this cycle
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
  q.it[M_it == "healthy with baby"] <- u.HB      # # QALYs accrued by being healthy with baby this cycle
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


#### 04 Model process ####
p <- Sys.time() # start simulation 

# initial cohort size
n.i_init <- n.i
# individuals 1 through n.i
i <- 1
while (i <= n.i) { # open individuals loop
  
  if (i > n.i_init) {
    # since babies do not enter at cycle 0, we need to start their time loop at the cycle they were born in
    tt <- min(which(!is.na(m.M[i, ])))  
  } else {
    tt <- 1 # normal individuals start at cycle 0 
    m.M[i, 1] <- v.M_Init[i]       # initial health state for individual i
    m.C[i, 1] <- Costs(m.M[i, 1])  # costs accrued individual i during cycle 0
    m.E[i, 1] <- Effs(m.M[i, 1])   # QALYs accrued individual i during cycle 0    
  }

  # cycles 1 thru n.t
  for (t in tt:n.t) { # open time loop
    
    # get transition probabilities based on previous health state
    v.p <- Probs(m.M[i, t])
    
    # sample the current health state based on transition probabilities v.p 
    m.M[i, t + 1] <- sample(x = v.n, prob = v.p, size = 1)
    
    # calculate costs and effects accrued by individual i for cycle t
    m.C[i, t + 1] <- Costs(m.M[i, t + 1])   # costs
    m.E[i, t + 1] <- Effs(m.M[i, t + 1])    # QALYs
    
    # if individual i enters the healthy with baby state, the baby enters the model
    if (m.M[i, t + 1] == "healthy with baby") {
      m.M <- rbind(m.M, rep(NA, ncol(m.M)))  # add the baby to the m.M matrix
      m.M[nrow(m.M), t + 1] <- 'healthy'     # baby begins in the healthy state
      m.C <- rbind(m.C, rep(0, ncol(m.C)))  # add the baby to the m.C matrix, previous costs assigned to 0
      m.C[nrow(m.C), t + 1] <- Costs(m.M[nrow(m.C), t + 1])    # cost in the healthy state
      m.E <- rbind(m.E, rep(0, ncol(m.E)))  # add the baby to the m.M matrix, previous utilities assigned to 0
      m.E[nrow(m.E), t + 1] <- Effs(m.M[nrow(m.E), t + 1])     # utility in the healthy state
      # add rownames
      rownames(m.M)[nrow(m.M)] <- rownames(m.C)[nrow(m.C)] <- rownames(m.E)[nrow(m.E)] <- 
        paste("ind ", as.character(nrow(m.M)), ",", " baby", sep = "")
      # increase the n.i
      n.i <- n.i + 1
    }
    
  } # close time loop
  # Display simulation progress
  if(i/(n.i/10) == round(i/(n.i/10), 0)) { # display progress every 10%
    cat('\r', paste(i/n.i * 100, "% done", sep = " "))
  }
  
  i <- i + 1
  
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

# examine the babies
tail(m.M, 3)
tail(m.C, 3)
tail(m.E, 3)


#### 06 Plot ####
#### 06.1 Histogram ####
# Histogram showing variability in individual total costs
plot(density(tc), main = paste("Total cost per person"), xlab = "Cost ($)")

# Histogram showing variability in individual total QALYs
plot(density(te), main = paste("Total QALYs per person"), xlab = "QALYs")

#### 06.2 Trace ####
# Plot the distribution of the population across health states over time (Trace)
# count the number of individuals in each health state at each cycle
m.TR <- t(apply(m.M, 2, function(x) table(factor(x, levels = v.n, ordered = TRUE)))) 
m.TR <- m.TR / n.i                                       # calculate the proportion of individuals 
colnames(m.TR) <- v.n                                    # name the rows of the matrix
rownames(m.TR) <- paste("Cycle", 0:n.t, sep = " ")       # name the columns of the matrix

# Plot trace 
matplot(0:n.t, m.TR, type = 'l', 
        ylab = "Proportion of cohort",
        xlab = "Cycle",
        main = "Health state trace")                 
legend("topright", v.n, col = 1:n.s, lty = 1:n.s, bty = "n")  # add a legend to the graph

