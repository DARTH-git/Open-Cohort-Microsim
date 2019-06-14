
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

# healthy, pregnant, death
# after 9 cycles you move back to healthy
# no need for utilities or costs
# add mom id beside the baby id
# sample age at baseline (15-40), look at CE16 examples on age and sex dependent probabilities
# time horizon: 2 years (24 months)
# at baseline, mother ids are NA (mother id means what is the id of your mother)
# dataset: ID, AGE, MOTHER.ID
# new baby: age 0, his/her mother id is the id of his/her mom
# new baby sex: 50% male
# track time in pregnancy, if time in pregancy == 9, 2 things will happen: 
#    1. either go back to healthy or die (dont stay)
#    2. n.i <- n.i + 1
# work with for loop version then samplev

rm(list = ls())  # Delete everything that is in R's memory

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))  #set working directory as the folder where the course material is stored

#### 01 Load packages ####
library(plyr)
library(reshape2)

#### 02 Load Functions ####

#### 03 Input Model Parameters ####
set.seed(1987)                            # set the seed 

# Model structure
v.n <- c("H", "P", "D")       # vector with state names
n.s <- length(v.n)                        # number of states
c.l <- 1/12                                 # cycle length
max.age <- 40                            # maximum age
n.t <- max.age / c.l                           # number of cycles
n.i <- 100                              # number of individuals
d.r <- 1.5                                # annual discount rate
v.M_Init <- rep("H", times = n.i)   # initial state for all individual at the start of the model 

# load age distribution in the population
p.male <- 0.5 
v.ageDist <- read.csv("agedist.csv")
# only allow 15 to 40 year olds 
v.ageDist <- v.ageDist[v.ageDist$Age >= 15 & v.ageDist$Age <= 40,]

# generate baseline characteristics
v.ages <- as.numeric(v.ageDist$Age)

v.ageFemale <- sample(x    = v.ages, replace = T,size = n.i * (1 - p.male), 
                      prob = c(v.ageDist[1,2],diff(v.ageDist[,2])))

v.ageMale <- sample(x    = v.ages, replace = T, size = n.i * p.male, 
                    prob = c(v.ageDist[1,3], diff(v.ageDist[,3])))

df.X = data.frame(ID = 1:n.i, Age = c(v.ageFemale, v.ageMale), 
                  Sex = c(rep('Female',length(v.ageFemale)), rep('Male',length(v.ageMale))), 
                  MotherID = rep(NA, n.i))
df.X$Sex <- as.character(df.X$Sex)
df.X <- df.X[sample(nrow(df.X)),]
df.X$ID <- rownames(df.X) <- 1:n.i

# load mortatlities
p.mort <- read.csv('mort.csv')
p.mort_m <- subset(p.mort, select = c(Age, Male))
p.mort_m1 <- melt(p.mort_m, id='Age')
p.mort_m1 <- p.mort_m1[rep(seq_len(nrow(p.mort_m1)), each=5),]
p.mort_m1 <- p.mort_m1[-c(1:4,6,107:110),]
p.mort_m1$Age <- 0:100
p.mort_f <- subset(p.mort, select = c(Age, Female))
p.mort_f1 <- melt(p.mort_f, id='Age')
p.mort_f1 <- p.mort_f1[rep(seq_len(nrow(p.mort_f1)), each=5),]
p.mort_f1 <- p.mort_f1[-c(1:4,6,107:110),]
p.mort_f1$Age <- 0:100
p.mort1 <- rbind(p.mort_m1, p.mort_f1)
colnames(p.mort1) <- c('Age', 'Sex', 'p.HD')

# Transition probabilities
p.HP <- 0.005         # probability healthy -> pregnant
p.PD <- 0.01         # probability pregnant -> death

# Cost inputs
c.H <- 1500           # cost of one cycle in healthy State
c.P <- 5000           # cost of one cycle in pregnant State
c.D <- 0              # cost of one cycle in death State

# utility inputs
u.H     <- 1          # utility of one cycle in healthy State
u.P     <- 0.85       # utility in pregnant State
u.D     <- 0          # utility in death State

#### 04.2 Dynamic characteristics 
# Specify the initial health state of the individuals 
# everyone begins in the healthy state (in this example)

#####

Probs <- function(M_it, tp) { 
  # M_it:    current health state
  
  p.it <- vector("numeric", length(v.n))     # intiatilize vector of state transition probabilities
  
  p.HD_all <- p.mort1[df.X$Age+1,]
  p.HD <- p.HD_all$p.HD[i]
  
  if (p.HD_all$Sex[i] == 'Male' | p.HD_all$Age[i] < 15 ) {p.HP <- 0}
  
  # update p.it with the appropriate probabilities   
  p.it[M_it == "H"] <- c(1 - p.HP - p.HD, p.HP, p.HD)      # transition probabilities when healthy 
  
  if (tp <= 9) {                                           # transition probabilities when pregnant
    p.it[M_it == "P"] <- c(0, 1-p.PD, p.PD) 
  } else {
    p.it[M_it == "P"] <- c(1 - p.PD, 0, p.PD)          
  }
 
   p.it[M_it == "D"]    <- c(0, 0, 1)                         # transition probabilities when dead      
  
  return(p.it)  # return transition probability
}  


#### 03.1.3 Initiate the matrices ####
# m.M: health state for each patient at each cycle
m.M  <-  matrix(nrow = n.i, ncol = n.t + 1, 
                dimnames = list(paste("ind", 1:n.i, sep = " "),     # name the rows ind1, ind2, ind3, etc.
                                paste("cycle", 0:n.t, sep = " ")))  # name the columns cycle0, cycle1, cycle2, cycle3, etc.

m.M[, 1] <- "H"       # initial health state for individual i


#### 04 Model process ####

start_time <- Sys.time()

df.X$tp <- rep(0, nrow(df.X))

for (t in 1:n.t) {
  
  for (i in 1:n.i) {
    
    if (m.M[i, t] == "P") {
      df.X$tp[i] <- df.X$tp[i] + 1
    }
    
    # get transition probabilities based on previous health state
    tp <- df.X$tp[i]
    v.p <- Probs(m.M[i, t], tp)
    
    # sample the current health state based on transition probabilities v.p 
    m.M[i, t + 1] <- sample(x = v.n, prob = v.p, size = 1)
    
    # new baby enters cohort when after 9 months of pregnancy
    if (tp > 9) {
    m.M <- rbind(m.M, rep(NA, ncol(m.M)))  # add the baby to the m.M matrix
    m.M[nrow(m.M), t + 1] <- 'H'     # baby begins in the healthy state
    # add rownames
    rownames(m.M)[nrow(m.M)] <- paste("ind ", as.character(nrow(m.M)), ",", " baby", sep = "")
    # add baby to the cohort matrix
    momID <- df.X$ID[i]  # store mother ID
    df.X <- rbind(df.X, c(nrow(df.X)+1, 0, 0, momID, 0))
    df.X$Sex[nrow(df.X)] <- sample(c('Female','Male'), prob=c(0.5,0.5), size=1)
    df.X$tp[i] <- 0 # time in pregnant state becomes 0 for the mother
    # increase n.i for the for loop
    n.i <- n.i + 1
    }
  }
  
  if (t %% c(1/c.l) == 0) {df.X$Age <- df.X$Age + 1}  # increase age by 1 every 12 months
  
  # Display simulation progress
  if(t/(n.t/10) == round(t/(n.t/10), 0)) { # display progress every 10%
    cat('\r', paste(t/n.t * 100, "% done", sep = " "))
  }
}

end_time <- Sys.time()

end_time - start_time

# 6.495949 secs at n.i = 100
# 4.099744 mins at n.i = 1000

