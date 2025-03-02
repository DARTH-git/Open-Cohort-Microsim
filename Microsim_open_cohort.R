
###############  open corhot 3-state microsimulation model #############

# Developed by:
# the Decision Analysis in R for Technologies in Health (DARTH) group
# Fernando Alarid-Escudero, PhD (1) 
# Eva A. Enns, MS, PhD (1)	
# M.G. Myriam Hunink, MD, PhD (2,3)
# Hawre J. Jalal, MD, PhD (4) 
# Eline M. Krijkamp, MSc (2)	
# Petros Pechlivanoglou, PhD (5) 
# Alan Yang, MSc (5)

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

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))  # Set working directory as the folder where the course material is stored

#### 01 Install and Load packages ####
if (!require(plyr))       install.packages('plyr')       ; library(plyr)
if (!require(reshape2))   install.packages('reshape2')   ; library(reshape2)
if (!require(data.table)) install.packages('data.table') ; library(data.table)

#### 02 Load Functions ####
source('samplev.R') # this script contains samplev

#### 03 Input Model Parameters ####
set.seed(1987)                            # set the seed 

# Model structure
v.n <- c("H", "P", "D")                   # vector with state names
n.s <- length(v.n)                        # number of states
c.l <- 1/12                               # cycle length
max.fu <- 40                              # follow-up years 
n.t <- max.fu / c.l                       # number of cycles
n.i <- 10000                              # number of individuals
d.r <- 1.5                                # annual discount rate
v.M_Init <- rep("H", times = n.i)         # initial state for all individual at the start of the model 

# Start 9 cycles (months) before to allow pregnant women to enter the cohort
n.t <- n.t + 9

# Specify gender distribution 
p.female <- 0.51 

# Load age distribution in the population
v.ageDist   <- read.csv("agedist.csv")
# Only sample 14 to 44 year olds 
v.ageDist   <- v.ageDist[v.ageDist$Age >= 14 & v.ageDist$Age <= 44,]

# Generate baseline characteristics
v.ages      <- as.numeric(v.ageDist$Age)

v.ageFemale <- sample(x = v.ages, replace = T, size = n.i * p.female, 
                      prob = c(v.ageDist[1,2], diff(v.ageDist[,2])))

v.ageMale   <- sample(x = v.ages, replace = T, size = n.i * (1 - p.female), 
                      prob = c(v.ageDist[1,3], diff(v.ageDist[,3])))

df.X        = data.frame(ID = 1:n.i, Age = c(v.ageFemale, v.ageMale), 
                         Sex = c(rep('Female',length(v.ageFemale)), rep('Male',length(v.ageMale))), 
                         MotherID = rep(NA, n.i))
df.X$Sex    <- as.character(df.X$Sex)
df.X        <- df.X[sample(nrow(df.X)), ]
df.X$ID     <- rownames(df.X) <- 1:n.i

# load and format age- and sex- dependent mortatlities
p.mort      <- read.csv('mort.csv')
p.mort_m    <- subset(p.mort, select = c(Age, Male))
p.mort_m1   <- melt(p.mort_m, id='Age')
p.mort_m1   <- p.mort_m1[rep(seq_len(nrow(p.mort_m1)), each=5), ]
p.mort_m1   <- p.mort_m1[-c(1:4,6,107:110),]
p.mort_m1$Age <- 0:100
p.mort_f    <- subset(p.mort, select = c(Age, Female))
p.mort_f1   <- melt(p.mort_f, id = 'Age')
p.mort_f1   <- p.mort_f1[rep(seq_len(nrow(p.mort_f1)), each=5), ]
p.mort_f1   <- p.mort_f1[-c(1:4,6,107:110),]
p.mort_f1$Age <- 0:100
p.mort1     <- rbind(p.mort_m1, p.mort_f1)
colnames(p.mort1) <- colnames(p.mort_f1) <- colnames(p.mort_m1) <- c('Age', 'Sex', 'p.HD')

# Transition probabilities
p.HP <- 0.0014         # probability healthy -> pregnant
p.PD <- 0.01           # probability pregnant -> death

#### 04 Probability functions ####
Probs <- function(M_it) { 
  # Arguments
  # M_it: current health state
  # Returns
  # m.p.it: matrix of probabilities for all individuals in a given cycle 
  
  m.p.it           <- matrix(0, nrow = n.s, ncol = n.i)  # create matrix of state transition probabilities
  rownames(m.p.it) <-  v.n                               # give the state names to the rows
  
  # Determine eligibility of getting into pregnant state
  sex_female       <- df.X$Sex == 'Female' 
  age_14_44        <- df.X$Age >= 14 & df.X$Age <= 44
  female_14_44     <- sex_female & age_14_44
  
  sex_male         <- df.X$Sex == 'Male'
  age_not_14_44    <- df.X$Age < 14 | df.X$Age > 44
  cant_preg        <- sex_male | age_not_14_44
  
  # Look up baseline probability and rate of dying based on individual characteristics
  p.HD_f           <- p.mort_f1[df.X$Age[df.X$Sex=="Female"] + 1,]
  p.HD_m           <- p.mort_m1[df.X$Age[df.X$Sex=="Male"] + 1,]
  p.HD_all         <- rbind(p.HD_f, p.HD_m)
  p.HD_P           <- p.HD_all[M_it == "H" & female_14_44, "p.HD"]
  p.HD_noP         <- p.HD_all[M_it == "H" & cant_preg, "p.HD"]
  
  # Transition probabilities when healthy
  m.p.it[, M_it == "H" & female_14_44]  <- rbind(1 - p.HP - p.HD_P, p.HP, p.HD_P)   
  m.p.it[, M_it == "H" & cant_preg]     <- rbind(1 - p.HP - p.HD_noP, p.HP, p.HD_noP)
  
  # Transition probabilities when pregnant
  m.p.it[, M_it == "P" & df.X$tp <= 9]  <- c(0, 1 - p.PD, p.PD)           
  m.p.it[, M_it == "P" & df.X$tp > 9]   <- c(1 - p.PD, 0, p.PD)
  
  # Transition probabilities when dead  
  m.p.it[, M_it == "D"]  <- c(0, 0, 1)                                                           
  
  return(t(m.p.it)) # return transition probability
}  

# Record computation time of the model process
start_time <- Sys.time()

#### 05 Model process ####
# Initiate the matrices 
# m.M: health state for each patient at each cycle
m.M  <-  matrix(nrow = n.i, ncol = n.t + 1, 
                dimnames = list(paste("ind", 1:n.i, sep = " "),     # name the rows ind1, ind2, ind3, etc.
                                paste("cycle", 0:n.t, sep = " ")))  # name the columns cycle0, cycle1, cycle2, cycle3, etc.

m.M[, 1] <- v.M_Init            # initial health state for individual i

df.X$tp <- rep(0, nrow(df.X))   # determine eligibilty of giving birth (track time in pregancy)
n.i_hist <- c(n.i)              # store total number of individuals at each cycle
n.i_babies <- c(0)              # store number of newborn babies at each cycle

# Loop through the time cycles
for (t in 1:n.t) {
  
    # Update time in pregnancy after each cycle
    df.X$tp[m.M[, t] == 'P'] <- df.X$tp[m.M[, t] == 'P'] + 1
    
    # Get transition probabilities based on health state at t
    v.p <- Probs(m.M[,t])
    
    # Sample the next health state based on transition probabilities v.p 
    m.M[, t + 1]  <- samplev(v.p, 1)

    # Add newborn babies to baseline charasteristics dataframe
    baby.momIDs <- df.X$ID[m.M[, t] == "P" & df.X$tp > 9]
    if (length(baby.momIDs) > 0) {
    baby.IDs   <- c(nrow(df.X)+1):c(nrow(df.X)+length(baby.momIDs))
    baby.Ages  <- rep(0, length(baby.momIDs))
    baby.Sex   <- sample(c(1,2), prob=c(0.5,0.5), replace=T, size=length(baby.momIDs))
    baby.tp    <- rep(0, length(baby.momIDs))
    babies     <- as.data.frame(cbind(baby.IDs, baby.Ages, baby.Sex, baby.momIDs,baby.tp))
    babies$baby.Sex[babies$baby.Sex==1] <- 'Female'
    babies$baby.Sex[babies$baby.Sex==2] <- 'Male'
    colnames(babies) <- colnames(df.X)
    df.X <- rbindlist(list(df.X, babies))
    
    # Add newborn babies to the health state matrix
    n.babies   <- length(baby.momIDs)   
    n.i_babies <- c(n.i_babies, n.babies) 
    new_babies_m.M <- matrix(NA, nrow=n.babies, ncol=n.t+1)
    new_babies_m.M[ , t+1] <- rep('H', n.babies)  
    rownames(new_babies_m.M) <- paste('ind ', baby.IDs, sep='')
    m.M <- rbind(m.M, new_babies_m.M)
    }
    
    n.i <- nrow(df.X)              # update the number of individuals 
    n.i_hist <- c(n.i_hist, n.i)
    df.X$tp[df.X$tp > 9] <- 0      # reset time in pregnancy to zero after giving birth
    
    if (t %% c(1/c.l) == 0) {df.X$Age <- df.X$Age + 1}  # increase age by 1 every 12 months

    # Display simulation progress
    if(abs(t/(n.t/10) - round(t/(n.t/10), 0)) < 0.02) { # display progress every 10%
       cat('\r', paste(round(t/n.t * 100, 0) , "% done", sep = " "))
    }
}

# Show computation time of the model process
end_time <- Sys.time()
run_time <- end_time - start_time
run_time

# 3.296799 secs at n.i = 1000
# 39.14869 secs at n.i = 10,000

#### 06 Trace ####
# Plot the distribution of the population across health states over time (Trace)

# Remove cycles 0 - 9 (first 9 months)
m.M      <- m.M[, -c(1:10)]
n.i_hist <- n.i_hist[-c(1:10)]

# Calculate the proportion of individuals in each health state at each cycle
m.TR           <- t(apply(m.M, 2, function(x) table(factor(x, levels = v.n, ordered = TRUE)))) 
m.TR           <- m.TR / n.i_hist                     
colnames(m.TR) <- v.n                                    
rownames(m.TR) <- paste("Cycle", 10:n.t, sep = " ")       

# Plot health state trace
matplot(10:n.t, m.TR, type = 'l', 
        ylab = "Proportion of cohort",
        xlab = "Cycle",
        main = "Health state trace")
legend("topright", v.n, col = 1:n.s, lty = 1:n.s, bty = "n")  

# Count the number of individuals in each health state at each cycle
m.TR_count           <- t(apply(m.M, 2, function(x) table(factor(x, levels = v.n, ordered = TRUE)))) 
colnames(m.TR_count) <- v.n                                    
rownames(m.TR_count) <- paste("Cycle", 10:n.t, sep = " ")      

# Plot health state trace
matplot(10:n.t, m.TR_count, type = 'l', 
        ylab = "Number of individuals",
        xlab = "Cycle",
        main = "Health state trace - total counts")
legend("topright", v.n, col = 1:n.s, lty = 1:n.s, bty = "n")  


