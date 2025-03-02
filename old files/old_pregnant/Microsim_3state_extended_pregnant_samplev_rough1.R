rm(list = ls())  # Delete everything that is in R's memory

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))  #set working directory as the folder where the course material is stored

#### 01 Load packages ####
library(plyr)
library(reshape2)
library(data.table)

#### 02 Load Functions ####
source('samplev.R')

#### 03 Input Model Parameters ####
set.seed(1987)                            # set the seed 

# Model structure
v.n <- c("H", "P", "D")       # vector with state names
n.s <- length(v.n)                        # number of states
c.l <- 1/12                                 # cycle length
max.age <- 40                            # maximum age
# n.t <- max.age / c.l                           # number of cycles
n.t <- max.age/c.l
n.i <- 1000                             # number of individuals
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

# load and format mortatlities
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
p.HP <- 0.005        # probability healthy -> pregnant
p.PD <- 0.01         # probability pregnant -> death

# Cost inputs
c.H <- 1500           # cost of one cycle in healthy State
c.P <- 5000           # cost of one cycle in pregnant State
c.D <- 0              # cost of one cycle in death State

# utility inputs
u.H     <- 1          # utility of one cycle in healthy State
u.P     <- 0.85       # utility in pregnant State
u.D     <- 0          # utility in death State


##########################################################################################

#### 03.1.3 Initiate the matrices ####
# m.M: health state for each patient at each cycle
m.M  <-  matrix(nrow = n.i, ncol = n.t + 1, 
                dimnames = list(paste("ind", 1:n.i, sep = " "),     # name the rows ind1, ind2, ind3, etc.
                                paste("cycle", 0:n.t, sep = " ")))  # name the columns cycle0, cycle1, cycle2, cycle3, etc.

m.M[, 1] <- "H"       # initial health state for individual i

# determine eligibilty of giving birth (track time in pregancy)
df.X$tp <- rep(0, nrow(df.X))

# store total number of individuals at each cycle
n.i_hist <- c(n.i)

#### 03.1.1 Probability functions ####
Probs <- function(M_it) { 
  # Arguments
  # M_it: current health state
  # Returns
  # m.p.it: matrix of probabilities for that individual at that cycle 
  
  m.p.it          <- matrix(0, nrow = n.s, ncol = n.i)  # create matrix of state transition probabilities
  rownames(m.p.it) <-  v.n                              # give the state names to the rows
  
  # determine eligibility of getting into pregnant state
  sex_female <- df.X$Sex == 'Female' 
  age_15_40 <- df.X$Age >= 15 & df.X$Age <= 40
  female_15_40 <- sex_female & age_15_40
  
  sex_male <- df.X$Sex == 'Male'
  age_not_15_40 <- df.X$Age < 15 | df.X$Age > 40
  cant_preg <- sex_male | age_not_15_40
  
  # lookup baseline probability and rate of dying based on individual characteristics
  p.HD_all <- p.mort1[df.X$Age+1,]
  p.HD_P     <- p.HD_all[m.M[, t] == "H" & female_15_40, "p.HD"]
  p.HD_noP  <- p.HD_all[m.M[, t] == "H" & cant_preg, "p.HD"]
  
  # update the v.p with the appropriate probabilitie, v.p will be fed to samplev 
  # (v.p is a n.i by n.s matrix, rows sum to 1)  
  
  # transition probabilities when healthy
  m.p.it[, m.M[, t] == "H" & female_15_40]  <- rbind(1 - p.HP - p.HD_P, p.HP, p.HD_P)   
  m.p.it[, m.M[, t] == "H" & cant_preg]  <- rbind(1 - p.HP - p.HD_noP, p.HP, p.HD_noP)
  
  # transition probabilities when pregnant
  m.p.it[, m.M[, t] == "P" & df.X$tp <= 9] <- c(0, 1 - p.PD, p.PD)           
  m.p.it[, m.M[, t] == "P" & df.X$tp > 9] <- c(1 - p.PD, 0, p.PD)
  
  # transition probabilities when dead  
  m.p.it[, m.M[, t] == "D"]  <- c(0, 0, 1)                                                           
  
  return(t(m.p.it))  # return transition probability
}  




start_time <- Sys.time()

for (t in 1:n.t) {
  
    # update tp after each cycle
    df.X$tp[m.M[, t] == 'P'] <- df.X$tp[m.M[, t] == 'P'] + 1
  
    ### inside Probs function ###
    m.p.it          <- matrix(0, nrow = n.s, ncol = n.i)  # create matrix of state transition probabilities
    rownames(m.p.it) <-  v.n                              # give the state names to the rows
    
    # determine eligibility of getting into pregnant state
    sex_female <- df.X$Sex == 'Female' 
    age_15_40 <- df.X$Age >= 15 & df.X$Age <= 40
    female_15_40 <- sex_female & age_15_40
    
    sex_male <- df.X$Sex == 'Male'
    age_not_15_40 <- df.X$Age < 15 | df.X$Age > 40
    cant_preg <- sex_male | age_not_15_40
    
    # lookup baseline probability and rate of dying based on individual characteristics
    p.HD_all <- p.mort1[df.X$Age+1,]
    p.HD_P     <- p.HD_all[m.M[, t] == "H" & female_15_40, "p.HD"]
    p.HD_noP  <- p.HD_all[m.M[, t] == "H" & cant_preg, "p.HD"]
    
    # update the v.p with the appropriate probabilitie, v.p will be fed to samplev 
    # (v.p is a n.i by n.s matrix, rows sum to 1)  
    
    # transition probabilities when healthy
    m.p.it[, m.M[, t] == "H" & female_15_40]  <- rbind(1 - p.HP - p.HD_P, p.HP, p.HD_P)   
    m.p.it[, m.M[, t] == "H" & cant_preg]  <- rbind(1 - p.HP - p.HD_noP, p.HP, p.HD_noP)
    
    # transition probabilities when pregnant
    m.p.it[, m.M[, t] == "P" & df.X$tp <= 9] <- c(0, 1 - p.PD, p.PD)           
    m.p.it[, m.M[, t] == "P" & df.X$tp > 9] <- c(1 - p.PD, 0, p.PD)
    
    # transition probabilities when dead  
    m.p.it[, m.M[, t] == "D"]  <- c(0, 0, 1)                                                           
    
    v.p <- t(m.p.it)
    
    ### outside Probs function, inside time loop ### 
    
    m.M[, t + 1]  <- samplev(v.p, 1)
    

    # add newborn baby to df.X
    baby.momIDs <- df.X$ID[m.M[, t] == "P" & df.X$tp > 9]
    if (length(baby.momIDs) > 0) {
    baby.IDs <- c(nrow(df.X)+1):c(nrow(df.X)+length(baby.momIDs))
    baby.Ages <- rep(0, length(baby.momIDs))
    baby.Sex <- sample(c(1,2), prob=c(0.5,0.5), replace=T, size=length(baby.momIDs))
    baby.tp <- rep(0, length(baby.momIDs))
    babies <- as.data.frame(cbind(baby.IDs, baby.Ages, baby.Sex, baby.momIDs,baby.tp))
    babies$baby.Sex[babies$baby.Sex==1] <- 'Female'
    babies$baby.Sex[babies$baby.Sex==2] <- 'Male'
    colnames(babies) <- colnames(df.X)
    df.X <- rbindlist(list(df.X, babies))
    
    # add newborn baby to m.M
    n.babies <- length(baby.momIDs)   # number of newborn babies in this cycle
    new_babies_m.M <- matrix(NA, nrow=n.babies, ncol=n.t+1)
    new_babies_m.M[ , t+1] <- rep('H', n.babies)  # at current cycle, newborns start in healthy state
    rownames(new_babies_m.M) <- paste('ind ', baby.IDs, sep='')
    m.M <- rbind(m.M, new_babies_m.M)
    # m.M[c(nrow(m.M)-n.babies+1):nrow(m.M), t+1] <- 'H'
    }
    
    # update n.i
    n.i <- nrow(df.X)
    n.i_hist <- c(n.i_hist, n.i)
    df.X$tp[df.X$tp > 9] <- 0
    
    # add one age
    if (t %% c(1/c.l) == 0) {df.X$Age <- df.X$Age + 1}  # increase age by 1 every 12 months

    # Display simulation progress
    if(t/(n.t/10) == round(t/(n.t/10), 0)) { # display progress every 10%
      cat('\r', paste(t/n.t * 100, "% done", sep = " "))
    }
}

end_time <- Sys.time()

end_time - start_time

# 1.693083 mins at 10,000

#### 06.2 Trace ####
# Plot the distribution of the population across healht states over time (Trace)

# count the number of individuals in each health state at each cycle
m.TR <- t(apply(m.M, 2, function(x) table(factor(x, levels = v.n, ordered = TRUE)))) 
m.TR <- m.TR / n.i_hist                                       # calculate the proportion of individuals 
colnames(m.TR) <- v.n                                    # name the rows of the matrix
rownames(m.TR) <- paste("Cycle", 0:n.t, sep = " ")       # name the columns of the matrix

# Plot trace of first health state
matplot(0:n.t, m.TR, type = 'l', 
        ylab = "Proportion of cohort",
        xlab = "Cycle",
        main = "Health state trace")                 
legend("topright", v.n, col = 1:n.s, lty = 1:n.s, bty = "n")  # add a legend to the graph

