rm(list = ls())  # Delete everything that is in R's memory

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))  #set working directory as the folder where the course material is stored

#### 01 Load packages ####
library(plyr)
library(reshape2)

#### 02 Load Functions ####
#### samplev solution ####
# samplev : efficient implementation of the rMultinom() function of the Hmisc package 
# probs:    the matrix of probabilities for each individual for each state at each time point
# m:      the number of values that need to be sampled at a time per individual

samplev <- function(m.Probs, m) {
  # Arguments
  # m.Probs: matrix with probabilities (n.i * n.s)
  # m:       number of states than need to be sampled per individual  
  # Return
  # ran:    n.i x m matrix filled with sampled health state(s) per individual
  
  d <- dim(m.Probs)  # dimensions of the matrix filled with the multinomical probabilities for the health states 
  n <- d[1]          # first dimension - number of rows (number of individuals to sample for)
  k <- d[2]          # second dimension - number of columns (number of health states considered)
  lev <- dimnames(m.Probs)[[2]]  # extract the names of the health states considered for sampling
  if (!length(lev))  # in case names for the health states are missing, use numbers to specify the health states
    lev <- 1:k       # create a sequence from 1:k (number of health states considered)
  # create a matrix 
  ran <- matrix(lev[1], ncol = m, nrow = n) # create the matrix ran, filled with the first health state of the levels 
  U <- t(m.Probs)    # transposed m.Probs matrix n.i x n.s --> n.s x n.i 
  
  for(i in 2:k) {    # start loop, from the 2nd health states
    U[i, ] <- U[i, ] + U[i - 1, ] # start summing the probabilities of the different health states per individual 
  }
  if (any((U[k, ] - 1) > 1e-05))  # sum of all probs per individual - 1 should be 0 (use 1e-05 for rounding issues), else print the error statement
    stop("error in multinom: probabilities do not sum to 1")
  
  for (j in 1:m) {   # start loop of the state that needs to be sampled (m)
    un <- rep(runif(n), rep(k, n))       # sample from a uniform distribution of length n*k
    ran[, j] <- lev[1 + colSums(un > U)] # store the health state at the jth column of the U matrix
  }
  ran # return the new health state per individual n.i x m
} # close the function 


#### 03 Input Model Parameters ####
set.seed(1987)                            # set the seed 

# Model structure
v.n <- c("H", "P", "D")       # vector with state names
n.s <- length(v.n)                        # number of states
c.l <- 1/12                                 # cycle length
max.age <- 40                            # maximum age
# n.t <- max.age / c.l                           # number of cycles
n.t <- 15
n.i <- 10                              # number of individuals
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
p.HP <- 0.05        # probability healthy -> pregnant
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

### inside Probs function ###

m.p.it          <- matrix(0, nrow = n.s, ncol = n.i)  # create matrix of state transition probabilities
rownames(m.p.it) <-  v.n                              # give the state names to the rows

### create and update M_it <- m.M[, t]
M_it <- m.M[,1]

# determine eligibility of getting into pregnant state
sex_female <- df.X$Sex == 'Female' 
age_15_40 <- df.X$Age >= 15 & df.X$Age <= 40
female_15_40 <- sex_female & age_15_40

# determine eligibilty of giving birth (track time in pregancy)
df.X$tp <- rep(0, nrow(df.X))

# lookup baseline probability and rate of dying based on individual characteristics
p.HD_all <- join(df.X, p.mort1, by = c("Sex", "Age"))
p.HD_P     <- p.HD_all[M_it == "H" & female_15_40, "p.HD"]
p.HD_noP  <- p.HD_all[M_it == "H" & !female_15_40, "p.HD"]

# update the v.p with the appropriate probabilitie, v.p will be fed to samplev 
# (n.i by n.s matrix, rows sum to 1)  

# transition probabilities when healthy
m.p.it[, M_it == "H" & female_15_40]  <- rbind(1 - p.HP - p.HD_P, p.HP, p.HD_P)   
m.p.it[, M_it == "H" & !female_15_40]  <- rbind(1 - p.HP - p.HD_noP, p.HP, p.HD_noP)

# transition probabilities when pregnant
m.p.it[, M_it == "P" & df.X$tp <= 9] <- c(0, 1 - p.PD, p.PD)           
m.p.it[, M_it == "P" & df.X$tp > 9] <- c(1 - p.PD, 0, p.PD)

# transition probabilities when dead  
m.p.it[, M_it == "D"]  <- c(0, 0, 1)                                                           

t(m.p.it)

### outside Probs function, inside time loop ###

# update tp after each cycle
# a <- 1:5
# b <- c(T,F,T,F,T)
# a[b] <- a[b] + 10
df.X$tp[m.M[, t] == 'P'] <- df.X$tp[m.M[, t] == 'P'] + 1
df.X$tp[df$X.tp == 9] <- 0

# add newborn baby to df.X

# test
M_it[1:2] <- 'P'
df.X$tp[1:2] <- 10

baby.momIDs <- df.X$ID[M_it == "P" & df.X$tp > 9]
baby.IDs <- c(nrow(df.X)+1):c(nrow(df.X)+length(baby.momIDs))
baby.Ages <- rep(0, length(baby.momIDs))
baby.Sex <- sample(c(1,2), prob=c(0.5,0.5), size=length(baby.momIDs))
baby.tp <- rep(0, length(baby.momIDs))
babies <- as.data.frame(cbind(baby.IDs, baby.Ages, baby.Sex, baby.momIDs,baby.tp))
babies$baby.Sex[babies$baby.Sex==1] <- 'Female'
babies$baby.Sex[babies$baby.Sex==2] <- 'Male'
colnames(babies) <- colnames(df.X)
df.X <- rbind(df.X, babies)

# add newborn baby to m.M
n.babies <- length(baby.momIDs)   # number of newborn babies in this cycle
new_babies_m.M <- matrix(NA, nrow=n.babies, ncol=n.t+1)
new_babies_m.M[,1] <- rep('H', n.babies)  # at current cycle, newborns start in healthy state
rownames(new_babies_m.M) <- paste('ind ', baby.IDs, sep='')
m.M <- rbind(m.M, new_babies_m.M)

# update n.i
n.i <- nrow(df.X)

# add one age
if (t %% c(1/c.l) == 0) {df.X$Age <- df.X$Age + 1}  # increase age by 1 every 12 months







