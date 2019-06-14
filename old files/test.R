### development


### individual i, cycle 0

p.it <- vector("numeric", length(v.n)) 

m.M[1, 1] <- v.M_Init[1]       # initial health state for individual i

m.M[1,1]

c.it <- c()
c.it[m.M[1, 1] == "dead"]    <- c.D     # costs at dead state
c.it[m.M[1, 1] == "healthy"] <- c.H     # costs accrued by being healthy this cycle
c.it[m.M[1, 1] == "sick"]    <- c.S     # costs accrued by being sick this cycle

q.it <- c() 
q.it[m.M[1, 1] == "dead"]    <- u.D        # QALYs at dead state
q.it[m.M[1, 1] == "healthy"] <- u.H        # QALYs accrued by being healthy this cycle
q.it[m.M[1, 1] == "sick"]    <- u.S        # QALYs accrued by being sick this cycle

m.C[1, 1] <- c.it
m.E[1, 1] <- q.it

p.it[m.M[1, 1] == "healthy"] <- c(1 - p.HD - p.HS, p.HS, p.HD)     # transition probabilities when healthy 
p.it[m.M[1, 1] == "sick"]    <- c(0, 1 - p.SD, p.SD)               # transition probabilities when sick 
p.it[m.M[1, 1] == "dead"]    <- c(0, 0, 1)                         # transition probabilities when dead 
p.it

v.p <- p.it


### individual i, cycle 1

v.n
v.p

m.M[1,2] <- sample(x = v.n, prob = v.p, size = 1)
m.M[1,2]

c.it <- c()
c.it[m.M[1, 2] == "dead"]    <- c.D     # costs at dead state
c.it[m.M[1, 2] == "healthy"] <- c.H     # costs accrued by being healthy this cycle
c.it[m.M[1, 2] == "sick"]    <- c.S     # costs accrued by being sick this cycle

q.it <- c() 
q.it[m.M[1, 2] == "dead"]    <- u.D        # QALYs at dead state
q.it[m.M[1, 2] == "healthy"] <- u.H        # QALYs accrued by being healthy this cycle
q.it[m.M[1, 2] == "sick"]    <- u.S        # QALYs accrued by being sick this cycle

m.C[1, 2] <- c.it
m.E[1, 2] <- q.it

# check cycles 0 and 1
m.M
m.C
m.E


####################################### now work on extended model #################################################

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
p.HHB <- 0.01        # probability healthy -> healthy with baby
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

# i=1, n.t=0

m.M[1, 1] <- v.M_Init[1]       # 1n1t1al health state for 1nd1v1dual 1
m.C[1, 1] <- Costs(m.M[1, 1])  # costs accrued 1nd1v1dual 1 dur1ng cycle 0
m.E[1, 1] <- Effs(m.M[1, 1])   # QALYs accrued 1nd1v1dual 1 dur1ng cycle 0

v.p <- Probs(m.M[1, 1])

# sample the current health state based on trans1t1on probab1l1t1es v.p 
m.M[1, 1 + 1] <- sample(x = v.n, prob = v.p, size = 1)
m.M[1, 1 + 1] <- 'healthy with baby'


if (m.M[1, 1 + 1] == "healthy with baby") {
  m.M <- rbind(m.M, rep(NA, ncol(m.M)))  # add the baby to the m.M matrix
  m.M[nrow(m.M), 1 + 1] <- 'healthy'     # baby begins in the healthy state
  m.C <- rbind(m.C, rep(0, ncol(m.C)))  # add the baby to the m.C matrix, previous costs assigned to 0
  m.C[nrow(m.C), 1 + 1] <- Costs(m.M[1, 1])    # cost in the healthy state
  m.E <- rbind(m.E, rep(0, ncol(m.E)))  # add the baby to the m.M matrix, previous utilities assigned to 0
  m.E[nrow(m.E), 1 + 1] <- Effs(m.M[1, 1])     # utility in the healthy state
  # add rownames
  rownames(m.M)[nrow(m.M)] <- rownames(m.C)[nrow(m.C)] <- rownames(m.E)[nrow(m.E)] <- 
    paste("ind ", as.character(nrow(m.M)), ",", " baby", sep = "")
}

# calculate costs and effects accrued by 1nd1v1dual 1 for cycle t
m.C[1, 1 + 1] <- Costs(m.M[1, 1 + 1])   # costs
m.E[1, 1 + 1] <- Effs(m.M[1, 1 + 1])    # QALYs

i <- 1
while (i <= 6) {
  print(i)
  i <- i + 1
}
