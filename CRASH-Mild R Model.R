
## CRASH-Mild - Economic model ## 



if(!require(dplyr)) install.packages('dplyr')
library(dplyr)



# Options

disc.c <- 0.035
disc.o <- 0.035

sims <- 10

age <- 57.73704
time.horizon = min(60, 100-ceiling(age))

## Age trace



## Clinical parameters

gen.clinical.characteristics <- function(){

treatment.effect <- 0.6672781
treatment.effect.sims <- exp(rnorm(sims, log(treatment.effect), 0.2472015))

hi.risk <- 13/440
hi.risk.sims <- rbeta(sims, 13, 440-13)

non.hi.risk <- 9/920
non.hi.risk.sims <- rbeta(sims, 9, 920-9)

smr.year1 <- 4
smr.year1.sims <- rnorm(sims, smr.year1, 0.41632922)

smr.year2 <- 30.99/13.72
smr.year2.sims <- rnorm(sims, smr.year2, 0.2350955) 

clin.names <- c("treatment.effect", "hi.risk", "non.hi.risk", "smr.year1", "smr.year2")


clin.char <- c(treatment.effect, hi.risk, non.hi.risk, smr.year1, smr.year2)
names(clin.char) <- clin.names

clin.char.sims <- data.frame(treatment.effect.sims, hi.risk.sims, non.hi.risk.sims, smr.year1.sims, smr.year2.sims)
colnames(clin.char.sims) <- clin.names

# Outcomes

outcome.names <- c("full", "good", "moderate", "severe", "vegetative")

full.placebo <- 0
good.placebo <- 1516
moderate.placebo <- 701
severe.placebo <- 227
vegetative.placebo <- 22

disability.placebo <- c(full.placebo, good.placebo, moderate.placebo, severe.placebo, vegetative.placebo)

full.txa <- 0
good.txa <- 1516
moderate.txa <- 701
severe.txa <- 227
vegetative.txa <- 22

disability.txa <- c(full.txa, good.txa, moderate.txa, severe.txa, vegetative.txa)

names(disability.placebo) <- outcome.names
names(disability.txa) <- outcome.names

# Dirichlet distributions for outcomes
disability.placebo.sims <- MCMCpack::rdirichlet(n = sims , disability.placebo)
disability.txa.sims <- MCMCpack::rdirichlet(n = sims , disability.placebo)

colnames(disability.placebo.sims) <- outcome.names
colnames(disability.txa.sims) <- outcome.names

return(list(clin.char,
            clin.char.sims,
            disability.placebo, 
            disability.txa, 
            disability.placebo.sims, 
            disability.txa.sims))
       
}


gen.clinical.characteristics()


clin.char <- gen.clinical.characteristics()[[1]]
clin.char.sims <- gen.clinical.characteristics()[[2]]
disability.placebo <- gen.clinical.characteristics()[[3]]
disability.txa <- gen.clinical.characteristics()[[4]]
disability.placebo.sims <- gen.clinical.characteristics()[[5]]
disability.txa.sims <- gen.clinical.characteristics()[[6]]



## Generate long-term risk of death

gen.acm <- function(male = 0.5){
  
  lifetable <- read.csv("inputs/ONS_life_tables_2017-19.csv", header=TRUE)
  lifetable <- as.data.frame(lifetable[c(1,2,8)])
  av <- lifetable
  av[,2] <- lifetable[,2] * male + lifetable[,3] * (1 - male)
  av <- av[,1:2]

  acm <- av %>% filter(x >= floor(age))
  colnames(acm) <- c("age", "rate")
  

  
  return(acm)
  
}

acm <- gen.acm()






## HRQoL

gen.utility <- function(dis.placebo = disability.placebo, dis.txa = disability.txa){

utility.full <- 1
utility.good <- 0.894
utility.moderate <- 0.675
utility.severe <- 0.382
utility.vegetative <- -0.178

utility.values <- c(utility.full, utility.good, utility.moderate, utility.severe, utility.vegetative)

utility.placebo <- sum(utility.values*dis.placebo)/sum(dis.placebo)
utility.txa <- sum(utility.values*dis.txa)/sum(dis.txa)

return(c(placebo = utility.placebo, 
       txa = utility.txa))

}

utility <- gen.utility()

gen.utility.sims <- function(dis.placebo = disability.placebo, dis.txa = disability.txa){
  
  utility.full <- rep(1, sims)
  utility.good <- rbeta(sims, 49.9585894, 5.9235016)
  utility.moderate <- rbeta(sims, 30.538, 14.7034815)
  utility.severe <- rbeta(sims, 10.9303087,	17.6830648)
  utility.vegetative <- rnorm(sims, -0.178, 0.19) ## needs normal distribution
  
  utility.values <- data.frame(utility.full, utility.good, utility.moderate, utility.severe, utility.vegetative)
  
  
  util.plac <- utility.values * matrix(dis.placebo, sims, length(dis.placebo), byrow= T)
  util.values.plac <- apply(util.plac,  1,  function(x) sum(x) / sum(dis.placebo))
  
  util.txa <- utility.values * matrix(dis.txa, sims, length(dis.txa), byrow= T)
  util.values.txa <- apply(util.txa,  1,  function(x) sum(x) / sum(dis.txa))
  
  
  return(data.frame(placebo = util.values.plac, 
                    txa = util.values.txa))
  
}

utility.sims <- gen.utility.sims()


## need to add utility decrements 

a <- data.frame(seq(0, 54, 1), 0)

b <- data.frame(seq(55,64, 1), 0.05)
c <- data.frame(seq(65, 74, 1), 0.07)
d <- data.frame(seq(75,100, 1), 0.12)
colnames(a) <- c("a","b")
colnames(b) <- c("a","b")
colnames(c) <- c("a","b")
colnames(d) <- c("a","b")

utility.decrement <- rbind(a, b, c, d)



## Costs 

# Treatment costs 

cost.txa.dose <- 6
cost.sodium <- 0.55 + 2.7 # 55p for 100ml, 270 for 500ml 
cost.needle <- 0.05 
cost.nurse <- 12.95

cost.treatment <- sum(cost.txa.dose, cost.sodium, cost.needle, cost.nurse)


# Hospital costs

los.placebo <-  12.44828
los.txa <- 12.44828

hospital.cost.initial <- 455.4503
hospital.cost.day <- 313.8758

hospital.cost.placebo <- hospital.cost.initial + hospital.cost.day * los.placebo
hospital.cost.txa <- hospital.cost.initial + hospital.cost.day * los.txa


# Monitoring costs 

cost.good.st <- 290.145026
cost.moderate.st <- 20745.36936
cost.severe.st <- 40982.984929

cost.good.lt <- 25.65601
cost.moderate.lt <- 1710.40065
cost.severe.lt <- 13362.505071


monitoring.costs.st <- sum(c(rep(cost.good.st,2), cost.moderate.st, rep(cost.severe.st, 2)) * disability.placebo) / sum(disability.placebo)
monitoring.costs.st <- sum(c(rep(cost.good.st,2), cost.moderate.st, rep(cost.severe.st, 2)) * disability.txa) / sum(disability.txa)

monitoring.costs.lt <- sum(c(rep(cost.good.lt,2), cost.moderate.lt, rep(cost.severe.lt, 2)) * disability.placebo) / sum(disability.placebo)
monitoring.costs.lt <- sum(c(rep(cost.good.lt,2), cost.moderate.lt, rep(cost.severe.lt, 2)) * disability.txa) / sum(disability.txa)



costs <- c(cost.treatment, hospital.cost.placebo, hospital.cost.txa, monitoring.costs.st, monitoring.costs.lt)



## Year 1 trace



# Generate year 1 risk of death 

gen.trace <- function(clinical){

  # this unpacks each element of the clinical parameter to use in calculations below
  for(i in 1:length(clinical)) assign(names(clinical[i]),clinical[i])
  
  trace.names <- c("alive placebo", "dead placebo", "alive txa", "dead txa")
  matrix <- matrix(0, 366, 4)
  matrix[1,] <- c(1,0,1,0)
  colnames(matrix) <- trace.names

  # Calculate the annual risk of death, and make sure it changes at correct point 
  age.trace <- seq(from = age, to = age+1, by= 1/365)
  age.trace.floor <- floor(age.trace)
  a <- as.vector(acm %>% filter(age == unique(age.trace.floor)[[1]])) 
  b <- as.vector(acm %>% filter(age == unique(age.trace.floor)[[2]]))
  risk.year1 <- c( rep(a[[2]], sum(age.trace.floor==a[[1]])), rep(b[[2]], sum(age.trace.floor==b[[1]])))


# Calculate the first year Markov trace (first 28 days vs. rest of year below)

for(t in 2:366){
  
  if(t<=28+1){
  
  matrix[t,2] <- matrix[t-1,2] + ((hi.risk + non.hi.risk)/28)
  matrix[t,1] <- 1 - matrix[t,2]
  
  matrix[t,4] <- matrix[t-1,4] + ((hi.risk * treatment.effect + non.hi.risk)/28)
  matrix[t,3] <- 1 - matrix[t,4]
  
  } else {
    
  matrix[t,2] <- matrix[t-1,2] + matrix[t-1,1] * (1 - exp(-((risk.year1[t]/365)*smr.year1)))
  matrix[t,1] <- 1 - matrix[t,2]
  
  matrix[t,4] <- matrix[t-1,4] + matrix[t-1,3] * (1 - exp(-((risk.year1[t]/365)*smr.year1)))
  matrix[t,3] <- 1 - matrix[t,4]
  }
  
}

  trace.y1 <- matrix



  # The main trace (using estimates from first year trace)

  trace.matrix <- matrix(0, time.horizon + 1, 4) 
  colnames(trace.matrix) <- trace.names
  trace.matrix[1,] <- c(1,0,1,0)
  trace.matrix[2,] <- as.numeric(tail(trace.y1, 1))
  
  for(i in 3:(time.horizon + 1)){
    trace.matrix[i,2] <- trace.matrix[i-1,2] + trace.matrix[i-1,1] * (1 - exp(-(acm[i,2] * smr.year2)))
    trace.matrix[i,1] <- 1 - trace.matrix[i,2] 
    trace.matrix[i,4] <- trace.matrix[i-1,4] + trace.matrix[i-1,3] * (1 - exp(-(acm[i,2] * smr.year2)))
    trace.matrix[i,3] <- 1 - trace.matrix[i,4] 
  }

  return(list(trace.y1, 
            trace.matrix))

}


# trace with deterministic
trace.results <- gen.trace(clin.char)



# trace with probabilistic simulation
 # gen.trace(unlist(clin.char.sims[1,]))


trace[[2]]

write.excel <- function(x,row.names=FALSE,col.names=TRUE,...) {
  write.table(x,"clipboard",sep="\t",row.names=row.names,col.names=col.names,...)
}

write.excel(trace[[2]])
 

# discount.costs = disc.c, discount.outcomes = disc.o





gen.outcomes <- function(trace, util = utility, cost = costs){
  
  # Costs
  
  # for(i in 1:length(cost)) assign(names(cost[i]),cost[i])
  
  
  
  # Utility
  
  utility.matrix <- matrix(0, time.horizon + 1, 2) # placebo / txa
  utility.matrix[2,] <- 1 # calculate utility based on first year trace
  
  utility.matrix[3:(time.horizon+1),1] <- trace[[2]][3:(time.horizon+1),1] * util[1]
  utility.matrix[3:(time.horizon+1),2] <- trace[[2]][3:(time.horizon+1),3] * util[2]
  
  utility.sum <- apply(utility.matrix, 2, sum) 
  
  return(list(utility.matrix, 
              utility.sum))
  
}



z <- round(gen.outcomes(trace.results)[[1]],5)

write.excel(z)




## Caclulating payoffs 




