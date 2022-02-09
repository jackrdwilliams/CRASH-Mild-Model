
## CRASH-Mild - Economic model ## 

if(!require(dplyr)) install.packages('dplyr')
if(!require(tidyr)) install.packages('tidyr')
if(!require(ggplot2)) install.packages('ggplot2')
if(!require(MCMCpack)) install.packages('MCMCpack')
if(!require(EnvStats)) install.packages('EnvStats')

library(tidyr)
library(dplyr)
library(ggplot2)
library(MCMCpack)
library(EnvStats)

save.results = TRUE

# as.integer(Sys.time())
set.seed(29716)

# Model Options

disc.c <- 0.035
disc.o <- 0.035

sims <- 10000

age <- 80
male = 0.5 
time.horizon = min(60, 100-ceiling(age))

lambda <- seq(from = 0, to = 50000, by = 100)

year1.cycle <- 1:12/12


discount.c <- matrix(1/(1+disc.c)^(0:time.horizon), time.horizon + 1, 2)
discount.o <- matrix(1/(1+disc.o)^(0:time.horizon), time.horizon + 1, 2)

## Age trace

## Clinical parameters

gen.clinical.characteristics <- function(){

# Mild ACM RR (CRASH-3)
treatment.effect <- 0.6972504
treatment.effect.log.se <- (log(1.076571) - log(0.45158))/(2*1.96)

treatment.effect.sims <- exp(rnorm(sims, log(treatment.effect), treatment.effect.log.se))

hi.risk <- 0.0208831
hi.risk.sims <- rtri(sims, min = 0.005, max = 0.05, mode = hi.risk)

non.hi.risk <- 0
non.hi.risk.sims <- rbeta(sims, 9, 920-9) * 0

smr.year1 <- 235.04/81.2
smr.year1.sims <- rnorm(sims, smr.year1, 0.301274691)

smr.year2 <- 61.47/39.45
smr.year2.sims <- rnorm(sims, smr.year2, 0.162178435) 

clin.names <- c("treatment.effect", "hi.risk", "non.hi.risk", "smr.year1", "smr.year2")


clin.char <- c(treatment.effect, hi.risk, non.hi.risk, smr.year1, smr.year2)
names(clin.char) <- clin.names

clin.char.sims <- data.frame(treatment.effect.sims, hi.risk.sims, non.hi.risk.sims, smr.year1.sims, smr.year2.sims)
colnames(clin.char.sims) <- clin.names

# Outcomes

outcome.names <- c("full", "good", "moderate", "severe", "vegetative")

full.placebo <- 0
good.placebo <- 208
moderate.placebo <- 49
severe.placebo <- 23
vegetative.placebo <- 5

dis.placebo <- c(full.placebo, good.placebo, moderate.placebo, severe.placebo, vegetative.placebo) 
disability.placebo <- dis.placebo / sum(dis.placebo)

  # full.txa <- 0
# good.txa <- 1516
# moderate.txa <- 701
# severe.txa <- 227
# vegetative.txa <- 22

#dis.txa <- c(full.txa, good.txa, moderate.txa, severe.txa, vegetative.txa)
#disability.txa <- dis.txa / sum(dis.txa)

disability.txa <- disability.placebo

names(disability.placebo) <- outcome.names
names(disability.txa) <- outcome.names

# Dirichlet distributions for outcomes
disability.placebo.sims <- MCMCpack::rdirichlet(n = sims , disability.placebo)
disability.txa.sims <- MCMCpack::rdirichlet(n = sims , disability.txa)

colnames(disability.placebo.sims) <- outcome.names
colnames(disability.txa.sims) <- outcome.names

return(list(clin.char,
            clin.char.sims,
            disability.placebo, 
            disability.txa, 
            disability.placebo.sims, 
            disability.txa.sims))
       
}


clin.char <- gen.clinical.characteristics()[[1]]
clin.char.sims <- gen.clinical.characteristics()[[2]]
disability.placebo <- gen.clinical.characteristics()[[3]]
disability.txa <- gen.clinical.characteristics()[[4]]
disability.placebo.sims <- gen.clinical.characteristics()[[5]]

## Assuming disability for TXA is equal to placebo
# disability.txa.sims <- gen.clinical.characteristics()[[6]]
disability.txa.sims <- disability.placebo.sims


# Adverse events 

gen.ae <- function(){
  
  # Risks from mild patients
  
  pe <- c(5, 1247)
  dvt <- c(1, 1251)
  stroke <- c(2, 1250)
  mi <- c(4, 1248)
  renal <- c(7, 1246)
  sepsis <- c(30, 1222)
  seizure <- c(15, 1237)
  gi <- c(6, 1246)
  
  ae <- data.frame(pe, dvt, stroke, mi, renal, sepsis, seizure, gi)
  
  # These RR are from all patients 

  pe.rr <-    c(0.9776573, 0.5093468, 1.876548)
  dvt.rr <-  c(1.222337, 0.572802, 2.608418)
  stroke.rr <- c(1.232966, 0.7144184, 2.127891)
  mi.rr <-    c(0.733402, 0.3093295, 1.738853)
  renal.rr <- c(1.275005, 0.9023189, 1.801623)
  sepsis.rr <- c(1.037228, 0.8853653, 1.215139)
  seizure.rr <- c(1.210695, 0.93925, 1.560589)
  gi.rr <-    c(0.7111777, 0.3740009, 1.352333)

  ae.rr <- data.frame(pe.rr, dvt.rr, stroke.rr, mi.rr, 
                      renal.rr, sepsis.rr, seizure.rr, gi.rr)  
  

  
  # Deterministic risks 
  
  placebo.risk <- ae[1,] / (ae[1,] + ae[2,])
  

  txa.risk  <- placebo.risk * as.numeric(ae.rr[1,])
  
  ae.risk <- data.frame(placebo.risk, txa.risk)
  
  # Simulations

  
  
  placebo.risk.sims <- matrix(NA, nrow = sims, ncol = length(placebo.risk))
  
  for(i in 1:length(placebo.risk)){
    placebo.risk.sims[,i] <- rbeta(sims, ae[1,i], ae[2,i]) 
  }
  
  
  log.rr <- as.numeric(log(ae.rr[1,]))
  log.se <- as.numeric((log(ae.rr[3,]) - log(ae.rr[2,]))/(2*1.96))
  
  rr.sims <- matrix(NA, nrow = sims, ncol = length(placebo.risk))
  
  for(i in 1:length(placebo.risk)){
    rr.sims[,i] <- exp(rnorm(sims, mean = log.rr[i] , sd = log.se[i])) 
  }
  
  
  txa.risk.sims <- placebo.risk.sims * rr.sims
  
  colnames(placebo.risk.sims) <- colnames(ae)
  colnames(txa.risk.sims) <- colnames(ae)

  return(list(placebo.risk, 
              txa.risk,
              placebo.risk.sims,
              txa.risk.sims))
  
}  

ae.placebo <- gen.ae()[[1]]
ae.txa <- gen.ae()[[2]]

ae.placebo.sims <- gen.ae()[[3]]
ae.txa.sims <- gen.ae()[[4]]



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



## Utility (now utility values only)

gen.utility <- function(){

utility.full <- 1
utility.good <- 0.894
utility.moderate <- 0.675
utility.severe <- 0.382
utility.vegetative <- -0.178

utility.values <- c(utility.full, utility.good, utility.moderate, utility.severe, utility.vegetative)

return(utility.values)

}


gen.utility.sims <- function(){
  
  utility.full <- rep(1, sims)
  utility.good <- rbeta(sims, 4331.027560, 513.522283)
  utility.moderate <- rbeta(sims, 2481.512500, 1194.802315)
  utility.severe <- rbeta(sims, 662.172521,	1071.263398)
  utility.vegetative <- rnorm(sims, -0.178, 0.0775672) ## needs normal distribution
  
  utility.values <- data.frame(utility.full, utility.good, utility.moderate, utility.severe, utility.vegetative)
  
# #  util.plac <- utility.values * matrix(dis.placebo, sims, length(dis.placebo), byrow= T)
#   util.plac <- utility.values * dis.placebo
#   util.values.plac <- apply(util.plac,  1,  sum )
#   
# #  util.txa <- utility.values * matrix(dis.txa, sims, length(dis.txa), byrow= T)
#   util.txa <- utility.values * dis.txa
#   util.values.txa <- apply(util.txa,  1,  sum)
  
  
  return(utility.values)
  
}


gen.utility.dec <- function(){

  # a <- data.frame(seq(0, 54, 1), 0)
  # b <- data.frame(seq(55,64, 1), 0.05)
  # c <- data.frame(seq(65, 74, 1), 0.07)
  # d <- data.frame(seq(75,100, 1), 0.12)
  # colnames(a) <- c("a","b")
  # colnames(b) <- c("a","b")
  # colnames(c) <- c("a","b")
  # colnames(d) <- c("a","b")
  # 
  # util.dec <- rbind(a, b, c, d)
  age.cycles <- c(50:100)
  est.util <- 0.9508566 + 0.0212126 * 0.5 - 0.0002587*age.cycles - 0.0000332*age.cycles^2
  est.dec <- est.util[1] - est.util 
  util.dec <- data.frame(age.cycle = age.cycles, dec = est.dec)
  
  utility.decrement <- util.dec %>% filter(age.cycle >= floor(age))
  
  return(utility.decrement)
  
} 



utility <- gen.utility()

utility.sims <- gen.utility.sims()

utility.decrement <- gen.utility.dec()


## Costs 

# Treatment costs 

gen.costs <- function(){

  inflation <- 314.915918 / c(301.150189, 282.5, 249.8)
  names(inflation) <- c("2018", "2012", "2007")
  
  cost.txa.dose <- 1.5
  cost.sodium <- 0 # 55p for 100ml, 270 for 500ml 
  cost.needle <- 0.12 # 12p needle syringe, £570 per year pre-drawn 
  cost.nurse <- 48 * (5 / 60) # 5 mins to draw TXA

  cost.treatment <- sum(cost.txa.dose, cost.sodium, cost.needle, cost.nurse)
  
  # Hospital costs

  los.placebo <-  4
  los.txa <- los.placebo
  prop.neuro <- 0.0345
  
  hospital.cost.initial <- 455.45033602 * inflation[1]
  hospital.cost.day <- 313.8757956 * inflation[1]
  neurosurgery.cost <- 7439.86435702 * inflation[1]

  hospital.cost.placebo <- hospital.cost.initial + hospital.cost.day * los.placebo + prop.neuro * neurosurgery.cost
  hospital.cost.txa <- hospital.cost.initial + hospital.cost.day * los.txa + prop.neuro * neurosurgery.cost


  # Monitoring costs 

  cost.good.st <- 240 * inflation[3]
  cost.moderate.st <- 17160  * inflation[3]
  cost.severe.st <- 33900 * inflation[3]

  cost.good.lt <- 24  * inflation[2]
  cost.moderate.lt <- 1600 * inflation[2]
  cost.severe.lt <- 12500 * inflation[2]


  # monitoring.costs.st <- sum(c(rep(cost.good.st,2), cost.moderate.st, rep(cost.severe.st, 2)) * dis.plac) / sum(dis.plac)
  # monitoring.costs.st <- sum(c(rep(cost.good.st,2), cost.moderate.st, rep(cost.severe.st, 2)) * dis.txa) / sum(dis.txa)
  # 
  # monitoring.costs.lt <- sum(c(rep(cost.good.lt,2), cost.moderate.lt, rep(cost.severe.lt, 2)) * dis.plac) / sum(dis.plac)
  # monitoring.costs.lt <- sum(c(rep(cost.good.lt,2), cost.moderate.lt, rep(cost.severe.lt, 2)) * dis.txa) / sum(dis.txa)
  
  ## Adverse events
  
  cost.pe <- 472.46526239 * inflation[1]
  cost.dvt <- 88.36794877 * inflation[1]
  cost.stroke <- 699.99139241 * inflation[1]
  cost.mi <- 1239.63378882 * inflation[1]
  cost.renal <- 444.09597767 * inflation[1]
  cost.sepsis <- 369.76900966 * inflation[1]
  cost.seizure <- 531.76481678 * inflation[1]
  cost.gi <- 347.38271318 * inflation[1]


  cost.names <- c("treatment", "hospital.placebo", "hospital.txa", "st.good", "st.mod", "st.sev", "lt.good", "lt.mod", "lt.sev", 
                  "pe.cost", "dvt.cost", "stroke.cost", "mi.cost", "renal.cost", "sepsis.cost", "seizure.cost", "gi.cost")
  costs <- c(cost.treatment, hospital.cost.placebo, hospital.cost.txa, cost.good.st, cost.moderate.st, cost.severe.st,
             cost.good.lt, cost.moderate.lt, cost.severe.lt, 
             cost.pe, cost.dvt, cost.stroke, cost.mi, cost.renal, cost.sepsis, cost.seizure, cost.gi)
  
  names(costs) <- cost.names
  
  cost.nurse.sims <- 48 * (runif(sims, 2, 8) / 60) 
  
  prob.neuro.sims <- rbeta(sims, 23.83398, 	667.00607)

  cost.treatment.sims <- rep(sum(cost.txa.dose, cost.sodium, cost.needle),sims) + cost.nurse.sims

  los.sims <- rgamma(sims, shape = 32, scale = 0.125) # same for both groups

  hospital.los.cost.placebo.sims <- los.sims * hospital.cost.day + hospital.cost.initial
  hospital.los.cost.txa.sims <- los.sims * hospital.cost.day + hospital.cost.initial

  hospital.cost.placebo.sims <- hospital.los.cost.placebo.sims * 1 + rep(neurosurgery.cost, sims) * prob.neuro.sims
  hospital.cost.txa.sims <- hospital.los.cost.txa.sims * 1 + rep(neurosurgery.cost, sims) * prob.neuro.sims
  
  # Monitoring costs # 

  cost.good.st.sims <- rgamma(sims, shape = (240^2 / 48^2), scale = 48^2 / 240 ) * inflation[3]
  cost.moderate.st.sims <- rgamma(sims, shape = (17160^2 / 3432^2), scale = 3432^2 / 17160 ) * inflation[3]
  cost.severe.st.sims <- rgamma(sims, shape = (33900^2 / 6780^2), scale = 6780^2 / 33900 ) * inflation[3]
  
  cost.good.lt.sims <- rgamma(sims, shape = (24^2 / 4.8^2), scale = 4.8^2 / 24 ) * inflation[2]
  cost.moderate.lt.sims <- rgamma(sims, shape = (1600^2 / 320^2), scale = 320^2 / 1600 ) * inflation[2]
  cost.severe.lt.sims <- rgamma(sims, shape = (12500^2 / 2500^2), scale = 2500^2 / 12500 ) * inflation[2]
  
  # Adverse event costs 
  
  pe.sims <- rep(cost.pe, sims) # cost.pe * runif(sims, 0.5, 1.5)
  dvt.sims <- rep(cost.dvt , sims)  # cost.dvt * runif(sims, 0.5, 1.5)
  stroke.sims <- rep(cost.stroke , sims)  # cost.stroke * runif(sims, 0.5, 1.5)
  mi.sims <- rep(cost.mi , sims)  # cost.mi * runif(sims, 0.5, 1.5)
  renal.sims <- rep(cost.renal , sims)  # cost.renal * runif(sims, 0.5, 1.5)
  sepsis.sims <- rep(cost.sepsis, sims)  # cost.sepsis * runif(sims, 0.5, 1.5)
  seizure.sims <- rep(cost.seizure , sims)  # cost.seizure * runif(sims, 0.5, 1.5)
  gi.sims <- rep(cost.gi , sims)  # cost.gi * runif(sims, 0.5, 1.5)
  
  costs.sims <- data.frame(cost.treatment.sims, hospital.cost.placebo.sims, hospital.cost.txa.sims,
                          cost.good.st.sims, cost.moderate.st.sims, cost.severe.st.sims,
                          cost.good.lt.sims, cost.moderate.lt.sims, cost.severe.lt.sims, 
                          pe.sims, dvt.sims, stroke.sims, mi.sims, renal.sims, sepsis.sims, seizure.sims, gi.sims)
  
  dim(costs.sims)
  colnames(costs.sims) <- cost.names
  
  return(list(costs,
              costs.sims))
  
}


costs <- gen.costs()[[1]]
costs.sims <- gen.costs()[[2]]



##  TRACE CALCULATIONS AND OUTCOMES  ## 

run.model <- function(clinical = clin.char, dis.plac = disability.placebo, dis.txa = disability.txa, util.values = utility, cost = costs, 
                      dec = utility.decrement, ae.p = ae.placebo, ae.t = ae.txa, d = discount.c, o = discount.o, output.type = 1) {
  

  for(i in 1:length(clinical)) assign(names(clinical[i]),clinical[i])
  
  trace.names <- c("alive placebo", "dead placebo", "alive txa", "dead txa")
  matrix <- matrix(0, 13, 4)
  matrix[1,] <- c(1,0,1,0)
  colnames(matrix) <- trace.names
  
  # Calculate the annual risk of death
  age.trace <- seq(from = age, to = age+1, by= 1/12)
  age.trace.floor <- floor(age.trace)
  a <- as.vector(acm %>% filter(age == unique(age.trace.floor)[[1]])) 
  b <- as.vector(acm %>% filter(age == unique(age.trace.floor)[[2]]))
  risk.year1 <- c( rep(a[[2]], sum(age.trace.floor==a[[1]])), rep(b[[2]], sum(age.trace.floor==b[[1]])))
  
  
  # Calculate the first year Markov trace
  matrix[2,2] <- matrix[1,2] + (hi.risk + non.hi.risk)
  matrix[2,1] <- 1 - matrix[2,2]
  matrix[2,4] <- matrix[1,4] + (hi.risk * treatment.effect + non.hi.risk)
  matrix[2,3] <- 1 - matrix[2,4]
  
  for(t in 3:13){
    matrix[t,2] <- matrix[t-1,2] + matrix[t-1,1] * (1 - exp(-((risk.year1[t]/12)*smr.year1)))
    matrix[t,1] <- 1 - matrix[t,2]
      
    matrix[t,4] <- matrix[t-1,4] + matrix[t-1,3] * (1 - exp(-((risk.year1[t]/12)*smr.year1)))
    matrix[t,3] <- 1 - matrix[t,4]
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
  
  trace <- list(trace.y1, trace.matrix)

  ## OUTCOMES 

  # Costs

  st.mon.costs <- c(cost[4], cost[5], cost[6])
  lt.mon.costs <- c(cost[7], cost[8], cost[9])
  
  st.mon.plac <- sum( st.mon.costs * c(dis.plac[1] + dis.plac[2], dis.plac[3], dis.plac[4] + dis.plac[5]))
  lt.mon.plac <- sum( lt.mon.costs * c(dis.plac[1] + dis.plac[2], dis.plac[3], dis.plac[4] + dis.plac[5]))

  st.mon.txa <- sum( st.mon.costs * c(dis.txa[1] + dis.txa[2], dis.txa[3], dis.txa[4] + dis.txa[5]))
  lt.mon.txa <- sum( lt.mon.costs * c(dis.txa[1] + dis.txa[2], dis.txa[3], dis.txa[4] + dis.txa[5]))
  
  # AE costs 
  
  cost.ae.p <- ae.p * cost[10:17]
  cost.ae.t <- ae.t * cost[10:17]
  cost.ae <- c(sum(cost.ae.p), sum(cost.ae.t))
  
  cost.matrix <- matrix(0, time.horizon + 1, 2) # placebo / txa
  
  cost.matrix[1,] <- cost[2:3] + cost.ae + (cost[1] * c(0,1)) 
  #cost.matrix[2,] <- trace[[2]][2,c(1,3)] * c(st.mon.plac, st.mon.txa)
  cost.matrix[2,] <- apply(trace.y1[2:13, c(1,3)], 2, mean) * c(st.mon.plac, st.mon.txa)
  
  cost.matrix[3:(time.horizon+1),] <- trace[[2]][3:(time.horizon+1),c(1,3)] * c(lt.mon.plac, lt.mon.txa)
  
  cost.matrix.d <- cost.matrix * d
  
  cost.sum <- apply(cost.matrix.d, 2, sum)
  
  
  # Utility

    #util.plac <- utility.values * matrix(dis.placebo, sims, length(dis.placebo), byrow= T)
  util.plac <- util.values * dis.plac
  util.values.plac <- sum(util.plac)
     
  #  util.txa <- utility.values * matrix(dis.txa, sims, length(dis.txa), byrow= T)
  util.txa <- util.values * dis.txa
  util.values.txa <- sum(util.txa)
  
  utility.matrix <- matrix(0, time.horizon + 1, 2) # placebo / txa
  
  # utility.matrix[2,1] <- mean(trace[[1]][2:13,1]) * (util.values.plac - dec[2,2])
  # utility.matrix[2,2] <- mean(trace[[1]][2:13,3]) * (util.values.txa - dec[2,2])
  
  dec.year1 <- dec[findInterval(age + year1.cycle, dec[,1]),2]

  utility.matrix[2,1] <- sum(trace[[1]][2:13,1] * (util.values.plac - dec.year1)/12)
  utility.matrix[2,2] <- sum(trace[[1]][2:13,3] * (util.values.txa -  dec.year1)/12)
  
  utility.matrix[3:(time.horizon+1),1] <- trace[[2]][3:(time.horizon+1),1] * (util.values.plac - dec[3:(time.horizon+1),2]) 
  utility.matrix[3:(time.horizon+1),2] <- trace[[2]][3:(time.horizon+1),3] * (util.values.txa - dec[3:(time.horizon+1),2])
  
  utility.matrix.d <- utility.matrix * o 
  
  utility.sum <- apply(utility.matrix.d, 2, sum) 
  
  
  ## ICER ## 
  
  icer <-   (cost.sum[2] - cost.sum[1]) / (utility.sum[2] - utility.sum[1])  
  
  if(output.type==1) {
    return(list(c(cost.placebo = cost.sum[1], 
                utility.placebo = utility.sum[1], 
                cost.txa = cost.sum[2], 
                utility.txa = utility.sum[2]), 
              icer = icer)) 
  } else if(output.type==9){
      return(trace)
  } else return(icer)
  
  
}



## Deterministic results 
run.model(clin.char, dis.plac = disability.placebo, dis.txa = disability.txa, util.values = utility,
          cost = costs, dec = utility.decrement, ae.p = ae.placebo, ae.t = ae.txa, d = discount.c, o = discount.o)

