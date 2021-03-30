
## CRASH-Mild - Economic model ## 

if(!require(dplyr)) install.packages('dplyr')
if(!require(tidyr)) install.packages('tidyr')
if(!require(ggplot2)) install.packages('ggplot2')

library(tidyr)
library(dplyr)
library(ggplot2)

# Model Options

disc.c <- 0.035
disc.o <- 0.035

sims <- 1000
outer.loops <- 100
inner.loops <- 100

age <- 70
male = 0.5 # plaecholder
time.horizon = min(60, 100-ceiling(age))



## Age trace

## Clinical parameters

gen.clinical.characteristics <- function(){

# Mild ACM RR (CRASH-3)
treatment.effect <- 0.6972504
treatment.effect.sims <- exp(rnorm(sims, log(treatment.effect), 0.2216285))

hi.risk <- 0.0142
hi.risk.sims <- rbeta(sims, 16.74564893, 1162.525402)

non.hi.risk <- 0
non.hi.risk.sims <- rbeta(sims, 9, 920-9) * 0

smr.year1 <- 235.04/81.2
smr.year1.sims <- rnorm(sims, smr.year1, 0.41632922)

smr.year2 <- 61.47/39.45
smr.year2.sims <- rnorm(sims, smr.year2, 0.2350955) 

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

disability.placebo <- c(full.placebo, good.placebo, moderate.placebo, severe.placebo, vegetative.placebo)

# full.txa <- 0
# good.txa <- 1516
# moderate.txa <- 701
# severe.txa <- 227
# vegetative.txa <- 22

#disability.txa <- c(full.txa, good.txa, moderate.txa, severe.txa, vegetative.txa)
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



## Utility

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


gen.utility.sims <- function(dis.placebo = disability.placebo, dis.txa = disability.txa){
  
  utility.full <- rep(1, sims)
  utility.good <- rbeta(sims, 49.9585894, 5.9235016)
  utility.moderate <- rbeta(sims, 30.538, 14.7034815)
  utility.severe <- rbeta(sims, 10.9303087,	17.6830648)
  utility.vegetative <- rnorm(sims, -0.178, 0.19) ## needs normal distribution
  
  utility.values <- data.frame(utility.full, utility.good, utility.moderate, utility.severe, utility.vegetative)
  
  
#  util.plac <- utility.values * matrix(dis.placebo, sims, length(dis.placebo), byrow= T)
  util.plac <- utility.values * dis.placebo
  util.values.plac <- apply(util.plac,  1,  sum )
  
#  util.txa <- utility.values * matrix(dis.txa, sims, length(dis.txa), byrow= T)
  util.txa <- utility.values * dis.txa
  util.values.txa <- apply(util.txa,  1,  sum)
  
  
  return(data.frame(placebo = util.values.plac, 
                    txa = util.values.txa))
  
}


gen.utility.dec <- function(){

  a <- data.frame(seq(0, 54, 1), 0)
  b <- data.frame(seq(55,64, 1), 0.05)
  c <- data.frame(seq(65, 74, 1), 0.07)
  d <- data.frame(seq(75,100, 1), 0.12)
  colnames(a) <- c("a","b")
  colnames(b) <- c("a","b")
  colnames(c) <- c("a","b")
  colnames(d) <- c("a","b")

  util.dec <- rbind(a, b, c, d)

  utility.decrement <- util.dec %>% filter(a >= floor(age))
  
  return(utility.decrement)
  
} 



utility <- gen.utility()

utility.sims <- gen.utility.sims(disability.placebo.sims, disability.placebo.sims)

utility.decrement <- gen.utility.dec()



## Costs 

# Treatment costs 

## need to add tx effect costs

gen.costs <- function(dis.plac, dis.txa, dis.plac.sims, dis.txa.sims){

  cost.txa.dose <- 6
  cost.sodium <- 0.55 + 2.7 # 55p for 100ml, 270 for 500ml 
  cost.needle <- 0.05151 
  cost.nurse <- 12.95

  cost.treatment <- sum(cost.txa.dose, cost.sodium, cost.needle, cost.nurse)
  

  # Hospital costs

  los.placebo <-  4
  los.txa <- los.placebo
  prop.neuro <- 0.0345
  
  hospital.cost.initial <- 455.4503
  hospital.cost.day <- 313.8758
  neurosurgery.cost <- 7439.864357

  hospital.cost.placebo <- hospital.cost.initial + hospital.cost.day * los.placebo + prop.neuro * neurosurgery.cost
  hospital.cost.txa <- hospital.cost.initial + hospital.cost.day * los.txa + prop.neuro * neurosurgery.cost


  # Monitoring costs 

  cost.good.st <- 290.145026
  cost.moderate.st <- 20745.36936
  cost.severe.st <- 40982.984929

  cost.good.lt <- 25.65601
  cost.moderate.lt <- 1710.40065
  cost.severe.lt <- 13362.505071


  monitoring.costs.st <- sum(c(rep(cost.good.st,2), cost.moderate.st, rep(cost.severe.st, 2)) * dis.plac) / sum(dis.plac)
  monitoring.costs.st <- sum(c(rep(cost.good.st,2), cost.moderate.st, rep(cost.severe.st, 2)) * dis.txa) / sum(dis.txa)

  monitoring.costs.lt <- sum(c(rep(cost.good.lt,2), cost.moderate.lt, rep(cost.severe.lt, 2)) * dis.plac) / sum(dis.plac)
  monitoring.costs.lt <- sum(c(rep(cost.good.lt,2), cost.moderate.lt, rep(cost.severe.lt, 2)) * dis.txa) / sum(dis.txa)


  cost.names <- c("treatment", "hospital.placebo", "hospital.txa", "st.mon", "lt.mon")
  costs <- c(cost.treatment, hospital.cost.placebo, hospital.cost.txa, monitoring.costs.st, monitoring.costs.lt)
  names(costs) <- cost.names
  

  # Simulations - to be input later on (PLACEHOLDER)
  
  prob.neuro.sims <- rbeta(sims, 23.83398, 	667.00607)
  
  cost.treatment.sims <- rep(cost.treatment, sims)  * 1 ## PLACEHOLDER
  hospital.cost.placebo.sims <- rep(hospital.cost.placebo, sims) * 1 + rep(neurosurgery.cost, sims) * prob.neuro.sims  ## PLACEHOLDER
  hospital.cost.txa.sims <- rep(hospital.cost.txa, sims) * 1 + rep(neurosurgery.cost, sims) * prob.neuro.sims ## PLACEHOLDER
  
  # Monitoring costs # 
  
  inf0607 <- 302/249.8
  
  cost.good.st.sims <- rgamma(sims, shape = (240^2 / 48^2), scale = 48^2 / 240 ) * inf0607
  cost.moderate.st.sims <- rgamma(sims, shape = (17160^2 / 3432^2), scale = 3432^2 / 17160 ) * inf0607
  cost.severe.st.sims <- rgamma(sims, shape = (33900^2 / 6780^2), scale = 6780^2 / 33900 ) * inf0607
  
  cost.good.lt.sims <- rgamma(sims, shape = (24^2 / 4.8^2), scale = 4.8^2 / 24 ) * inf0607
  cost.moderate.lt.sims <- rgamma(sims, shape = (1600^2 / 320^2), scale = 320^2 / 1600 ) * inf0607
  cost.severe.lt.sims <- rgamma(sims, shape = (12500^2 / 2500^2), scale = 2500^2 / 12500 ) * inf0607
  
  df.st <- data.frame(cost.good.st.sims, cost.good.st.sims, cost.moderate.st.sims, cost.severe.st.sims, cost.severe.st.sims) 
  df.lt <- data.frame(cost.good.lt.sims, cost.good.lt.sims, cost.moderate.lt.sims, cost.severe.lt.sims, cost.severe.lt.sims)
  
  
  monitoring.costs.st <- df.st * dis.plac.sims # These are currently the same but structured in case these need to differ
  monitoring.costs.st <- df.st * dis.txa.sims  # These are currently the same but structured in case these need to differ
  
  m.costs.st <- apply(monitoring.costs.st, 1, sum)
  
  monitoring.costs.lt <- df.lt * dis.plac.sims  # These are currently the same but structured in case these need to differ
  monitoring.costs.lt <- df.lt * dis.txa.sims  # These are currently the same but structured in case these need to differ
  
  m.costs.lt <- apply(monitoring.costs.lt, 1, sum)
  
  costs.sims <- data.frame(cost.treatment.sims, hospital.cost.placebo.sims, hospital.cost.txa.sims, 
                           m.costs.st, m.costs.lt)
  
  colnames(costs.sims) <- cost.names
  
  return(list(costs,
              costs.sims))
  
}


costs <- gen.costs(disability.placebo, disability.txa, disability.placebo.sims, disability.txa.sims)[[1]]
costs.sims <- gen.costs(disability.placebo, disability.txa, disability.placebo.sims, disability.txa.sims)[[2]]



##  TRACE CALCULATIONS AND OUTCOMES  ## 


gen.trace <- function(clinical){

  # this unpacks each element of the clinical parameter to use in calculations below
  for(i in 1:length(clinical)) assign(names(clinical[i]),clinical[i])
  
  trace.names <- c("alive placebo", "dead placebo", "alive txa", "dead txa")
  matrix <- matrix(0, 13, 4)
  matrix[1,] <- c(1,0,1,0)
  colnames(matrix) <- trace.names

  # Calculate the annual risk of death, and make sure it changes at correct point 
  age.trace <- seq(from = age, to = age+1, by= 1/12)
  age.trace.floor <- floor(age.trace)
  a <- as.vector(acm %>% filter(age == unique(age.trace.floor)[[1]])) 
  b <- as.vector(acm %>% filter(age == unique(age.trace.floor)[[2]]))
  risk.year1 <- c( rep(a[[2]], sum(age.trace.floor==a[[1]])), rep(b[[2]], sum(age.trace.floor==b[[1]])))


# Calculate the first year Markov trace

for(t in 2:13){
  
  if(t==2){
  
  matrix[t,2] <- matrix[t-1,2] + (hi.risk + non.hi.risk)
  matrix[t,1] <- 1 - matrix[t,2]
  
  matrix[t,4] <- matrix[t-1,4] + (hi.risk * treatment.effect + non.hi.risk)
  matrix[t,3] <- 1 - matrix[t,4]
  
  } else {
    
  matrix[t,2] <- matrix[t-1,2] + matrix[t-1,1] * (1 - exp(-((risk.year1[t]/12)*smr.year1)))
  matrix[t,1] <- 1 - matrix[t,2]
  
  matrix[t,4] <- matrix[t-1,4] + matrix[t-1,3] * (1 - exp(-((risk.year1[t]/12)*smr.year1)))
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


# deterministic trace
trace.results <- gen.trace(clin.char)





## The costs in here are mostly named, but should really be changed to core components for the PSA part

gen.outcomes <- function(trace, util = utility, dec = utility.decrement, cost = costs, discount.c = disc.c, discount.o = disc.o){
  
  # discount
  
  d <- matrix( rep(1/((1+discount.c)^seq(0, time.horizon, 1)),2), time.horizon + 1, 2)
  o <- matrix( rep(1/((1+discount.o)^seq(0, time.horizon, 1)),2), time.horizon + 1, 2)
  
  # Costs
  
  # cost names here: cost.names
  
  cost.matrix <- matrix(0, time.horizon + 1, 2) # placebo / txa
  
  cost.matrix[1,] <- cost[2:3] + (cost[1] * c(0,1)) 
  cost.matrix[2,] <- trace[[2]][2,c(1,3)] * cost[4]
  cost.matrix[3:(time.horizon+1),] <- trace[[2]][3:(time.horizon+1),c(1,3)] * cost[5] 
  
  cost.matrix.d <- cost.matrix * d
  
  cost.sum <- apply(cost.matrix.d, 2, sum)
  

  # Utility

  utility.matrix <- matrix(0, time.horizon + 1, 2) # placebo / txa
  
  utility.matrix[2,1] <- mean(trace[[1]][2:13,1]) * (util[1] - dec[2,2])
  utility.matrix[2,2] <- mean(trace[[1]][2:13,3]) * (util[1] - dec[2,2])
  
  utility.matrix[3:(time.horizon+1),1] <- trace[[2]][3:(time.horizon+1),1] * (util[1] - dec[3:(time.horizon+1),2]) 
  utility.matrix[3:(time.horizon+1),2] <- trace[[2]][3:(time.horizon+1),3] * (util[2] - dec[3:(time.horizon+1),2])
  
  utility.matrix.d <- utility.matrix * o 
  
  utility.sum <- apply(utility.matrix.d, 2, sum) 

  
  ## ICER ## 
  
  icer <-  if((cost.sum[2] - cost.sum[1]) <= 0 & (utility.sum[2] - utility.sum[1]) > 0) "Intervention dominates" else 
      if((cost.sum[2] - cost.sum[1]) > 0 & (utility.sum[2] - utility.sum[1]) <= 0 ) "Control dominates" else
        (cost.sum[2] - cost.sum[1]) / (utility.sum[2] - utility.sum[1])  

    
  return(list(c(cost.placebo = cost.sum[1], 
                utility.placebo = utility.sum[1], 
                cost.txa = cost.sum[2], 
                utility.txa = utility.sum[2]), 
              icer = icer)) 

  
}

outcomes <- gen.outcomes(trace.results)


## PSA ## 

# Generate matrix to store results 
psa.results <- matrix(0, sims, 4)
colnames(psa.results) <- c("cost.placebo", "utility.placebo","cost.txa","utility.txa")



for(p in 1:sims){
  
  # Subset and assign existing sims
  clin.sim <- unlist(clin.char.sims[p,])
  utility.sim <- unlist(utility.sims[p,])
  cost.sim <- unlist(costs.sims[p,])
  
  trace.results.sim <- gen.trace(clin.sim)
  psa.results[p,] <- gen.outcomes(trace.results.sim, util = utility.sim, cost = cost.sim)[[1]] 
  
}

inc.costs <- psa.results[,3] - psa.results[,1]
inc.qaly <- psa.results[,4] - psa.results[,2]
icer <- inc.costs / inc.qaly
  
head(psa.results)
# Generate CEAC table

gen.ceac.table <- function(results, lambda.inc = 500){
  
  lambda <- seq(from = 0, to = 40000, by = lambda.inc)
    
  ## CEAC ## 
  
  lambda.table <- matrix(lambda, ncol = length(lambda), nrow = dim(results)[1], byrow = TRUE) 
  inmb.count <- ((results[,4] * lambda.table) - results[,3]) - ((results[,2] * lambda.table) - results[,1]) > 0 
  
  prob.ce <- apply(inmb.count, 2, mean) 
  ceac.table <- data.frame(lambda, prob.ce) 
  
  ## EVPI ## 
  
  evpi.table <- matrix(lambda, ncol = length(lambda), nrow = dim(results)[1], byrow = TRUE) 

  nmb.p <- ((results[,2] * evpi.table) - results[,1])  
  nmb.t <- ((results[,4] * evpi.table) - results[,3]) 
    
  av.nmb.p <- apply(nmb.p, 2, mean)
  av.nmb.t <- apply(nmb.t, 2, mean)

  evpi.mat <- matrix(0, dim(results)[1], length(av.nmb.p))
  
  for(i in 1:ncol(nmb.p)){
      if(av.nmb.p[i] >= av.nmb.t[i]) evpi.mat[,i] <- nmb.t[,i] - nmb.p[,i]
      if(av.nmb.t[i] >= av.nmb.p[i]) evpi.mat[,i] <- nmb.p[,i] - nmb.t[,i]
  }
  
  evpi.mat[evpi.mat<0] <- 0
  evpi.m <- apply(evpi.mat, 2, mean)
  
  evpi <- data.frame(lambda, evpi.m)

  return(list(ceac.table,
              evpi))
  
}



ceac <- gen.ceac.table(psa.results, 100)[[1]]
evpi <- gen.ceac.table(psa.results, 100)[[2]]


# Graphics  - TBC (take from other sources)




# CEAC # 

gen.ceac.graph = function(psa, save = FALSE) {
  
  z = ggplot(psa) + geom_line(aes(x=lambda, y=prob.ce), size=0.6) + 
    labs(x = "Willingness to pay (£)", text = element_text(size=4)) + 
    labs (y = "Probability cost-effective", text = element_text(size=4)) + theme_classic() +
    theme(legend.title = element_blank(), axis.title=element_text(face="bold"), 
          axis.title.x = element_text(margin = margin(t = 7, r = 0, b = 3, l = 0)), 
          axis.title.y = element_text(margin = margin(t = 0, r = 7, b = 0, l = 3)), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          legend.key.width=unit(1.8,"line"), text = element_text(size=7),
          plot.margin=unit(c(1.2,0.5,0,1.2),"cm")) + 
    scale_x_continuous(labels = scales::comma, breaks = c(seq(0,100000,5000)), limits = c(0,40000), expand = c(0, 0.1)) + 
    scale_y_continuous(limits = c(0,1), breaks=seq(0,1,0.1), expand = c(0, 0)) + 
    geom_vline(xintercept = 20000, linetype="dotted", size=0.25)
  
  
  if(save == TRUE) ggsave(paste("figures\\CEAC",Sys.Date(),".png"), z, width=107, height=70, dpi=300, units='mm')
  
  
  return(z)
  
}

gen.evpi.graph = function(evpi, save = FALSE) {
  
  z = ggplot(evpi) + geom_line(aes(x=lambda, y=evpi.m), size=0.6) + 
    labs(x = "Willingness to pay (£)", text = element_text(size=4)) + 
    labs (y = "EVPI (£)", text = element_text(size=4)) + theme_classic() +
    theme(legend.title = element_blank(), axis.title=element_text(face="bold"), 
          axis.title.x = element_text(margin = margin(t = 7, r = 0, b = 3, l = 0)), 
          axis.title.y = element_text(margin = margin(t = 0, r = 7, b = 0, l = 3)), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          legend.key.width=unit(1.8,"line"), text = element_text(size=7),
          plot.margin=unit(c(1.2,0.5,0,1.2),"cm")) + 
    scale_x_continuous(labels = scales::comma, breaks = c(seq(0,100000,5000)), limits = c(0,40000), expand = c(0, 0.1)) + 
    scale_y_continuous(expand = c(0, 0)) + 
    geom_vline(xintercept = 20000, linetype="dotted", size=0.25)
  
  if(save == TRUE) ggsave(paste("figures\\EVPI",Sys.Date(),".png"), z, width=107, height=70, dpi=300, units='mm')
  
  return(z)
  
}



gen.ceac.graph(ceac)
gen.evpi.graph(evpi)




##------------------------------##
##           EVPPI              ## 
##------------------------------##

# inner.loops <- 300
# outer.loops <- 300

# Select lambda values to be considered 
lambda <- seq(from = 0, to = 30000, by = 500)



# Sample all probabilistic parameters
clin.char.sims <- gen.clinical.characteristics()[[2]]
disability.placebo.sims <- gen.clinical.characteristics()[[5]]
disability.txa.sims <- gen.clinical.characteristics()[[6]]
utility.sims <- gen.utility.sims()
costs.sims <- gen.costs(disability.placebo, disability.txa, disability.placebo.sims, disability.txa.sims)[[2]]




gen.nmb <- function(results, lam = lambda){
  

  nmb.table <- matrix(c(lam), ncol = length(lam), nrow = dim(results)[1],  byrow = TRUE) 
  
  p <- ((results[,2] * nmb.table) - results[,1])  
  t <- ((results[,4] * nmb.table) - results[,3])

  nmb.p <- apply(p, 2, mean)
  nmb.t <- apply(t, 2, mean) 
    
  #colnames(nmb.p) <- as.character(lam)
  #colnames(nmb.t) <- as.character(lam)
  
  
  return(list(nmb.t, nmb.p))
  
}


# Generate matrices for EVPPI results to be stored in

inner.results <- matrix(0, inner.loops, 4)

evppi.results.placebo <- matrix(0, ncol = length(lambda), nrow = outer.loops)
colnames(evppi.results.placebo) <- as.character(lambda)
evppi.results.txa <- evppi.results.placebo

# EVPPI functions

gen.evppi.results <- function(evppi.results1 = evppi.results.placebo, evppi.results2 = evppi.results.txa, lam = lambda){
  
  ## calculate the mean NMB for placebo and txa, at each lambda 
  current.info1 <- apply(evppi.results1, 2, mean)
  current.info2 <- apply(evppi.results2, 2, mean)

  current.info <- pmax(current.info1, current.info2)

  evppi.array <- array(0, dim = c(outer.loops, length(lam), 2))
  evppi.array[,,1] <- evppi.results1
  evppi.array[,,2] <- evppi.results2

  perf.info.sims <- apply(evppi.array, c(1,2), max)
  perf.info <- apply(perf.info.sims, 2, mean)

  evppi.results <- c(perf.info - current.info)

  evppi <- data.frame(lambda, evppi.results)

  return(evppi)

}






## EVPPI Loops - 'Double Monte Carlo loop method' 



## EVPPI loops - Head injury and TXA treatment effect

for(a in 1:outer.loops){

## 1. Select the 'partial' parameter from the outer loop 
clin.sim <- unlist(clin.char.sims[a,])

  for(b in 1:inner.loops){
  
  # Select traditional parameters, minus the outer loop parameter
  
  #clin.sim <- unlist(clin.char.sims[b,])
  clin.sim[4:5] <- unlist(clin.char.sims[b,4:5]) # Keep SMRs in PSA 
  utility.sim <- unlist(utility.sims[b,])
  cost.sim <- unlist(costs.sims[b,])
  
  trace.results.sim <- gen.trace(clin.sim)
  inner.results[b,] <- gen.outcomes(trace.results.sim, util = utility.sim, cost = cost.sim)[[1]] 
  }

  #after each inner loop PSA, calculate the mean NMB for each tx and store the results
  nmb <- gen.nmb(inner.results)
  evppi.results.placebo[a,] <- nmb[[1]]
  evppi.results.txa[a,] <- nmb[[2]]
  
}

# Calculate the EVPPI 
evppi.head.injury <- gen.evppi.results(evppi.results.placebo, evppi.results.txa, lambda)




## EVPPI loops - Head injury and TXA treatment effect

for(a in 1:outer.loops){
  
  ## 1. Select the 'partial' parameter from the outer loop 
  clin.sim <- unlist(clin.char.sims[a,])
  
  for(b in 1:inner.loops){
    
    # Select traditional parameters, minus the outer loop parameter
    
    #clin.sim <- unlist(clin.char.sims[b,])
    clin.sim[1:3] <- unlist(clin.char.sims[b,1:3]) # only SMRs excl from PSA 
    utility.sim <- unlist(utility.sims[b,])
    cost.sim <- unlist(costs.sims[b,])
    
    trace.results.sim <- gen.trace(clin.sim)
    inner.results[b,] <- gen.outcomes(trace.results.sim, util = utility.sim, cost = cost.sim)[[1]] 
  }
  
  #after each inner loop PSA, calculate the mean NMB for each tx and store the results
  nmb <- gen.nmb(inner.results)
  evppi.results.placebo[a,] <- nmb[[1]]
  evppi.results.txa[a,] <- nmb[[2]]
  
}

# Calculate the EVPPI 
evppi.smr <- gen.evppi.results(evppi.results.placebo, evppi.results.txa, lambda)



## EVPPI loops - Utility  

for(a in 1:outer.loops){
  
  ## 1. Select the 'partial' parameter from the outer loop 
  utility.sim <- unlist(utility.sims[a,])
  
  for(b in 1:inner.loops){
    
    # Select traditional parameters, minus the outer loop parameter
    
    clin.sim <- unlist(clin.char.sims[b,])
    #utility.sim <- unlist(utility.sims[b,])
    cost.sim <- unlist(costs.sims[b,])
    
    trace.results.sim <- gen.trace(clin.sim)
    inner.results[b,] <- gen.outcomes(trace.results.sim, util = utility.sim, cost = cost.sim)[[1]] 
  }
  
  #after each inner loop PSA, calculate the mean NMB for each tx and store the results
  nmb <- gen.nmb(inner.results)
  evppi.results.placebo[a,] <- nmb[[1]]
  evppi.results.txa[a,] <- nmb[[2]]
  
}

# Calculate the EVPPI 
evppi.utility <- gen.evppi.results(evppi.results.placebo, evppi.results.txa, lambda)




## EVPPI loops - Costs  

for(a in 1:outer.loops){
  
  ## 1. Select the 'partial' parameter from the outer loop 
  cost.sim <- unlist(costs.sims[a,])
  
  for(b in 1:inner.loops){
    
    # Select traditional parameters, minus the outer loop parameter
    
    clin.sim <- unlist(clin.char.sims[b,])
    utility.sim <- unlist(utility.sims[b,])
    #cost.sim <- unlist(costs.sims[b,])
    
    trace.results.sim <- gen.trace(clin.sim)
    inner.results[b,] <- gen.outcomes(trace.results.sim, util = utility.sim, cost = cost.sim)[[1]] 
  }
  
  #after each inner loop PSA, calculate the mean NMB for each tx and store the results
  nmb <- gen.nmb(inner.results)
  evppi.results.placebo[a,] <- nmb[[1]]
  evppi.results.txa[a,] <- nmb[[2]]
  
}

# Calculate the EVPPI 
evppi.costs <- gen.evppi.results(evppi.results.placebo, evppi.results.txa, lambda)


## Reshaping? 

evppi.wide <- data.frame(evppi.head.injury,
                                 evppi.smr[,2],
                                 evppi.utility[,2],
                                 evppi.costs[,2])

colnames(evppi.wide) <- c('lambda', 'death following head injury and treatment effect', 'SMR', 'utility', 'costs')

evppi.long <- evppi.wide %>% gather(Parameters, VoI, 2:5)
#evppi.long.pop <- reshape2::melt(evppi.wide.pop, id.vars = c("lambda"))

# Plots 

gen.evppi.graph = function(evppi, save = FALSE) {
  
  z = ggplot(evppi) + geom_line(aes(x=lambda, y=VoI, colour = Parameters), size=0.6) + 
    labs(x = "Willingness to pay (£)", text = element_text(size=4)) + 
    labs (y = "EVPPI (£)", text = element_text(size=4)) + theme_classic() +
    theme(legend.title = element_blank(), axis.title=element_text(face="bold"), 
          axis.title.x = element_text(margin = margin(t = 7, r = 0, b = 3, l = 0)), 
          axis.title.y = element_text(margin = margin(t = 0, r = 7, b = 0, l = 3)), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          legend.key.width=unit(1.8,"line"), text = element_text(size=7),
          plot.margin=unit(c(1.2,0.5,0,1.2),"cm")) + 
    scale_x_continuous(labels = scales::comma, breaks = c(seq(0,100000,5000)), limits = c(0,40000), expand = c(0, 0.1)) + 
    scale_y_continuous(expand = c(0, 0)) + 
    geom_vline(xintercept = 20000, linetype="dotted", size=0.25)
  
  
  
  if(save == TRUE) ggsave(paste("figures\\EVPPI",Sys.Date(),".png"), z, width=107, height=70, dpi=300, units='mm')
  
  return(z)  
  
}

gen.evppi.graph(evppi.long)


