
## CRASH-Mild - Economic model ## 



if(!require(dplyr)) install.packages('dplyr')
library(dplyr)



# Options

time.horizon = 60
disc.c <- 0.035
disc.o <- 0.035

sims <- 10

age <- 57.73704



## Clinical parameters

treatment.effect <- 0.6672781
treatment.effect.sims <- exp(rnorm(sims, log(treatment.effect), 0.2472015))

hi.risk <- 13/440
hi.risk.sims <- rbeta(sims, 13, 440-13)

non.hi.risk <- 9/920
non.hi.risk.sims <- rbeta(sims, 9, 920-9)

smr.year1 <- 4
smr.year1.sims <- rnorm(sims, smr.year1, 0.41632922)

smr.year2 <- 30.99/13.72
smr.year2.sims <- rnorn(sims, smr.year2, 0.2350955) 




## Generate long-term risk of death

gen.acm <- function(age = age_start, male = 0.5){
  
  lifetable <- read.csv("inputs/ONS_life_tables_2017-19.csv", header=TRUE)
  lifetable <- as.data.frame(lifetable[c(1,2,8)])
  av <- lifetable
  av[,2] <- lifetable[,2] * male + lifetable[,3] * (1 - male)
  colnames(av) <- c("age", "rate")

  return(av[,1:2])
  
}

acm <- gen.acm()




## HRQoL

gen.utility <- function(){

utility.full <- 1
utility.good <- 0.894
utility.moderate <- 0.675
utility.severe <- 0.382
utility.vegetative <- -0.178

utility.values <- c(utility.full, utility.good, utility.moderate, utility.severe, utility.vegetative)

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

utility.placebo <- sum(utility.values*disability.placebo)/sum(disability.placebo)
utility.txa <- sum(utility.values*disability.txa)/sum(disability.txa)

return(c(placebo = utility.placebo, 
       txa = utility.txa))

}

gen.utility()



gen.utility.sims <- function(){
  
  utility.full <- rep(1, sims)
  utility.good <- rbeta(sims, 49.9585894, 5.9235016)
  utility.moderate <- rbeta(sims, 30.538, 14.7034815)
  utility.severe <- rbeta(sims, 10.9303087,	17.6830648)
  utility.vegetative <- rnorm(sims, -0.178, 0.19) ## needs normal distribution
  
  utility.values <- data.frame(utility.full, utility.good, utility.moderate, utility.severe, utility.vegetative)
  
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
  
  
  util.plac <- utility.values * matrix(disability.placebo, sims, length(disability.placebo), byrow= T)
  util.values.plac <- apply(util.plac,  1,  function(x) sum(x) / sum(disability.placebo))
  
  util.txa <- utility.values * matrix(disability.placebo, sims, length(disability.placebo), byrow= T)
  util.values.txa <- apply(util.txa,  1,  function(x) sum(x) / sum(disability.placebo))
  
  
  return(data.frame(placebo = util.values.plac, 
                    txa = util.values.txa))
  
}

gen.utility.sims()


## need to add utility decrements 





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







## Year 1 trace


matrix <- matrix(0, 366, 4)
matrix[1,] <- c(1,0,1,0)
age.trace <- seq(from = age, to = age+1, by= 1/365)


# Generate year 1 risk of death 

age.trace.floor <- floor(age_trace)
a <- as.vector(acm %>% filter(age == unique(age.trace.floor)[[1]])) 
b <- as.vector(acm %>% filter(age == unique(age.trace.floor)[[2]]))
risk.year1 <- c( rep(a[[2]], sum(age.trace.floor==a[[1]])), rep(b[[2]], sum(age.trace.floor==b[[1]])))

## 28 days  + 1 as first day is day 0

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

first.year.trace.placebo <- as.data.frame(matrix[,1:2])
first.year.trace.txa <- as.data.frame(matrix[,3:4])



# Clinical trace 

trace.matrix <- matrix(0, time.horizon + 1, 4) 
dim(trace.matrix)





