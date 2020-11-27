
## CRASH-Mild - Economic model ## 



if(!require(dplyr)) install.packages('dplyr')
library(dplyr)




## Clinical parameters

treatment.effect <- 0.6672781

treatment.effect.lci <- 0.411
treatment.effect.uci <- 1.083

hi.risk <- 13/440
non.hi.risk <- 9/920

age <- 57.73704

smr.year1 <- 4
smr.year2 <- 30.99/13.72




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



matrix <- matrix(0, 366, 4)
matrix[1,] <- c(1,0,1,0)
age.trace <- seq(from = age, to = age+1, by= 1/365)




## HRQoL

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


## need to add utility decrements 





## Costs 

# Treamtent costs 
cost.txa.dose <- 6
cost.sodium <- 0.55 + 2.7 # 55p for 100ml, 270 for 500ml 
cost.needle <- 0.05 
cost.nurse <- 12.95

cost.treatment <- sum(cost.txa.dose, cost.sodium, cost.needle, cost.nurse)

# hospital costs

los.placebo <-  12.44828
los.txa <- 12.44828

hospital.cost.initial <- 455.4503
hospital.cost.day <- 313.8758

hospital.cost.placebo <- hospital.cost.initial + hospital.cost.day * los.placebo
hospital.cost.txa <- hospital.cost.initial + hospital.cost.day * los.txa

# monitoring costs 

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

tail(matrix)






