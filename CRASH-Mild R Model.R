
## CRASH-Mild - Economic model ## 



if(!require(dplyr)) install.packages('dplyr')
library(dplyr)



treatment.effect <- 0.667

treatment.effect.lci <- 0.411
treatment.effect.uci <- 1.083

hi.risk <- 9/920
non.hi.risk <- 13/440

age <- 57.73704

smr.year1 <- 4
smr.year2 <- 2






gen_acm <- function(age = age_start, male = 0.5){
  
  lifetable <- read.csv("inputs/ONS_life_tables_2017-19.csv", header=TRUE)
  lifetable <- as.data.frame(lifetable[c(1,2,8)])
  av <- lifetable
  av[,2] <- lifetable[,2] * male + lifetable[,3] * (1 - male)
  colnames(av) <- c("age", "rate")

  return(av[,1:2])
  
}

acm <- gen_acm()



matrix <- matrix(0, 366, 2)
matrix[1,] <- c(1,0)
age_trace <- seq(from = age, to = age+1, by= 1/365)

# Generate year 1 risk of death 

age_trace_floor <- floor(age_trace)
a <- as.vector(acm %>% filter(age == unique(age_trace_floor)[[1]])) 
b <- as.vector(acm %>% filter(age == unique(age_trace_floor)[[2]]))
risk_year1 <- c( rep(a[[2]], sum(age_trace_floor==a[[1]])), rep(b[[2]], sum(age_trace_floor==b[[1]])))

## Year 1 trace

## 28 days  + 1 as first day is day 0

for(t in 2:366){
  
  if(t<=28+1){
  
  matrix[t,2] <- matrix[t-1,2] + ((hi.risk + non.hi.risk)/28)
  matrix[t,1] <- 1 - matrix[t,2]
  
  } else {
    
  matrix[t,2] <- matrix[t-1,2] + matrix[t-1,1] * (1 - exp(-(risk_year1[t]/365)))
  matrix[t,1] <- 1 - matrix[t,2]
  }
  
    
}

first.year.trace <- as.data.frame(matrix)


