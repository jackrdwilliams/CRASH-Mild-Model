## VoI effective population estimate ## 

data.m <- read.csv("inputs/age_standardised_incidence_male.csv")
data.f <- read.csv("inputs/age_standardised_incidence_female.csv")
ons.pop <- read.csv("inputs/ONS_pop.csv", header=TRUE)
colnames(ons.pop) <- c("age", "male", "female")

# Settings and inputs for analysis
age.range <- c(70:90)
incidence.ratio <- 260/369
discount <- 0.035
voi.th <- 10

# Linear approximation and incidence tables
est.male   <-  approx(data.m$Age, data.m$Male/100000, xout=c(age.range), method="linear")
est.female <-  approx(data.f$Age, data.f$Female/100000, xout=c(age.range), method="linear")
incidence <- data.frame(age = est.male[[1]], male = est.male[[2]], female = est.female[[2]])


## Multiply table with ratio of UK to global incidence (assumes UK incidence follows same trend)
subset(ons.pop, age >= age.range[1]) 
population <- subset(ons.pop, age >= age.range[1]) 
age.counts <- population[,2:3] * incidence[,2:3] * incidence.ratio  
total <- sum(apply(age.counts, 2, sum))

## Multiply UK age-adjusted incidence by population 
effective.population <- sum( total * (1/(1+discount)^(0:(voi.th-1))))
effective.population

