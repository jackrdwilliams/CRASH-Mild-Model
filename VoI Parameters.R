## VoI effective population estimate ## 


ons.pop <- read.csv("inputs/ONS_pop.csv", header=TRUE)
ons.pop.proj <- read.csv("inputs/ons.pop.proj.csv", header = TRUE)
colnames(ons.pop) <- c("age", "male", "female")

# Settings and inputs for analysis
age.range <- c(70:90)
age.range.alt <- c(60:90)
discount <- 0.035
voi.th <- 20
prop.mild <- 0.9 

## Linear approximation and incidence tables (Global burden of disease)

# data.m <- read.csv("inputs/age_standardised_incidence_male.csv")
# data.f <- read.csv("inputs/age_standardised_incidence_female.csv")
# 
# incidence.ratio <- 260/369
# est.male   <-  approx(data.m$Age, data.m$Male/100000, xout=c(age.range), method="linear")
# est.female <-  approx(data.f$Age, data.f$Female/100000, xout=c(age.range), method="linear")
# incidence <- data.frame(age = est.male[[1]], male = est.male[[2]], female = est.female[[2]]) * incidence.ratio 

## Alternative incidence - Tennent 2005
incidence <- 410.8/100000

# Population (fixed per year)
population <- sum(subset(ons.pop, age >= age.range[1])) 

# Population projection (ONS 2018 projections)
# Estimate from ONS projections
rownames(ons.pop.proj) <- ons.pop.proj[,1]
ons.pop.proj[,1] <- seq(0, 100, 5)
pop <- ons.pop.proj[,-c(2,3,4,5)] ## Removing 2018 to 2021
pop.sums <- apply(subset(pop, Ages>=70), 2, sum) 
population.vec <- as.numeric(pop.sums[2:(voi.th+1)] * 1000) # adjusting for per 1000 data


## Multiply UK age-adjusted incidence by population 
total <- population.vec * incidence * prop.mild

effective.population <- sum( total * (1/(1+discount)^(0:(voi.th-1))))
effective.population

# 60 + population ('alt')

population.alt <- sum(subset(ons.pop, age >= age.range.alt[1])) 
total.alt <- population.alt * incidence * prop.mild
effective.population.alt <- sum( total.alt * (1/(1+discount)^(0:(voi.th-1))))



