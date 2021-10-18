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

## Alternative incidence - Tennent 2005 (75 and over)
incidence <- 410.8/100000


# Population (fixed per year)
# population <- sum(subset(ons.pop, age >= age.range[1])) 

# Population projection (ONS 2018 projections)
# Estimate from ONS projections
rownames(ons.pop.proj) <- ons.pop.proj[,1]
ons.pop.proj[,1] <- seq(0, 100, 5)
pop <- ons.pop.proj[,c(1, 6:(5+voi.th))] ## Removing 2018 to 2021, and years beyond TH
pop.matrix <- subset(pop, Ages>=70) 
pop.matrix.inc <- pop.matrix[,-1] * 1000 * incidence # adjusting for per 1000 data and incidence 

population.vec <- colSums(pop.matrix.inc)


## Multiply UK age-adjusted incidence by population 
total <- as.numeric(population.vec) * prop.mild 

effective.population <- sum(total * (1/(1+discount)^(0:(voi.th-1))))
effective.population


# 60 + population ('alt')
incidence.alt <- incidence * 0.5
pop.matrix.alt <- subset(pop, Ages>=60 & Ages<70)
pop.matrix.inc.alt <- pop.matrix.alt[,-1] * 1000 * incidence.alt # adjusting for per 1000 data and incidence 
population.vec.alt <- colSums(pop.matrix.inc.alt)

total.alt <- as.numeric(population.vec.alt) * prop.mild  
effective.population.alt <- effective.population + sum( total.alt * (1/(1+discount)^(0:(voi.th-1))))


## shorter/longer TH

#rownames(ons.pop.proj) <- ons.pop.proj[,1]
#ons.pop.proj[,1] <- seq(0, 100, 5)
pop30 <- ons.pop.proj[,c(1, 6:(5+30))] ## Removing 2018 to 2021, and years beyond TH
pop.matrix30 <- subset(pop30, Ages>=70) 
pop.matrix.inc30 <- colSums(pop.matrix30[,-1] * 1000 * incidence) # adjusting for per 1000 data and incidence 
total30 <- as.numeric(pop.matrix.inc30) * prop.mild 
effective.population30 <- sum(total30 * (1/(1+discount)^(0:(30-1))))
effective.population5 <- sum(total30[1:5] * (1/(1+discount)^(0:(5-1))))
effective.population10 <- sum(total30[1:10] * (1/(1+discount)^(0:(10-1))))



