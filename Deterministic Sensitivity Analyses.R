## Deterministic analyses ## 

#source("CRASH-Mild R Model.R")

## DSA ## 


## Head injury risk x Treatment effect threshold 

clin.char.dsa <- clin.char
tx.effect.goalseek <- function(x){
  clin.char.dsa[1] <- x 
  z <- run.model(clin.char.dsa, output.type = 2) - 20000
  return(z)
}
hi.risk.thresholds <- seq(from = 0.0001, to = 0.005, by = 0.0001)
threshold.mat <- matrix(0, ncol = 4, nrow = length(hi.risk.thresholds))
threshold.mat[,1] <- hi.risk.thresholds
colnames(threshold.mat) <- c("Risk", "70","80","90")

for(a in 1:3){
  age.v <- c(70,80,90)
  age <- age.v[a]
  time.horizon = min(60, 100-ceiling(age))
  acm <- gen.acm()
  utility.decrement <- gen.utility.dec()
  
    for(i in 1:length(hi.risk.thresholds)){
      clin.char.dsa[2] <- hi.risk.thresholds[i] 
      res <- uniroot(tx.effect.goalseek, c(0.0001,0.9999999))
      threshold.mat[i,a+1] <- res$root
  }
} 

# Age based parameters

threshold.long <- as.data.frame(threshold.mat)  %>% gather(Age, Threshold, 2:4)
ggplot(threshold.long) + geom_line( aes(x = Risk, y = Threshold, group = Age))






age <- 90
acm <- gen.acm()
utility.decrement <- gen.utility.dec()
clin.char.dsa
run.model(clinical = clin.char.dsa)

