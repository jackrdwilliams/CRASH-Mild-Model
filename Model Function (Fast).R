## Alternative run model (faster)


run.model2 <- function(clinical = clin.char, dis.plac = disability.placebo, dis.txa = disability.txa, util.values = utility, cost = costs, 
                      dec = utility.decrement, discount.c = disc.c, discount.o = disc.o) {
  
  
  for(i in 1:length(clinical)) assign(names(clinical[i]),clinical[i])
  
  matrix <- matrix(0, 13, 4)
  matrix[1,] <- c(1,0,1,0)

  age.trace <- seq(from = age, to = age+1, by= 1/12)
  age.trace.floor <- floor(age.trace)
  a <- as.vector(acm %>% filter(age == unique(age.trace.floor)[[1]])) 
  b <- as.vector(acm %>% filter(age == unique(age.trace.floor)[[2]]))
  risk.year1 <- c(rep(a[[2]], sum(age.trace.floor==a[[1]])), rep(b[[2]], sum(age.trace.floor==b[[1]])))
  
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
  
  trace.matrix <- matrix(0, time.horizon + 1, 4) 
  trace.matrix[1,] <- c(1,0,1,0)
  trace.matrix[2,] <- as.numeric(tail(trace.y1, 1))
  
  for(i in 3:(time.horizon + 1)){
    trace.matrix[i,2] <- trace.matrix[i-1,2] + trace.matrix[i-1,1] * (1 - exp(-(acm[i,2] * smr.year2)))
    trace.matrix[i,1] <- 1 - trace.matrix[i,2] 
    trace.matrix[i,4] <- trace.matrix[i-1,4] + trace.matrix[i-1,3] * (1 - exp(-(acm[i,2] * smr.year2)))
    trace.matrix[i,3] <- 1 - trace.matrix[i,4] 
  }
  
  trace <- list(trace.y1, trace.matrix)

  

  d <- matrix( rep(1/((1+discount.c)^(0:time.horizon)),2), time.horizon + 1, 2)
  o <- matrix( rep(1/((1+discount.o)^(0:time.horizon)),2), time.horizon + 1, 2)

  
  st.mon.costs <- c(cost[4], cost[5], cost[6])
  lt.mon.costs <- c(cost[7], cost[8], cost[9])
  st.mon.plac <- sum( st.mon.costs * c(dis.plac[1] + dis.plac[2], dis.plac[3], dis.plac[4] + dis.plac[5]))
  lt.mon.plac <- sum( lt.mon.costs * c(dis.plac[1] + dis.plac[2], dis.plac[3], dis.plac[4] + dis.plac[5]))
  st.mon.txa <- sum( st.mon.costs * c(dis.txa[1] + dis.txa[2], dis.txa[3], dis.txa[4] + dis.txa[5]))
  lt.mon.txa <- sum( lt.mon.costs * c(dis.txa[1] + dis.txa[2], dis.txa[3], dis.txa[4] + dis.txa[5]))
  
  cost.matrix <- matrix(0, time.horizon + 1, 2) 
  cost.matrix[1,] <- cost[2:3] + (cost[1] * c(0,1)) 
  cost.matrix[2,] <- trace[[2]][2,c(1,3)] * c(st.mon.plac, st.mon.txa)
  cost.matrix[3:(time.horizon+1),] <- trace[[2]][3:(time.horizon+1),c(1,3)] * c(lt.mon.plac, lt.mon.txa)
  
  cost.matrix.d <- cost.matrix * d
  
  cost.sum <- apply(cost.matrix.d, 2, sum)
  
  
  util.values.plac <- sum(util.values * dis.plac)
  util.values.txa <- sum(util.values * dis.txa)
  
  utility.matrix <- matrix(0, time.horizon + 1, 2) 
  utility.matrix[2,1] <- mean(trace[[1]][2:13,1]) * (util.values.plac - dec[2,2])
  utility.matrix[2,2] <- mean(trace[[1]][2:13,3]) * (util.values.txa - dec[2,2])
  utility.matrix[3:(time.horizon+1),1] <- trace[[2]][3:(time.horizon+1),1] * (util.values.plac - dec[3:(time.horizon+1),2]) 
  utility.matrix[3:(time.horizon+1),2] <- trace[[2]][3:(time.horizon+1),3] * (util.values.txa - dec[3:(time.horizon+1),2])
  utility.matrix.d <- utility.matrix * o 
  
  utility.sum <- apply(utility.matrix.d, 2, sum) 


  return(c(cost.placebo = cost.sum[1], 
           utility.placebo = utility.sum[1], 
           cost.txa = cost.sum[2], 
           utility.txa = utility.sum[2])) 
  
}

sims <- 10000

time1 <- Sys.time()
for(p in 1:sims){
  
  # Subset and assign existing sims
  clin.sim <- unlist(clin.char.sims[p,])
  dis.placebo.sim <- unlist(disability.placebo.sims[p,])
  dis.txa.sim <- unlist(disability.txa.sims[p,])
  utility.sim <- unlist(utility.sims[p,])
  cost.sim <- unlist(costs.sims[p,])
  
  psa.results[p,] <- run.model(clin.sim, dis.placebo.sim, dis.txa.sim, utility.sim, cost.sim)[[1]] 
}

time2 <- Sys.time()
for(p in 1:sims){
  
  # Subset and assign existing sims
  clin.sim <- unlist(clin.char.sims[p,])
  dis.placebo.sim <- unlist(disability.placebo.sims[p,])
  dis.txa.sim <- unlist(disability.txa.sims[p,])
  utility.sim <- unlist(utility.sims[p,])
  cost.sim <- unlist(costs.sims[p,])
  
  psa.results[p,] <- run.model2(clin.sim, dis.placebo.sim, dis.txa.sim, utility.sim, cost.sim)[[1]]
}
time3 <- Sys.time()

time2 - time1
time3 - time2