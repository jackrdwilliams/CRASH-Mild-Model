##------------------------------##
##           EVPPI              ## 
##------------------------------##

source("CRASH-Mild R Model.R")

# inner.loops <- 300
# outer.loops <- 300

# Sample all probabilistic parameters

# clin.char.sims <- gen.clinical.characteristics()[[2]]
# disability.placebo.sims <- gen.clinical.characteristics()[[5]]
# disability.txa.sims <- disability.placebo.sims # same as placebo (equal for both arms)
# #disability.txa.sims <- gen.clinical.characteristics()[[6]]
# utility.sims <- gen.utility.sims()
# costs.sims <- gen.costs(disability.placebo, disability.txa, disability.placebo.sims, disability.txa.sims)[[2]]

evppi.start.time <- Sys.time()

# Generate matrices for EVPPI results to be stored in

inner.results <- matrix(0, inner.loops, 4)

evppi.results.placebo <- matrix(0, ncol = length(lambda), nrow = outer.loops)
colnames(evppi.results.placebo) <- as.character(lambda)
evppi.results.txa <- evppi.results.placebo


# EVPPI functions

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
  
  evppi <- data.frame(lam, evppi.results)
  
  return(evppi)
  
}



## EVPPI Loops - 'Double Monte Carlo loop method' 

## EVPPI loops - Head injury and TXA treatment effect
for(a in 1:outer.loops){
  
  ## 1. Select the 'partial' parameter from the outer loop 
  
  clin.sim <- unlist(clin.char.sims[a,])
  
  for(b in 1:inner.loops){
    
    # Select traditional parameters, minus the outer loop parameter
    
    
    clin.sim[4:5] <- unlist(clin.char.sims[b,4:5]) # Keep SMRs in PSA 
    dis.placebo.sim <- unlist(disability.placebo.sims[b,])
    dis.txa.sim <- unlist(disability.txa.sims[b,])
    utility.sim <- unlist(utility.sims[b,])
    cost.sim <- unlist(costs.sims[b,])
    
    inner.results[b,] <- run.model(clin.sim, dis.placebo.sim, dis.txa.sim, utility.sim, cost.sim)[[1]] 
  }
  
  #after each inner loop PSA, calculate the mean NMB for each tx and store the results
  nmb <- gen.nmb(inner.results)
  evppi.results.placebo[a,] <- nmb[[1]]
  evppi.results.txa[a,] <- nmb[[2]]
  
}
evppi.head.injury <- gen.evppi.results(evppi.results.placebo, evppi.results.txa, lambda)


## EVPPI loops - Head injury and TXA treatment effect
for(a in 1:outer.loops){
  
  ## 1. Select the 'partial' parameter from the outer loop 
  clin.sim <- unlist(clin.char.sims[a,])
  
  for(b in 1:inner.loops){
    
    # Select traditional parameters, minus the outer loop parameter
    
    #clin.sim <- unlist(clin.char.sims[b,])
    clin.sim[1:3] <- unlist(clin.char.sims[b,1:3]) # only SMRs excl from PSA 
    dis.placebo.sim <- unlist(disability.placebo.sims[b,])
    dis.txa.sim <- unlist(disability.txa.sims[b,])
    utility.sim <- unlist(utility.sims[b,])
    cost.sim <- unlist(costs.sims[b,])
    
    
    inner.results[b,] <- run.model(clin.sim, dis.placebo.sim, dis.txa.sim, utility.sim, cost.sim)[[1]] 
  }
  
  #after each inner loop PSA, calculate the mean NMB for each tx and store the results
  nmb <- gen.nmb(inner.results)
  evppi.results.placebo[a,] <- nmb[[1]]
  evppi.results.txa[a,] <- nmb[[2]]
  
}
evppi.smr <- gen.evppi.results(evppi.results.placebo, evppi.results.txa, lambda)


## EVPPI loops - Disability / outcomes   
for(a in 1:outer.loops){
  
  ## 1. Select the 'partial' parameter from the outer loop 
  dis.placebo.sim <- unlist(disability.placebo.sims[a,])
  dis.txa.sim <- unlist(disability.txa.sims[a,])
  
  for(b in 1:inner.loops){
    
    # Select traditional parameters, minus the outer loop parameter
    
    clin.sim <- unlist(clin.char.sims[b,])
    # dis.placebo.sim <- unlist(disability.placebo.sims[b,])
    # dis.txa.sim <- unlist(disability.txa.sims[b,])
    utility.sim <- unlist(utility.sims[b,])
    cost.sim <- unlist(costs.sims[b,])
    
    inner.results[b,] <- run.model(clin.sim, dis.placebo.sim, dis.txa.sim, utility.sim, cost.sim)[[1]] 
  }
  
  #after each inner loop PSA, calculate the mean NMB for each tx and store the results
  nmb <- gen.nmb(inner.results)
  evppi.results.placebo[a,] <- nmb[[1]]
  evppi.results.txa[a,] <- nmb[[2]]
  
}
evppi.disability <- gen.evppi.results(evppi.results.placebo, evppi.results.txa, lambda)


## EVPPI loops - Utility  
for(a in 1:outer.loops){
  
  ## 1. Select the 'partial' parameter from the outer loop 
  utility.sim <- unlist(utility.sims[a,])
  
  for(b in 1:inner.loops){
    
    # Select traditional parameters, minus the outer loop parameter
    
    clin.sim <- unlist(clin.char.sims[b,])
    dis.placebo.sim <- unlist(disability.placebo.sims[b,])
    dis.txa.sim <- unlist(disability.txa.sims[b,])
    #utility.sim <- unlist(utility.sims[b,])
    cost.sim <- unlist(costs.sims[b,])
    
    inner.results[b,] <- run.model(clin.sim, dis.placebo.sim, dis.txa.sim, utility.sim, cost.sim)[[1]] 
  }
  
  #after each inner loop PSA, calculate the mean NMB for each tx and store the results
  nmb <- gen.nmb(inner.results)
  evppi.results.placebo[a,] <- nmb[[1]]
  evppi.results.txa[a,] <- nmb[[2]]
  
}
evppi.utility <- gen.evppi.results(evppi.results.placebo, evppi.results.txa, lambda)


## EVPPI loops - Costs  
for(a in 1:outer.loops){
  
  ## 1. Select the 'partial' parameter from the outer loop 
  cost.sim <- unlist(costs.sims[a,])
  
  for(b in 1:inner.loops){
    
    # Select traditional parameters, minus the outer loop parameter
    
    clin.sim <- unlist(clin.char.sims[b,])
    dis.placebo.sim <- unlist(disability.placebo.sims[b,])
    dis.txa.sim <- unlist(disability.txa.sims[b,])
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
evppi.costs <- gen.evppi.results(evppi.results.placebo, evppi.results.txa, lambda)


## Reshaping and plotting ## 

evppi.wide <- data.frame(evppi.head.injury,
                         evppi.smr[,2],
                         evppi.disability[,2],
                         evppi.utility[,2],
                         evppi.costs[,2])

colnames(evppi.wide) <- c('lambda', 'death following head injury and treatment effect', 'SMR', 'disability', 'utility', 'costs')

evppi.long <- evppi.wide %>% gather(Parameters, VoI, 2:6)
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

gen.evpi.graph(evppi.head.injury)
gen.evppi.graph(evppi.long)




#### EVPPI for trial parameters 

for(a in 1:outer.loops){
  
  ## 1. Select trial parameters
  clin.sim[1:3] <- unlist(clin.char.sims[a,1:3]) # HI risk and tx effect
  dis.placebo.sim <- unlist(disability.placebo.sims[a,])
  dis.txa.sim <- unlist(disability.txa.sims[a,])  
  
  
  for(b in 1:inner.loops){
    
    # Select traditional parameters, minus the outer loop parameter
    clin.sim[4:5] <- unlist(clin.char.sims[b,4:5]) # Keep SMRs in PSA 
    utility.sim <- unlist(utility.sims[b,]) ## Utility values
    cost.sim[4:9] <- unlist(costs.sims[b,4:9])
    
    inner.results[b,] <- run.model(clin.sim, dis.placebo.sim, dis.txa.sim, utility.sim, cost.sim)[[1]] 
  }
  
  #after each inner loop PSA, calculate the mean NMB for each tx and store the results
  nmb <- gen.nmb(inner.results)
  evppi.results.placebo[a,] <- nmb[[1]]
  evppi.results.txa[a,] <- nmb[[2]]
  
}

evppi.trial.parms <- gen.evppi.results(evppi.results.placebo, evppi.results.txa, lambda)

## Reshape / plot

evppi.trial <- cbind(evpi, evppi.trial.parms[,2])

dim(evppi.trial)
colnames(evppi.trial) <- c('lambda', 'EVPI', 'EVPPI for trial parameters')
evppi.trial.long <- evppi.trial %>% gather(Parameters, VoI, 2:3)
gen.evppi.graph(evppi.trial.long)


gen.evppi.trial.graph = function(evppi, save = FALSE) {
  
  z = ggplot(evppi) + geom_line(aes(x=lambda, y=VoI, colour = Parameters, linetype = Parameters), size=0.6) + 
    labs(x = "Willingness to pay (£)", text = element_text(size=4)) + 
    labs(y = "EVPPI (£)", text = element_text(size=4)) + theme_classic() +
    theme(legend.title = element_blank(), axis.title=element_text(face="bold"), 
          axis.title.x = element_text(margin = margin(t = 7, r = 0, b = 3, l = 0)), 
          axis.title.y = element_text(margin = margin(t = 0, r = 7, b = 0, l = 3)), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          legend.key.width=unit(1.8,"line"), text = element_text(size=7),
          plot.margin=unit(c(1.2,0.5,0,1.2),"cm")) + 
    scale_x_continuous(labels = scales::comma, breaks = c(seq(0,100000,5000)), limits = c(0,40000), expand = c(0, 0.1)) + 
    scale_y_continuous(expand = c(0, 0)) + 
    scale_linetype_manual(values=c("solid", "longdash")) + 
    scale_color_manual(values=c("black", "black")) + 
    geom_vline(xintercept = 20000, linetype="dotted", size=0.25)
  
  
  if(save == TRUE) ggsave(paste("figures\\EVPPI-Trial",Sys.Date(),".png"), z, width=180, height=100, dpi=300, units='mm')
  
  return(z)  
  
}
gen.evppi.trial.graph(evppi.trial.long)

evppi.stop.time <- Sys.time()
evppi.stop.time - evppi.start.time

## Save EVPPI's

save(evppi.long, file=paste("stored results/evppi.",Sys.Date(),".Rda", sep=""))
save(evppi.trial.long, file=paste("stored results/evppi.trial.",Sys.Date(),".Rda", sep=""))






# ## EVSI - WOrk in progress (unsure if JAGS/OpenBUGS needed?) ##
# 
{# 
  # inner.loops <- 30
  # outer.loops <- 30
  # 
  # sample.size.vec <- c(1000, 2000, 4000, 8000, 12000, 16000, 20000)
  # 
  # # Select lambda values to be considered
  # evsi.lambda <- seq(from = 0, to = 40000, by = 2000)
  # 
  # # Sample all probabilistic parameters
  # clin.char.sims <- gen.clinical.characteristics()[[2]]
  # disability.placebo.sims <- gen.clinical.characteristics()[[5]]
  # disability.txa.sims <- disability.placebo.sims # same between groups
  # utility.sims <- gen.utility.sims()
  # costs.sims <- gen.costs()[[2]]
  # 
  # 
  # # Develop results matrices
  # inner.results <- matrix(0, inner.loops, 4)
  # 
  # evsi.results.placebo <- matrix(0, ncol = length(evsi.lambda), nrow = outer.loops)
  # colnames(evsi.results.placebo) <- as.character(evsi.lambda)
  # evsi.results.txa <- evsi.results.placebo
  # 
  # evsi.results <- matrix(0, nrow = length(sample.size.vec), ncol = length(evsi.lambda + 1))
  # 
  # # EVSI functions
  # 
  # gen.evsi.results <- function(evsi.results1 = evsi.results.placebo, evsi.results2 = evsi.results.txa, lam = evsi.lambda){
  # 
  #   ## calculate the mean NMB for placebo and txa, at each lambda
  #   current.info1 <- apply(evsi.results1, 2, mean)
  #   current.info2 <- apply(evsi.results2, 2, mean)
  # 
  #   current.info <- pmax(current.info1, current.info2)
  # 
  #   evsi.array <- array(0, dim = c(outer.loops, length(lam), 2))
  #   evsi.array[,,1] <- evsi.results1
  #   evsi.array[,,2] <- evsi.results2
  # 
  #   perf.info.sims <- apply(evsi.array, c(1,2), max)
  #   perf.info <- apply(perf.info.sims, 2, mean)
  # 
  #   evsi.results <- c(perf.info - current.info)
  # 
  #   #evsi <- data.frame(lam, evsi.results)
  #   evsi <- c(evsi.results)
  #   names(evsi.results) <- lam
  #   return(evsi)
  # 
  # }
  # 
  # 
  # est.trial.parameters <- function(n, hi.risk.mean, treatment.effect.mean, sample.size){
  #   
  #   ## Hi risk 
  #   se <- sqrt((hi.risk.mean * (1 - hi.risk.mean)) / sample.size)
  #   alpha <- hi.risk.mean * (hi.risk.mean * (1 - hi.risk.mean) / (se^2) - 1 )
  #   beta <- ( alpha /  hi.risk.mean) - alpha 
  #   hi.risk.sim <- rbeta(n, alpha, beta)
  #   
  #   # Treatment effect (based on probabiity of head injury death, and times tx effect selected)
  #   
  #   treatment.effect.mean <- rnorm(1,0.7,0.1)
  #   
  #   control.deaths <- hi.risk.sim * sample.size
  #   intervention.deaths <- 
  #   
  #   
  #   return(list(hi.risk.sim, 
  #               treatment.effect.sim))
  #   
  # }
  # 
  # est.norm <- function(n, mean, sample.size){
  #   
  #   se <- sqrt((mean * (1 - mean)) / sample.size)
  #   alpha <- mean * (mean * (1 - mean) / (se^2) - 1 )
  #   beta <- ( alpha /  mean) - alpha 
  #   sim <- rbeta(n, alpha, beta)
  #   
  #   return(c(sim))
  #   
  # }
  # 
  # 
  # # Sample size loops  
  # 
  # for(z in 1:length(sample.size.vec)){
  #   
  #   sample.size <- sample.size.vec[z]  
  #   
  # 
  #   
  # 
  #   
  #   # Outer loop - Head injury risk, TXA treatment effect, Clinical outcomes, costs of treatment and hospital 
  #   
  #   for(a in 1:outer.loops){
  #   
  #     ## 1. Select study parameters 
  #     clin.sim <- unlist(clin.char.sims[a,])
  #     dis.placebo.sim <- unlist(disability.placebo.sims[a,])
  #     dis.txa.sim <- unlist(disability.txa.sims[a,])
  #     cost.sim <- unlist(costs.sims[a,])
  #     
  #     # Sample a value and create dataset around it 
  #     
  #     
  #     # New hi.death.risk 
  #     hi.risk.sample.prior <- unlist(clin.char.sims[a,2])
  #     hi.risk.posterior <- est.beta(n = inner.loops, hi.risk.sample.prior, sample.size)
  # 
  #     # New hi.death.risk 
  #     treatment.effect.sample.prior <- unlist(clin.char.sims[a,1])
  #     treatment.effect.posterior <- #est.beta(n = inner.loops, hi.risk.sample.mean, sample.size)
  #     
  #       
  #     # create dataset 
  #     #new.sims <- rbeta(1,1,1)
  #     
  #     
  #       # Inner loop
  #       
  #       for(b in 1:inner.loops){
  #     
  #         # New samples from posterior dist
  #         #clin.sim <- unlist(clin.char.sims[b,])
  #         clin.sim[2] <- hi.risk.posterior[b]
  #         
  #         # Standard old samples
  #         clin.sim[4:5] <- unlist(clin.char.sims[b,4:5]) 
  #         utility.sim <- unlist(utility.sims[b,])
  #         cost.sim[4:9] <- unlist(costs.sims[b,4:9])
  #         
  #         ## new distribution - draw samples
  #         
  #         inner.results[b,] <- run.model(clin.sim, dis.placebo.sim, dis.txa.sim, utility.sim, cost.sim)[[1]] 
  #     }
  #   
  #     #after each inner loop PSA, calculate the mean NMB for each tx and store the results
  #     nmb <- gen.nmb(inner.results, lam = evsi.lambda)
  #     evsi.results.placebo[a,] <- nmb[[1]]
  #     evsi.results.txa[a,] <- nmb[[2]]
  #   
  #   }
  # 
  # # Calculate the EVSI for given sample size
  # evsi.sim.results <- gen.evsi.results(evsi.results.placebo, evsi.results.txa, evsi.lambda)
  # 
  # # Save the results of appropriate sample and loop back around
  # evsi.results[z,] <- c(sample.size, evsi.sim.results)
  # 
  # }
  # 
  # 
  # head(evsi.results)
  # 
  # 
  # 
}