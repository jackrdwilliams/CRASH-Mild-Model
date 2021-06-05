##------------------------------##
##           EVPPI              ## 
##------------------------------##

source("PSA.R")

inner.loops <- 1000
outer.loops <- 1000

# Generate matrices for EVPPI results to be stored in

inner.results <- matrix(0, inner.loops, 4)
evppi.results.placebo <- matrix(0, ncol = length(lambda), nrow = outer.loops)
colnames(evppi.results.placebo) <- as.character(lambda)
evppi.results.txa <- evppi.results.placebo


# EVPPI functions # 

gen.nmb <- function(results, lam = lambda){
  
  
  nmb.table <- matrix(c(lam), ncol = length(lam), nrow = dim(results)[1],  byrow = TRUE) 
  
  p <- ((results[,2] * nmb.table) - results[,1])  
  t <- ((results[,4] * nmb.table) - results[,3])
  
  nmb.p <- apply(p, 2, mean)
  nmb.t <- apply(t, 2, mean) 
  
  #colnames(nmb.p) <- as.character(lam)
  #colnames(nmb.t) <- as.character(lam)
  
  
  return(list(nmb.p, nmb.t))
  
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

gen.evppi.graph = function(evppi, save = FALSE) {
  
  
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length=n+1)
    hcl(h=hues, l=65, c=100)[1:n]
  }
  
  adj_names = sort(setdiff(unique(evppi$Parameters), "EVPI"))
  values = gg_color_hue(length(adj_names))
  names(values) = adj_names
  values = c(values, c(EVPI="black"))
  
  z = ggplot(evppi) + geom_line(aes(x=lambda, y=VoI, colour = Parameters), size=0.6) + 
    labs(x = "Willingness to pay (£)", text = element_text(size=4)) + 
    labs(y = "EVPPI (£)", text = element_text(size=4)) + theme_classic() +
    theme(legend.title = element_blank(), axis.title=element_text(face="bold"), 
          axis.title.x = element_text(margin = margin(t = 7, r = 0, b = 3, l = 0)), 
          axis.title.y = element_text(margin = margin(t = 0, r = 7, b = 0, l = 3)), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          legend.key.width=unit(1.8,"line"), text = element_text(size=7),
          plot.margin=unit(c(0.2,0.5,0,0.2),"cm")) + 
    scale_x_continuous(labels = scales::comma, breaks = c(seq(0,100000,5000)), limits = c(0,max(evppi$lambda)), expand = c(0, 0.1)) + 
    scale_y_continuous(labels = scales::comma, breaks = c(seq(0,400000000,4000000)), limits = c(0,max(evppi$VoI)*1.08), expand = c(0, 0)) + 
    scale_colour_manual(values=values)
  
    geom_vline(xintercept = 20000, linetype="dotted", size=0.25)
  
    if(save == TRUE) ggsave(paste("figures\\EVPPI",Sys.Date(),".png"), z, width=165, height=90, dpi=300, units='mm')
  
  return(z)  
  
}

gen.evppi.trial.graph = function(evppi, save = FALSE) {
  
  z = ggplot(evppi) + geom_line(aes(x=lambda, y=VoI, colour = Parameters, linetype = Parameters), size=0.6) + 
    labs(x = "Willingness to pay (£)", text = element_text(size=4)) + 
    labs(y = "EVPPI (£)", text = element_text(size=4)) + theme_classic() +
    theme(legend.title = element_blank(), axis.title=element_text(face="bold"), 
          axis.title.x = element_text(margin = margin(t = 7, r = 0, b = 3, l = 0)), 
          axis.title.y = element_text(margin = margin(t = 0, r = 7, b = 0, l = 3)), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          legend.key.width=unit(1.8,"line"), text = element_text(size=7),
          plot.margin=unit(c(0.5,0.5,0,0.5),"cm")) + 
    scale_x_continuous(labels = scales::comma, breaks = c(seq(0,100000,5000)), limits = c(0,max(evppi$lambda)), expand = c(0, 0.1)) + 
    scale_y_continuous(labels = scales::comma, breaks = seq(0, 100000000, 4000000), limits = c(0, max(evppi$VoI)*1.08), expand = c(0, 0)) + 
    scale_linetype_manual(values=c("solid", "longdash")) + 
    scale_color_manual(values=c("black", "black")) + 
    geom_vline(xintercept = 20000, linetype="dotted", size=0.25)
  
  
  if(save == TRUE) ggsave(paste("figures\\EVPPI-Trial",Sys.Date(),".png"), z, width=180, height=100, dpi=300, units='mm')
  
  return(z)  
  
}





#### EVPPI for trial parameters ####

pb.trial = txtProgressBar(min = 0, max = outer.loops, initial = 0, style = 3)
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
  setTxtProgressBar(pb.trial,a)
}

evppi.trial.parms <- gen.evppi.results(evppi.results.placebo, evppi.results.txa, lambda)


## Reshape / plot

evppi.trial <- cbind(evpi, evppi.trial.parms[,2])
colnames(evppi.trial) <- c('lambda', 'EVPI', 'EVPPI for trial parameters')
evppi.trial.long <- evppi.trial %>% gather(Parameters, VoI, 2:3)
evppi.trial.long.pop <- evppi.trial.long
evppi.trial.long.pop$VoI <- evppi.trial.long$VoI * effective.population
subset(evppi.trial.long.pop, lambda==20000)
# Plots

#gen.evppi.graph(evppi.long.pop)
gen.evppi.trial.graph(evppi.trial.long.pop, TRUE)
save(evppi.trial.long, file=paste("stored results/evppi.trial.",outer.loops, ".", inner.loops, Sys.Date(),".Rda", sep=""))
subset(evppi.trial.long.pop, lambda==20000)




#### EVPPI for individual loops - 'Double Monte Carlo loop method' 

pb = txtProgressBar(min = 0, max = outer.loops*5, initial = 0, style = 3)

## EVPPI loops - TXA treatment effect
for(a in 1:outer.loops){
  
  ## 1. Select the 'partial' parameter from the outer loop 
  clin.sim <- unlist(clin.char.sims[a,])
  
  for(b in 1:inner.loops){
    
    # Select traditional parameters, minus the outer loop parameter
    
    clin.sim[2:5] <- unlist(clin.char.sims[b,2:5]) # Keep SMRs in PSA 
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
  setTxtProgressBar(pb,a)
}
evppi.tx.effect <- gen.evppi.results(evppi.results.placebo, evppi.results.txa, lambda)


## EVPPI loops - Head injury
for(a in 1:outer.loops){
  
  ## 1. Select the 'partial' parameter from the outer loop 
  clin.sim <- unlist(clin.char.sims[a,])
  
  for(b in 1:inner.loops){
    
    # Select traditional parameters, minus the outer loop parameter
    
    clin.sim[c(1,3,4,5)] <- unlist(clin.char.sims[b,c(1,3,4,5)]) # Keep SMRs in PSA 
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
  setTxtProgressBar(pb,a)
}
evppi.head.injury <- gen.evppi.results(evppi.results.placebo, evppi.results.txa, lambda)

## EVPPI loops - SMRs
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
  setTxtProgressBar(pb,outer.loops+a)
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
  setTxtProgressBar(pb,outer.loops*2 + a)
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
  setTxtProgressBar(pb,outer.loops*3 + a)
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
    
    inner.results[b,] <- run.model(clin.sim, dis.placebo.sim, dis.txa.sim, utility.sim, cost.sim)[[1]] 
  }
  
  #after each inner loop PSA, calculate the mean NMB for each tx and store the results
  nmb <- gen.nmb(inner.results)
  evppi.results.placebo[a,] <- nmb[[1]]
  evppi.results.txa[a,] <- nmb[[2]]
  setTxtProgressBar(pb,outer.loops*4 + a)
}
evppi.costs <- gen.evppi.results(evppi.results.placebo, evppi.results.txa, lambda)


## Reshaping and plotting ## 

evppi.wide <- data.frame(evpi,
                         evppi.tx.effect[,2], 
                         evppi.head.injury[,2],
                         evppi.smr[,2],
                         evppi.disability[,2],
                         evppi.utility[,2],
                         evppi.costs[,2])
colnames(evppi.wide) <- c('lambda', 'EVPI', 'Treatment effect',"Mortality risk", 'SMR', 'Outcomes post-TBI', 'Utility values', 'Costs')

evppi.long <- evppi.wide %>% gather(Parameters, VoI, 2:8)
evppi.long.pop <- evppi.long
evppi.long.pop$VoI <- evppi.long$VoI * effective.population

evppi.long.pop$Parameters <- factor(evppi.long.pop$Parameters, levels = unique(evppi.long.pop$Parameters))
## Save EVPPI's

save(evppi.long, file=paste("stored results/evppi.parms.new.",Sys.Date(),".Rda", sep=""))
gen.evppi.graph(evppi.long.pop, TRUE)

