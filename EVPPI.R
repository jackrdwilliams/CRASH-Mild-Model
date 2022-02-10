##------------------------------##
##           EVPPI              ## 
##------------------------------##

source("PSA.R")
source("CRASH-Mild R Model.R") # return to base case

inner.loops <- 1000
outer.loops <- 10000

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

gen.evppi.graph = function(evppi, save = save.results) {
  
  
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
    scale_y_continuous(labels = scales::comma, breaks = c(seq(0,400000000,5000000)), limits = c(0,max(evppi$VoI)*1.08), expand = c(0, 0)) + 
    scale_colour_manual(values=values) + 
    geom_vline(xintercept = 20000, linetype="dotted", size=0.25)
  
    if(save == TRUE) ggsave(paste("figures\\EVPPI",Sys.Date(),".png"), z, width=165, height=90, dpi=300, units='mm')
  
  return(z)  
  
}

gen.evppi.trial.graph = function(evppi, save = save.results) {
  
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
    scale_y_continuous(labels = scales::comma, breaks = seq(0, 100000000, 5000000), limits = c(0, max(evppi$VoI)*1.02), expand = c(0, 0)) + 
    scale_linetype_manual(values=c("solid", "longdash")) + 
    scale_color_manual(values=c("black", "black")) + 
    geom_vline(xintercept = 20000, linetype="dotted", size=0.25)
  
  
  if(save == TRUE) ggsave(paste("figures\\EVPPI-Trial",Sys.Date(),".png"), z, width=174, height=105, dpi=300, units='mm')
  
  return(z)  
  
}





#### EVPPI for trial parameters 

pb.trial = txtProgressBar(min = 0, max = outer.loops, initial = 0, style = 3)
for(a in 1:outer.loops){

  ## 1. Select trial parameters
  clin.sim[1:3] <- unlist(clin.char.sims[a,1:3]) # HI risk and tx effect
  dis.placebo.sim <- unlist(disability.placebo.sims[a,])
  dis.txa.sim <- unlist(disability.txa.sims[a,])
  ae.p.sim <- unlist(ae.placebo.sims[a,])
  ae.t.sim <- unlist(ae.txa.sims[a,])
  cost.sim[1:3] <- unlist(costs.sims[a,1:3])
  
  for(b in 1:inner.loops){

    # Select traditional parameters, minus the outer loop parameter
    clin.sim[4:5] <- unlist(clin.char.sims[b,4:5]) # Keep SMRs in PSA
    utility.sim <- unlist(utility.sims[b,]) ## Utility values
    cost.sim[4:9] <- unlist(costs.sims[b,4:9])

    inner.results[b,] <- run.model(clin.sim, dis.placebo.sim, dis.txa.sim, utility.sim, cost.sim, ae.p = ae.p.sim, ae.t = ae.t.sim)[[1]]
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
subset(evppi.trial.long, lambda==20000)
subset(evppi.trial.long.pop, lambda==20000)
# Plots

#gen.evppi.graph(evppi.long.pop)
gen.evppi.trial.graph(evppi.trial.long.pop, TRUE)
save(evppi.trial.long, file=paste("stored results/evppi.trial.",outer.loops, ".", inner.loops, Sys.Date(),".Rda", sep=""))





#### EVPPI for individual loops - 'Double Monte Carlo loop method' 


inner.loops <- 1000
outer.loops <- 1000

inner.results <- matrix(0, inner.loops, 4)
evppi.results.placebo <- matrix(0, ncol = length(lambda), nrow = outer.loops)
colnames(evppi.results.placebo) <- as.character(lambda)
evppi.results.txa <- evppi.results.placebo


parameter.groups <- 7
evppi.wide <- data.frame(lambda = lambda,
                         evpi = evpi[,2], 
                         tx = NA,
                         mort = NA,
                         smr = NA,
                         outcomes = NA,
                         utility = NA, 
                         costs = NA, 
                         ae = NA)
colnames(evppi.wide) <- c('lambda', 'EVPI', 'Treatment effect',"Mortality risk", 'SMR', 'Outcomes post-TBI', 'Utility values', 'Costs', "Adverse events")


pb = txtProgressBar(min = 0, max = outer.loops*7, initial = 0, style = 3)

for(j in 1:parameter.groups){
  ## EVPPI loops - TXA treatment effect
  for(a in 1:outer.loops){
    
    
    for(b in 1:inner.loops){
      
      if(j==1) ev.tx <- a else tx <- b
      if(j==2) ev.mort <- a else ev.mort <- b 
      if(j==3) ev.smr <- a else ev.smr <- b   
      if(j==4) ev.outcomes <- a else ev.outcomes <- b      
      if(j==5) ev.utility <- a else ev.utility <- b
      if(j==6) ev.costs <- a else ev.costs <- b
      if(j==7) ev.ae <- a else ev.ae <- b  
    
      clin.sim[1] <- unlist(clin.char.sims[ev.tx,1])
      clin.sim[2] <- unlist(clin.char.sims[ev.mort,2])
      clin.sim[4:5] <- unlist(clin.char.sims[ev.smr,4:5])
      
      dis.placebo.sim <- unlist(disability.placebo.sims[ev.outcomes,])
      dis.txa.sim <- unlist(disability.txa.sims[ev.outcomes,])
      
      utility.sim <- unlist(utility.sims[ev.utility,])
      cost.sim <- unlist(costs.sims[ev.costs,])
      ae.p.sim <- unlist(ae.placebo.sims[ev.ae,])
      ae.t.sim <- unlist(ae.txa.sims[ev.ae,])
      
      inner.results[b,] <- run.model(clin.sim, dis.placebo.sim, dis.txa.sim, utility.sim, cost.sim, ae.p = ae.p.sim, ae.t = ae.t.sim)[[1]] 
    }
    
    #after each inner loop PSA, calculate the mean NMB for each tx and store the results
    nmb <- gen.nmb(inner.results)
    evppi.results.placebo[a,] <- nmb[[1]]
    evppi.results.txa[a,] <- nmb[[2]]
    setTxtProgressBar(pb, (j-1)*outer.loops + a) 
  }
  
  evppi.wide[,j+2] <- gen.evppi.results(evppi.results.placebo, evppi.results.txa, lambda)[,2]

}



## Reshaping and plotting ## 

evppi.long <- evppi.wide %>% gather(Parameters, VoI, 2:9)
evppi.long.pop <- evppi.long
evppi.long.pop$VoI <- evppi.long$VoI * effective.population

evppi.long.pop$Parameters <- factor(evppi.long.pop$Parameters, levels = unique(evppi.long.pop$Parameters))

## Save EVPPI's

save(evppi.long, file=paste("stored results/evppi.parms.new1.",Sys.Date(),".Rda", sep=""))
gen.evppi.graph(evppi.long.pop)

