## Probabilistic sensitivity analysis ## 

source("CRASH-Mild R Model.R")
source("VoI Parameters.R")

## PSA - Base case ## 

psa.results <- matrix(0, sims, 4)
colnames(psa.results) <- c("cost.placebo", "utility.placebo","cost.txa","utility.txa")
pb = txtProgressBar(min = 0, max = sims, initial = 0, style = 3)

for(p in 1:sims){
  
  # Subset and assign existing sims
  clin.sim <- unlist(clin.char.sims[p,])
  dis.placebo.sim <- unlist(disability.placebo.sims[p,])
  dis.txa.sim <- unlist(disability.txa.sims[p,])
  utility.sim <- unlist(utility.sims[p,])
  cost.sim <- unlist(costs.sims[p,])
  
  psa.results[p,] <- run.model(clin.sim, dis.placebo.sim, dis.txa.sim, utility.sim, cost.sim)[[1]] 
  setTxtProgressBar(pb,p)
}

psa.incr <- data.frame(cost = psa.results[,3] - psa.results[,1], 
                       utility = psa.results[,4] - psa.results[,2])

# Proportion of simulations where TXA less effective
mean((psa.results[,4] - psa.results[,2])<0)


## PSA - Sensitivity ## 

# psa.results.sens <- matrix(0, sims, 4)
# colnames(psa.results.sens) <- c("cost.placebo", "utility.placebo","cost.txa","utility.txa")
# tx.effect.alt <- treatment.effect.sims <- exp(rnorm(sims, log(0.90), 0.2216285))  
# pb = txtProgressBar(min = 0, max = sims, initial = 0, style = 3)
# 
# for(p in 1:sims){
#   
#   # Subset and assign existing sims
#   clin.sim <- unlist(clin.char.sims[p,])
#   clin.sim[1] <- tx.effect.alt[p]
#   dis.placebo.sim <- unlist(disability.placebo.sims[p,])
#   dis.txa.sim <- unlist(disability.txa.sims[p,])
#   utility.sim <- unlist(utility.sims[p,])
#   cost.sim <- unlist(costs.sims[p,])
#   
#   psa.results.sens[p,] <- run.model(clin.sim, dis.placebo.sim, dis.txa.sim, utility.sim, cost.sim)[[1]] 
#   setTxtProgressBar(pb,p)
# }
# 



# Generate CEAC table

gen.ceac.table <- function(results, lam = lambda){
  
  
  ## CEAC ## 
  
  lambda.table <- matrix(lam, ncol = length(lam), nrow = dim(results)[1], byrow = TRUE) 
  inmb.count <- ((results[,4] * lambda.table) - results[,3]) - ((results[,2] * lambda.table) - results[,1]) > 0 
  
  prob.ce <- apply(inmb.count, 2, mean) 
  ceac.table <- data.frame(lam, prob.ce) 
  
  ## EVPI ## 
  
  evpi.table <- matrix(lam, ncol = length(lam), nrow = dim(results)[1], byrow = TRUE) 
  
  nmb.p <- ((results[,2] * evpi.table) - results[,1])  
  nmb.t <- ((results[,4] * evpi.table) - results[,3]) 
  
  av.nmb.p <- apply(nmb.p, 2, mean)
  av.nmb.t <- apply(nmb.t, 2, mean)
  
  evpi.mat <- matrix(0, dim(results)[1], length(av.nmb.p))
  
  for(i in 1:ncol(nmb.p)){
    if(av.nmb.p[i] >= av.nmb.t[i]) evpi.mat[,i] <- nmb.t[,i] - nmb.p[,i]
    if(av.nmb.t[i] >= av.nmb.p[i]) evpi.mat[,i] <- nmb.p[,i] - nmb.t[,i]
  }
  
  evpi.mat[evpi.mat<0] <- 0
  evpi.m <- apply(evpi.mat, 2, mean)
  
  evpi <- data.frame(lam, evpi.m)
  
  return(list(ceac.table,
              evpi))
  
}

ceac <- gen.ceac.table(psa.results)[[1]]
evpi <- gen.ceac.table(psa.results)[[2]]

evpi.pop <- evpi[,2] * effective.population

# Graphics  - TBC (take from other sources)


# CEAC # 

gen.ceac.graph = function(psa, save = FALSE) {
  
  z = ggplot(psa) + geom_line(aes(x=lam, y=prob.ce), size=0.6) + 
    labs(x = "Willingness to pay (£)", text = element_text(size=4)) + 
    labs (y = "Probability cost-effective", text = element_text(size=4)) + theme_classic() +
    theme(legend.title = element_blank(), axis.title=element_text(face="bold"), 
          axis.title.x = element_text(margin = margin(t = 7, r = 0, b = 3, l = 0)), 
          axis.title.y = element_text(margin = margin(t = 0, r = 7, b = 0, l = 3)), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          legend.key.width=unit(1.8,"line"), text = element_text(size=7),
          plot.margin=unit(c(0.5,0.5,0,0.5),"cm")) + 
    scale_x_continuous(labels = scales::comma, breaks = c(seq(0,100000,5000)), limits = c(0,max(psa$lam)), expand = c(0, 0.1)) + 
    scale_y_continuous(limits = c(0,1), breaks=seq(0,1,0.1), expand = c(0, 0)) + 
    geom_vline(xintercept = 20000, linetype="dotted", size=0.25)
  
  
  if(save == TRUE) ggsave(paste("figures\\CEAC",Sys.Date(),".png"), z, width=107, height=70, dpi=300, units='mm')
  
  
  return(z)
  
}
gen.evpi.graph = function(evpi, save = FALSE) {
  
  z = ggplot(evpi) + geom_line(aes(x=lam, y=evpi.m), size=0.6) + 
    labs(x = "Willingness to pay (£)", text = element_text(size=4)) + 
    labs (y = "EVPI (£)", text = element_text(size=4)) + theme_classic() +
    theme(legend.title = element_blank(), axis.title=element_text(face="bold"), 
          axis.title.x = element_text(margin = margin(t = 7, r = 0, b = 3, l = 0)), 
          axis.title.y = element_text(margin = margin(t = 0, r = 7, b = 0, l = 3)), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          legend.key.width=unit(1.8,"line"), text = element_text(size=7),
          plot.margin=unit(c(0.5,0.5,0,0.5),"cm")) + 
    scale_x_continuous(labels = scales::comma, breaks = c(seq(0,100000,5000)), limits = c(0,max(evpi$lam)), expand = c(0, 0.1)) + 
    scale_y_continuous(expand = c(0, 0)) + 
    geom_vline(xintercept = 20000, linetype="dotted", size=0.25)
  
  if(save == TRUE) ggsave(paste("figures\\EVPI",Sys.Date(),".png"), z, width=107, height=70, dpi=300, units='mm')
  
  return(z)
  
}

gen.ceac.graph(ceac, TRUE)
gen.evpi.graph(evpi, TRUE)

# CEAC - Sensitivity analysis 
# ceac.sens <- gen.ceac.table(psa.results.sens)[[1]]
# gen.ceac.graph(ceac.sens)
