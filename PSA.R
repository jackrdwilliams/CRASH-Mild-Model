## Probabilistic sensitivity analysis ## 

source("CRASH-Mild R Model.R")
source("VoI Parameters.R")

## PSA - Base case ## 

sims <- 10000

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

gen.ceac.graph(ceac)
gen.evpi.graph(evpi)

# CEAC - Sensitivity analysis 
# ceac.sens <- gen.ceac.table(psa.results.sens)[[1]]
# gen.ceac.graph(ceac.sens)





# PSA - Sensitivity ##
# 
# gen.ceac.graph.sens = function(psa, save = FALSE) {
# 
#   z = ggplot(psa) + geom_line(aes(x=lambda, y=Probability, colour = Treatment.Effect), size=0.6) +
#     labs(x = "Willingness to pay (£)", text = element_text(size=4)) +
#     labs (y = "Probability cost-effective", text = element_text(size=4)) + theme_classic() +
#     theme(legend.title = element_blank(), axis.title=element_text(face="bold"),
#           axis.title.x = element_text(margin = margin(t = 7, r = 0, b = 3, l = 0)),
#           axis.title.y = element_text(margin = margin(t = 0, r = 7, b = 0, l = 3)),
#           panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#           legend.key.width=unit(1.8,"line"), text = element_text(size=7),
#           plot.margin=unit(c(0.5,0.5,0,0.5),"cm")) +
#     scale_x_continuous(labels = scales::comma, breaks = c(seq(0,100000,5000)), limits = c(0,max(psa$lam)), expand = c(0, 0.1)) +
#     scale_y_continuous(limits = c(0,1), breaks=seq(0,1,0.1), expand = c(0, 0)) +
#     geom_vline(xintercept = 20000, linetype="dotted", size=0.25)
# 
# 
#   if(save == TRUE) ggsave(paste("figures\\CEAC.Sensitivity.Analysis",Sys.Date(),".png"), z, width=140, height=80, dpi=300, units='mm')
# 
# 
#   return(z)
# 
# }
# gen.evpi.graph.sens = function(evpi, save = FALSE) {
#   
#   z = ggplot(evpi) + geom_line(aes(x=lambda, y=VoI, colour = Treatment.Effect), size=0.6) + 
#     labs(x = "Willingness to pay (£)", text = element_text(size=4)) + 
#     labs (y = "EVPI (£)", text = element_text(size=4)) + theme_classic() +
#     theme(legend.title = element_blank(), axis.title=element_text(face="bold"), 
#           axis.title.x = element_text(margin = margin(t = 7, r = 0, b = 3, l = 0)), 
#           axis.title.y = element_text(margin = margin(t = 0, r = 7, b = 0, l = 3)), 
#           panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
#           legend.key.width=unit(1.8,"line"), text = element_text(size=7),
#           plot.margin=unit(c(0.5,0.5,0,0.5),"cm")) + 
#     scale_x_continuous(labels = scales::comma, breaks = c(seq(0,100000,5000)), limits = c(0,max(evpi$lam)), expand = c(0, 0.1)) + 
#     scale_y_continuous(labels = scales::comma, breaks = c(seq(0,100000000,5000000)), expand = c(0, 0)) + 
#     geom_vline(xintercept = 20000, linetype="dotted", size=0.25)
#   
#   if(save == TRUE) ggsave(paste("figures\\EVPI.Sensitivity",Sys.Date(),".png"), z, width=140, height=80, dpi=300, units='mm')
#   
#   return(z)
#   
# }
# 
# psa.results.sens <- matrix(0, sims, 4)
# colnames(psa.results.sens) <- c("cost.placebo", "utility.placebo","cost.txa","utility.txa")
# pb = txtProgressBar(min = 0, max = sims, initial = 0, style = 3)
# 
# tx.effect.alt <- data.frame(a = exp(rnorm(sims, log(0.80), 0.2216285)),
#                             b = exp(rnorm(sims, log(0.90), 0.2216285)),
#                             c = exp(rnorm(sims, log(0.95), 0.2216285)))
# res <- data.frame(lambda = lambda, a = NA, b = NA, c = NA)
# evpi.res <- res
# 
# for(w in 1:3){
#   for(p in 1:sims){
#     # Subset and assign existing sims
#     clin.sim[1] <- tx.effect.alt[p,w]
#     clin.sim[2:5] <- unlist(clin.char.sims[p,2:5])
#     dis.placebo.sim <- unlist(disability.placebo.sims[p,])
#     dis.txa.sim <- unlist(disability.txa.sims[p,])
#     utility.sim <- unlist(utility.sims[p,])
#     cost.sim <- unlist(costs.sims[p,])
# 
#     psa.results.sens[p,] <- run.model(clin.sim, dis.placebo.sim, dis.txa.sim, utility.sim, cost.sim)[[1]]
#     setTxtProgressBar(pb,p)
#   }
# 
# res[,w+1] <- gen.ceac.table(psa.results.sens)[[1]][,2]
# evpi.res[,w+1] <- gen.ceac.table(psa.results.sens)[[2]][,2]
# }
# 
# colnames(res) <- c("lambda", "Risk ratio = 0.8", "Risk ratio = 0.9", "Risk ratio = 0.95")
# colnames(evpi.res) <- c("lambda", "Risk ratio = 0.8", "Risk ratio = 0.9", "Risk ratio = 0.95")
# 
# psa.sens <- res %>% gather(Treatment.Effect, Probability, 2:4)
# gen.ceac.graph.sens(psa.sens, TRUE)
# 
# evpi.s <- evpi.res %>% gather(Treatment.Effect, VoI, 2:4) 
# evpi.sens <- evpi.s 
# evpi.sens[,3] <- evpi.sens[,3] * effective.population
# 
# gen.evpi.graph.sens(evpi.sens, TRUE)
# 

## ANCOVA ## 

sim.parameters <- cbind(clin.char.sims[,c(1,2,4,5)], disability.placebo.sims[,2:5], utility.sims[,2:5], costs.sims[,2:9])

psa.results <- as.data.frame(psa.results)
incr.cost <- psa.results$cost.txa - psa.results$cost.placebo
incr.qaly <- psa.results$utility.txa - psa.results$utility.placebo
incr.nmb <- (incr.qaly * 20000) - incr.cost
icer <- incr.cost / incr.qaly

#### FUNCTIONS #### 
gen.anova.results <- function(model){
  
  anova <- anova(model)
  
  mean.sq <- anova$`Mean Sq` 
  total.sq <- sum(mean.sq) + deviance(model)
  
  prop.total.sq <- mean.sq / total.sq
  prop.vars.sq <- mean.sq / sum(mean.sq)
  
  rownames <- rownames(anova) 
  names <- c(rownames, "Total")
  group <- c("Treatment effect", "Mortality risk", rep("SMR",2), 
             rep("Outcomes post-TBI", 3), rep("Utility",4), "Hospital cost", 
             rep("Post-discharge costs", 6), "Residuals", "Total")

    
  data <- data.frame(var = names, 
                     group = group, 
                     mean.squares = c(mean.sq, total.sq), 
                     proportion.total = c(prop.total.sq,0),
                     proportion.var = c(prop.vars.sq,0))

  dplyr::arrange(data, proportion.var)
  
  rank <- data %>% arrange(desc(proportion.var))
  
  return(suppressWarnings(list(data, rank)))
  
}


## Incremental Costs
data.costs <- cbind(incr.cost, sim.parameters)
model.costs <- lm(incr.cost~. , data.costs)

## Incremental Outcomes
data.qaly <- cbind(incr.qaly, sim.parameters)
model.qaly <- lm(incr.qaly~. , data.qaly)

## Incremental NMB
data.nmb <- cbind(incr.nmb, sim.parameters)
model.nmb <- lm(incr.nmb~. , data.nmb)

## ICER
data.icer <- cbind(icer, sim.parameters)
model.icer <- lm(icer~. , data.icer)

## Results/Output

gen.anova.results(model.costs)[[1]]
gen.anova.results(model.qaly)[[1]]
gen.anova.results(model.nmb)[[1]]

results.cost <- gen.anova.results(model.costs)[[1]]
results.qaly <- gen.anova.results(model.qaly)[[1]]
results.nmb <- gen.anova.results(model.nmb)[[1]]


combine.rank <- function(results) {
  
  res <- results %>% dplyr::select(group, proportion.var) %>% filter(!(group == c("Residual", "Total"))) %>% 
          group_by(group) %>% summarize(proportion.var=sum(proportion.var)) %>% arrange(desc(proportion.var))

  return(res)

  }


rank.cost <- combine.rank(results.cost)
rank.qaly <- combine.rank(results.qaly)
rank.nmb <- combine.rank(results.nmb)


gen.ancova.plot <- function(result){
  
  res <- as.data.frame(result)
  
  res$group <- factor(res$group, levels = rev(res$group))
  ggplot(res, aes(x = proportion.var, y = group)) + geom_bar(stat="identity", width=0.5,  fill = "steelblue") +
    labs(x = "Proportion of variance explained by parameter", text = element_text(size=4)) + 
    theme_classic() + 
    theme(axis.title= element_text(face="bold"), 
          axis.title.x = element_text(margin = margin(t = 7, r = 0, b = 4, l = 0)), 
          axis.title.y = element_blank(),
          plot.margin=unit(c(0.5,0.5,0,0.5),"cm")) +
    scale_x_continuous(labels = scales::percent, limits = c(0,max(res$proportion.var*1.05)), expand = c(0, 0.001)) 

}


gen.ancova.plot(rank.cost)  
gen.ancova.plot(rank.qaly) 
gen.ancova.plot(rank.nmb)  
