## Deterministic analyses ## 

source("CRASH-Mild R Model.R")


## DSA  ## 

clin.char.dsa <- clin.char
costs.dsa <- costs

base.case.icer <- run.model(clin.char, disability.placebo, disability.txa, utility, costs, utility.decrement, ae.placebo, ae.txa)[[2]]


# Treatment effect
dsa.res <- rep(NA, 6)
dsa.names <- dsa.res
disability.placebo.dsa <- disability.placebo
disability.txa <- disability.placebo.dsa


dsa.names[1] <- "Risk ratio = 0.95"
clin.char.dsa[1] <- 0.95
dsa.res[1] <- run.model(clin.char.dsa, disability.placebo, disability.txa, utility, costs, utility.decrement, ae.placebo, ae.txa)[[2]]

dsa.names[2] <- "Risk ratio = 0.98"
clin.char.dsa[1] <- 0.98
dsa.res[2] <- run.model(clin.char.dsa, disability.placebo, disability.txa, utility, costs, utility.decrement, ae.placebo, ae.txa)[[2]]

dsa.names[3] <- "Head injury risk = Howard et al (3.7%)"
clin.char.dsa <- clin.char
clin.char.dsa[2] <- 0.0372449
dsa.res[3] <- run.model(clin.char.dsa, disability.placebo, disability.txa, utility, costs, utility.decrement, ae.placebo, ae.txa)[[2]]
clin.char.dsa <- clin.char

dsa.names[4] <- "High cost of injecting equipment (£50)"
costs.dsa[1] <- (1.62) + 50
dsa.res[4] <- run.model(clin.char, disability.placebo, disability.txa, utility, costs.dsa, utility.decrement, ae.placebo, ae.txa)[[2]]
costs.dsa <- costs

dsa.names[5] <- "Standardised mortality ratio following mild TBI applied for 1 year only"
clin.char.dsa <- clin.char
clin.char.dsa[5] <- 1
dsa.res[5] <- run.model(clin.char.dsa, disability.placebo, disability.txa, utility, costs, utility.decrement, ae.placebo, ae.txa)[[2]]

dsa.names[6] <- "Post-discharge costs applied for 1 year only"
costs.dsa[7:9] <- 0
dsa.res[6] <- run.model(clin.char, disability.placebo, disability.txa, utility, costs.dsa, utility.decrement, ae.placebo, ae.txa)[[2]]

dsa.results <- data.frame(Scenario = dsa.names, ICER = dsa.res)



## Head injury risk x Treatment effect threshold 

clin.char.dsa <- clin.char

tx.effect.goalseek <- function(x){
  clin.char.dsa[1] <- x 
  z <- run.model(clin.char.dsa, output.type = 2) - 20000
  return(z)
}

hi.risk.thresholds <- seq(from = 0.0005, to = 0.01, by = 0.00001)
threshold.mat <- matrix(0, ncol = 5, nrow = length(hi.risk.thresholds))
threshold.mat[,1] <- hi.risk.thresholds
colnames(threshold.mat) <- c("Risk", "60 years old", "70 years old","80 years old","90 years old")

for(a in 1:4){
  age.v <- c(60,70,80,90)
  age <- age.v[a]
  time.horizon = min(60, 100-ceiling(age))
  acm <- gen.acm()
  utility.decrement <- gen.utility.dec()
  discount.c <- matrix(1/(1+disc.c)^(0:time.horizon), time.horizon + 1, 2)
  discount.o <- matrix(1/(1+disc.o)^(0:time.horizon), time.horizon + 1, 2)
  
  
    for(i in 1:length(hi.risk.thresholds)){
      clin.char.dsa[2] <- hi.risk.thresholds[i] 
      res <- uniroot(tx.effect.goalseek, c(0.0001,0.9999999))
      threshold.mat[i,a+1] <- res$root
  }
} 

threshold.long <- as.data.frame(threshold.mat)  %>% gather(Age, Threshold, 2:5)


# Plot 

gen.threshold.dsa <- function(results, save = save.results) {

  plot = ggplot(results) + geom_line(aes(x=Risk, y=Threshold, colour = Age), size=0.6) + 
    labs(x = "Mortality risk following mild TBI (28-day risk)", text = element_text(size=4)) + 
    labs(y = "TXA mortality risk ratio", text = element_text(size=4)) + theme_classic() +
    theme(legend.title = element_blank(), axis.title=element_text(face="bold"), 
          axis.title.x = element_text(margin = margin(t = 7, r = 0, b = 3, l = 0)), 
          axis.title.y = element_text(margin = margin(t = 0, r = 7, b = 0, l = 3)), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          legend.key.width=unit(1.8,"line"), text = element_text(size=7),
          plot.margin=unit(c(0.5,0.5,0,0.5),"cm")) + 
    scale_x_continuous(limits = c(min(hi.risk.thresholds), max(hi.risk.thresholds)), breaks = seq(0, 0.1, 0.0005), expand = c(0, 0), 
                       labels = scales::percent_format(accuracy = 0.01L)) + 
    scale_y_continuous(limits = c(0.6,1), breaks = seq(0, 1, 0.05), expand = c(0, 0)) 
    
  if(save == TRUE) ggsave(paste("figures\\Threshold-Risk.TxEffect",Sys.Date(),".png"), plot, width=180, height=100, dpi=300, units='mm')
  
  return(plot)  
  
}

gen.threshold.dsa(threshold.long, FALSE)


