## Deterministic analyses ## 

source("CRASH-Mild R Model.R")


## Head injury risk x Treatment effect threshold 

clin.char.dsa <- clin.char
tx.effect.goalseek <- function(x){
  clin.char.dsa[1] <- x 
  z <- run.model(clin.char.dsa, output.type = 2) - 20000
  return(z)
}
hi.risk.thresholds <- seq(from = 0.0001, to = 0.005, by = 0.00001)
threshold.mat <- matrix(0, ncol = 4, nrow = length(hi.risk.thresholds))
threshold.mat[,1] <- hi.risk.thresholds
colnames(threshold.mat) <- c("Risk", "70 years old","80 years old","90 years old")

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

threshold.long <- as.data.frame(threshold.mat)  %>% gather(Age, Threshold, 2:4)


# Plot 

gen.threshold.dsa <- function(results, save = FALSE) {

  plot = ggplot(results) + geom_line(aes(x=Risk, y=Threshold, colour = Age), size=0.6) + 
    labs(x = "Mortality risk following mild TBI (28-day risk)", text = element_text(size=4)) + 
    labs(y = "TXA mortality risk ratio", text = element_text(size=4)) + theme_classic() +
    theme(legend.title = element_blank(), axis.title=element_text(face="bold"), 
          axis.title.x = element_text(margin = margin(t = 7, r = 0, b = 3, l = 0)), 
          axis.title.y = element_text(margin = margin(t = 0, r = 7, b = 0, l = 3)), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          legend.key.width=unit(1.8,"line"), text = element_text(size=7),
          plot.margin=unit(c(0.5,0.5,0,0.5),"cm")) + 
    scale_x_continuous(labels = scales::percent, limits = c(0, 0.005), expand = c(0, 0)) + 
    scale_y_continuous(limits = c(0.6,1), expand = c(0, 0)) 
    
  if(save == TRUE) ggsave(paste("figures\\Threshold-Risk.TxEffect",Sys.Date(),".png"), plot, width=180, height=100, dpi=300, units='mm')
  
  return(plot)  
  
  
}

gen.threshold.dsa(threshold.long, save = TRUE)


