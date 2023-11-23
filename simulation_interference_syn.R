#code to produce Fig 6

source('crt.R')
require(R.utils)



fisher_combiner <- function(p_values){
  
  # Compute p-values and combined p-values
    
  num_p = length(p_values)
  combined_ps = c()
  for (j in 1:(num_p-1)){
    p_fisher <- fisher(p_values[j:num_p])$p
    combined_ps <- c(combined_ps,p_fisher)
  }
  combined_ps <- c(combined_ps,p_values[num_p])
  return(combined_ps)
}


# Set level, #assignments and #iterations
alpha <- 0.1
num_assignments <- 1000
max_iterations <- 1000

# Treated proportion
prob <- 0.3

# Radius for spillover effect
radius_gap = 1.0
radius_list <- seq(1,5,radius_gap) 
max_radius <- length(radius_list)

# Signal size
signal_list = c(1.0,0.8,0.6,0.4,0.2)*2.0 



# Varying sample sizes


for (N in c(500,400,300,200,100)){
  ours_list <- c()
  fisher_list <- c()
  for (seed in 1:max_iterations){
    p_vec <- c()
    out = NULL
    out <- setup(seed,N,prob,radius_gap,radius_list,signal_list)
    
    for (j in 1:max_radius){
      
      
      Z_status <- subdesign_matrix(out$Z, num_assignments, out$K[[j]])
      p_value <- run_seq_test(Z_status, out$Y, out$I[[j]],out$J[[j]],out$G[[j]], alpha)
      
      if (!is.na(p_value)){
        p_vec <- c(p_vec, p_value)
        
      }
      
    }
    
    if (length(p_vec)==max_radius){
      
      ours_list  <- rbind(ours_list,p_vec)
      fisher_list <- rbind(fisher_list,fisher_combiner(p_vec))
      print(paste(c("indep:", round(colMeans(ours_list<=0.1),3)), collapse="  "))
      print(paste(c("fisher:", round(colMeans(fisher_list<=0.1),3)), collapse="  "))
      cat("\n")
      save(ours_list,fisher_list, file = paste("result_interference/pval-sample-1-",N,".RData",sep = ''))
    }
  }
}
  





# Varying treated proportions

N = 300

for (prob in seq(0.1,0.5,0.1)){
  ours_list <- c()
  fisher_list <- c()
  for (seed in 1:max_iterations){
    p_vec <- c()
    out <- setup(seed,N,prob,radius_gap,radius_list,signal_list)
  
    
    for (j in 1:max_radius){
      Z_status <- subdesign_matrix(out$Z, num_assignments, out$K[[j]])
      p_value <- run_seq_test(Z_status, out$Y, out$I[[j]],out$J[[j]],out$G[[j]], alpha)
      if (!is.na(p_value)){
        p_vec <- c(p_vec, p_value)
      }
      
    }
    
    if (length(p_vec)==max_radius){
      ours_list  <- rbind(ours_list,p_vec)
      fisher_list <- rbind(fisher_list,fisher_combiner(p_vec))
      print(paste(c("indep:", round(colMeans(ours_list<=0.1),3)), collapse="  "))
      print(paste(c("fisher:", round(colMeans(fisher_list<=0.1),3)), collapse="  "))
      cat("\n")
      save(ours_list,fisher_list,  file = paste("result_interference/pval-prob-1-",prob,".RData",sep = ''))
    }

  }

}




# Varying signal sizes

N = 300
prob <- 0.3

for (signal in c(0, 0.5,1,1.5,2)){
  
  signal_list = c(1.0,0.8,0.6,0.4,0.2)*signal*2
  ours_list <- c()
  fisher_list <- c()
  for (seed in 1:max_iterations){
    p_vec <- c()

    out = NULL
    out <- setup(seed,N,prob,radius_gap,radius_list,signal_list)
    
    for (j in 1:max_radius){
      Z_status <- subdesign_matrix(out$Z, num_assignments, out$K[[j]])
      p_value <- run_seq_test(Z_status, out$Y, out$I[[j]],out$J[[j]],out$G[[j]], alpha)
      if (!is.na(p_value)){
        p_vec <- c(p_vec, p_value)
      }
      
    }
    
    if (length(p_vec)==max_radius){
      ours_list  <- rbind(ours_list,p_vec)
      fisher_list <- rbind(fisher_list,fisher_combiner(p_vec))
      print(paste(c("indep:", round(colMeans(ours_list<=0.1),3)), collapse="  "))
      print(paste(c("fisher:", round(colMeans(fisher_list<=0.1),3)), collapse="  "))
      cat("\n")
      save(ours_list,fisher_list, file = paste("result_interference/pval-signal-",signal,".RData",sep = ''))
    }
      

  }
  
}
    























####plot the results####




library(tidyverse)


nam1 = "1"
nam2 = "2"
nam3 = "3"
nam4 = "4"
nam5 = "5"


p_mat  <- c()
for (N in seq(100,500,100)){
  load(paste("result_interference/pval-sample-1-",N,".RData",sep = ''))
  p_mat = cbind(p_mat, colMeans(ours_list<=0.1))
}

data_ <- data.frame(Numbers = seq(100,500,100),
                    e1 = p_mat[1,],
                    e2 = p_mat[2,],
                    e3 = p_mat[3,],
                    e4 = p_mat[4,],
                    e5 = p_mat[5,])

p <- ggplot(data_, aes(x=Numbers))+
  geom_hline(yintercept = 0.1,linewidth=2) +
  scale_y_continuous(breaks = seq(0,0.6,0.1), limits = c(0.0,0.6)) + 
  geom_line(aes(y = e1, color =  nam1),linewidth=4) + geom_line(aes(y = e2, color =  nam2),linewidth=4) +
  geom_line(aes(y = e3, color =  nam3),linewidth=4) + geom_line(aes(y = e4, color =  nam4),linewidth=4) +
  geom_line(aes(y = e5, color =  nam5),linewidth=4) +
  geom_point(aes(y = e1, color = nam1),size=8) + geom_point(aes(y = e2, color =  nam2),size=8) +
  geom_point(aes(y = e3, color =  nam3),size=8) + geom_point(aes(y = e4, color =  nam4),size=8) +
  geom_point(aes(y = e5, color =  nam5),size=8) +
  xlab("Sample size") + ylab("Power") 

p <- p + theme_bw(base_line_size =1,base_rect_size = 0.5)
p <-p + theme(
  legend.position="null",
  text = element_text(size=45, face="bold"),
  axis.title.x = element_text(color="black", size=45, face="bold",vjust = 0.1),
  axis.title.y = element_text(color="black", size=45, face="bold",vjust = 1.9),
  plot.margin=unit(c(1.0,1.0,1.0,1.0),"cm"))

ggsave(p,  width=10,height=10, dpi = 300, filename =paste("result_interference/pval-sample-1.png"))    




p_mat  <- c()
for (prob in seq(0.1,0.5,0.1)){
  load(paste("result_interference/pval-prob-1-",prob,".RData",sep = ''))
  p_mat = cbind(p_mat, colMeans(ours_list<=0.1))
}


data_ <- data.frame(Numbers = seq(0.1,0.5,0.1),
                    e1 = p_mat[1,],
                    e2 = p_mat[2,],
                    e3 = p_mat[3,],
                    e4 = p_mat[4,],
                    e5 = p_mat[5,])

p <- ggplot(data_, aes(x=Numbers))+
  geom_hline(yintercept = 0.1,linewidth=2) +
  scale_y_continuous(breaks = seq(0,0.4,0.1), limits = c(0.0,0.4)) + 
  geom_line(aes(y = e1, color =  nam1),linewidth=4) + geom_line(aes(y = e2, color =  nam2),linewidth=4) +
  geom_line(aes(y = e3, color =  nam3),linewidth=4) + geom_line(aes(y = e4, color =  nam4),linewidth=4) +
  geom_line(aes(y = e5, color =  nam5),linewidth=4) +
  geom_point(aes(y = e1, color = nam1),size=8) + geom_point(aes(y = e2, color =  nam2),size=8) +
  geom_point(aes(y = e3, color =  nam3),size=8) + geom_point(aes(y = e4, color =  nam4),size=8) +
  geom_point(aes(y = e5, color =  nam5),size=8) +
  xlab("Treated proporition") + ylab("Power") 


p <- p + theme_bw(base_line_size =1,base_rect_size = 0.5)
p <-p + theme(
  legend.title = element_blank(),
  legend.position="NULL",
  text = element_text(size=45, face="bold"),
  axis.title.x = element_text(color="black", size=45, face="bold",vjust = 0.1),
  axis.title.y = element_text(color="black", size=45, face="bold",vjust = 1.9),
  plot.margin=unit(c(1.0,1.0,1.0,1.0),"cm"))

ggsave(p,  width=10,height=10, dpi = 300, filename =paste("result_interference/pval-prob-1.png"))    



p_mat  <- c()
for (signal in c(0,0.5,1,1.5,2.0)){
  load(paste("result_interference/pval-signal-",signal,".RData",sep = ''))
  p_mat = cbind(p_mat, colMeans(ours_list<=0.1))
}



data_ <- data.frame(Numbers = c(0,0.5,1,1.5,2.0),
                    e1 = p_mat[1,],
                    e2 = p_mat[2,],
                    e3 = p_mat[3,],
                    e4 = p_mat[4,],
                    e5 = p_mat[5,])

p <- ggplot(data_, aes(x=Numbers))+
  labs(col = "k",size=4) +
  geom_hline(yintercept = 0.1,linewidth=2) +
  scale_y_continuous(breaks = seq(0,0.5,0.1), limits = c(0.0,0.5)) + 
  geom_line(aes(y = e1, color =  nam1),linewidth=4) + geom_line(aes(y = e2, color =  nam2),linewidth=4) +
  geom_line(aes(y = e3, color =  nam3),linewidth=4) + geom_line(aes(y = e4, color =  nam4),linewidth=4) +
  geom_line(aes(y = e5, color =  nam5),linewidth=4) +
  geom_point(aes(y = e1, color = nam1),size=8) + geom_point(aes(y = e2, color =  nam2),size=8) +
  geom_point(aes(y = e3, color =  nam3),size=8) + geom_point(aes(y = e4, color =  nam4),size=8) +
  geom_point(aes(y = e5, color =  nam5),size=8) +
  xlab("Signal strength") + ylab("Power") 


p <- p + theme_bw(base_line_size =1,base_rect_size = 0.5)
p <-p + theme(
  #legend.title = element_blank(),
  legend.position="right",
  text = element_text(size=45, face="bold"),
  axis.title.x = element_text(color="black", size=45, face="bold",vjust = 0.1),
  axis.title.y = element_text(color="black", size=45, face="bold",vjust = 1.9),
  plot.margin=unit(c(1.0,1.0,1.0,1.0),"cm"))

ggsave(p,  width=11.3,height=10, dpi = 300, filename =paste("result_interference/pval-signal-1.png"))    













nam1 = "1"
nam2 = "2"
nam3 = "3"
nam4 = "4"
nam5 = "5"


p_mat  <- c()
for (N in seq(100,500,100)){
  load(paste("result_interference/pval-sample-1-",N,".RData",sep = ''))
  p_mat = cbind(p_mat, colMeans(fisher_list<=0.1))
}

data_ <- data.frame(Numbers = seq(100,500,100),
                    e1 = p_mat[1,],
                    e2 = p_mat[2,],
                    e3 = p_mat[3,],
                    e4 = p_mat[4,],
                    e5 = p_mat[5,])

p <- ggplot(data_, aes(x=Numbers))+
  geom_hline(yintercept = 0.1,linewidth=2) +
  scale_y_continuous(breaks = seq(0,0.8,0.1), limits = c(0.0,0.8)) + 
  geom_line(aes(y = e1, color =  nam1),linewidth=4) + geom_line(aes(y = e2, color =  nam2),linewidth=4) +
  geom_line(aes(y = e3, color =  nam3),linewidth=4) + geom_line(aes(y = e4, color =  nam4),linewidth=4) +
  geom_line(aes(y = e5, color =  nam5),linewidth=4) +
  geom_point(aes(y = e1, color = nam1),size=8) + geom_point(aes(y = e2, color =  nam2),size=8) +
  geom_point(aes(y = e3, color =  nam3),size=8) + geom_point(aes(y = e4, color =  nam4),size=8) +
  geom_point(aes(y = e5, color =  nam5),size=8) +
  xlab("Sample size") + ylab("Power") 

p <- p + theme_bw(base_line_size =1,base_rect_size = 0.5)
p <-p + theme(
  legend.position="null",
  text = element_text(size=45, face="bold"),
  axis.title.x = element_text(color="black", size=45, face="bold",vjust = 0.1),
  axis.title.y = element_text(color="black", size=45, face="bold",vjust = 1.9),
  plot.margin=unit(c(1.0,1.0,1.0,1.0),"cm"))

ggsave(p,  width=10,height=10, dpi = 300, filename =paste("result_interference/pval-sample-2.png"))    




p_mat  <- c()
for (prob in seq(0.1,0.5,0.1)){
  load(paste("result_interference/pval-prob-1-",prob,".RData",sep = ''))
  p_mat = cbind(p_mat, colMeans(fisher_list<=0.1))
}


data_ <- data.frame(Numbers = seq(0.1,0.5,0.1),
                    e1 = p_mat[1,],
                    e2 = p_mat[2,],
                    e3 = p_mat[3,],
                    e4 = p_mat[4,],
                    e5 = p_mat[5,])

p <- ggplot(data_, aes(x=Numbers))+
  geom_hline(yintercept = 0.1,linewidth=2) +
  scale_y_continuous(breaks = seq(0,0.6,0.1), limits = c(0.0,0.6)) + 
  geom_line(aes(y = e1, color =  nam1),linewidth=4) + geom_line(aes(y = e2, color =  nam2),linewidth=4) +
  geom_line(aes(y = e3, color =  nam3),linewidth=4) + geom_line(aes(y = e4, color =  nam4),linewidth=4) +
  geom_line(aes(y = e5, color =  nam5),linewidth=4) +
  geom_point(aes(y = e1, color = nam1),size=8) + geom_point(aes(y = e2, color =  nam2),size=8) +
  geom_point(aes(y = e3, color =  nam3),size=8) + geom_point(aes(y = e4, color =  nam4),size=8) +
  geom_point(aes(y = e5, color =  nam5),size=8) +
  xlab("Treated proporition") + ylab("Power") 


p <- p + theme_bw(base_line_size =1,base_rect_size = 0.5)
p <-p + theme(
  legend.title = element_blank(),
  legend.position="NULL",
  text = element_text(size=45, face="bold"),
  axis.title.x = element_text(color="black", size=45, face="bold",vjust = 0.1),
  axis.title.y = element_text(color="black", size=45, face="bold",vjust = 1.9),
  plot.margin=unit(c(1.0,1.0,1.0,1.0),"cm"))

ggsave(p,  width=10,height=10, dpi = 300, filename =paste("result_interference/pval-prob-2.png"))    



p_mat  <- c()
for (signal in c(0,0.5,1,1.5,2.0)){
  load(paste("result_interference/pval-signal-",signal,".RData",sep = ''))
  p_mat = cbind(p_mat, colMeans(fisher_list<=0.1))
}



data_ <- data.frame(Numbers = c(0,0.5,1,1.5,2.0),
                    e1 = p_mat[1,],
                    e2 = p_mat[2,],
                    e3 = p_mat[3,],
                    e4 = p_mat[4,],
                    e5 = p_mat[5,])

p <- ggplot(data_, aes(x=Numbers))+
  labs(col = "k",size=4) +
  geom_hline(yintercept = 0.1,linewidth=2) +
  scale_y_continuous(breaks = seq(0,0.8,0.1), limits = c(0.0,0.8)) + 
  geom_line(aes(y = e1, color =  nam1),linewidth=4) + geom_line(aes(y = e2, color =  nam2),linewidth=4) +
  geom_line(aes(y = e3, color =  nam3),linewidth=4) + geom_line(aes(y = e4, color =  nam4),linewidth=4) +
  geom_line(aes(y = e5, color =  nam5),linewidth=4) +
  geom_point(aes(y = e1, color = nam1),size=8) + geom_point(aes(y = e2, color =  nam2),size=8) +
  geom_point(aes(y = e3, color =  nam3),size=8) + geom_point(aes(y = e4, color =  nam4),size=8) +
  geom_point(aes(y = e5, color =  nam5),size=8) +
  xlab("Signal strength") + ylab("Power") 


p <- p + theme_bw(base_line_size =1,base_rect_size = 0.5)
p <-p + theme(
  #legend.title = element_blank(),
  legend.position="right",
  text = element_text(size=45, face="bold"),
  axis.title.x = element_text(color="black", size=45, face="bold",vjust = 0.1),
  axis.title.y = element_text(color="black", size=45, face="bold",vjust = 1.9),
  plot.margin=unit(c(1.0,1.0,1.0,1.0),"cm"))

ggsave(p,  width=11.3,height=10, dpi = 300, filename =paste("result_interference/pval-signal-2.png"))    













