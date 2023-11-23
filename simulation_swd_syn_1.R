#code to produce Fig 4




source('crt.R')

# Number of runs
num_iters<- 1000


run_exps <- function(num_individuals,num_steps,lag,effect_vector){
  
    #run one experiment given sample size, number of steps, lag of the effect and true effect parameters
    
    vec1 <- NULL
    vec2 <- NULL
    vec3 <- NULL
    vec4 <- NULL
    x <- rep(x = NA, times = num_iters)
    pb <- txtProgressBar(0, length(x), style = 3)
    for (seed in 1:num_iters){
      setTxtProgressBar(pb, seed)
      x[seed] 
      Sys.sleep(time = 0.001)
      
      data_xyz <- synthetic(num_individuals,num_steps,effect_vector,seed,normalize=TRUE,include_x=TRUE) 

      data_list_1 <- data_swd(data_xyz,lag,sample_split=FALSE)
      p_list_1 <- crt_swd(data_list_1,type='p_values')
      result1 <- global_test(data_list_1, p_list_1, method='holms')
      result2 <- global_test(data_list_1, p_list_1, method='mcrts+')

      
      data_list_2 <-  data_swd(data_xyz,lag,sample_split=TRUE)
      p_list_2 <- crt_swd(data_list_2,type='p_values')
      result3 <- global_test(data_list_2, p_list_2, method='mcrts')
      result4 <- global_test(data_list_2, p_list_2, method='mcrts+')
      
      vec1 <- c(vec1,result1) 
      vec2 <- c(vec2,result2) 
      vec3 <- c(vec3,result3) 
      vec4 <- c(vec4,result4) 
    }
  return(list(vec1,vec2,vec3,vec4))
}








#vary sample size, fix number of steps, effect lag and effect size
sample_list <- c(100, 200,300,400,500)
num_steps <- 8
lag <- 2
effect_vector <- rep(0, num_steps)

for (i in 1:2){
  
  if (i==2){
    effect_vector[lag+1] = 0.03
  }
  mat1 <- matrix(0, num_iters,5)
  mat2 <- matrix(0, num_iters,5)
  mat3 <- matrix(0, num_iters,5)
  mat4 <- matrix(0, num_iters,5)
  for (j in 1:5){
    num_individuals <- sample_list[j]
    cat("\n","num_individuals=",num_individuals,"\n")
    out <- run_exps(num_individuals,num_steps,lag,effect_vector)
    mat1[,j] <- out[[1]]
    mat2[,j] <- out[[2]]
    mat3[,j] <- out[[3]]
    mat4[,j] <- out[[4]]
    cat("\n")
    cat("rejection_rate_1:",colMeans(mat1))
    cat(" rejection_rate_2:",colMeans(mat2))
    cat(" rejection_rate_3:",colMeans(mat3))
    cat(" rejection_rate_4:",colMeans(mat4))
  }
  save(mat1,mat2,mat3,mat4, file= paste0("result_swd/syn_1_samples_",i, ".RData"))
}







#vary number of steps, fix sample size, effect lag and effect size
num_individuals <- 300
steps_list <- c(4,6,8,10,12)
lag <- 2
effect_vector <- rep(0, num_steps)

for (i in 1:2){
  
  mat1 <- matrix(0, num_iters,5)
  mat2 <- matrix(0, num_iters,5)
  mat3 <- matrix(0, num_iters,5)
  mat4 <- matrix(0, num_iters,5)
  for (j in 1:5){
    num_steps<- steps_list[j]
    cat("\n","num_steps=",num_steps,"\n")
    if (i==2){
      effect_vector <- rep(0, num_steps)
      effect_vector[lag+1] = 0.03
    }
    out <- run_exps(num_individuals,num_steps,lag,effect_vector)
    mat1[,j] <- out[[1]]
    mat2[,j] <- out[[2]]
    mat3[,j] <- out[[3]]
    mat4[,j] <- out[[4]]
    cat("\n")
    cat("rejection_rate_1:",colMeans(mat1))
    cat(" rejection_rate_2:",colMeans(mat2))
    cat(" rejection_rate_3:",colMeans(mat3))
    cat(" rejection_rate_4:",colMeans(mat4))
  }
  save(mat1,mat2,mat3,mat4, file= paste0("result_swd/syn_1_steps_",i, ".RData"))
}








#vary effect lag, fix sample size, number of steps and effect size
num_individuals <- 300
num_steps <- 8
lag_list <- 0:4
effect_vector <- rep(0, num_steps)

for (i in 1:2){
  
  mat1 <- matrix(0, num_iters,5)
  mat2 <- matrix(0, num_iters,5)
  mat3 <- matrix(0, num_iters,5)
  mat4 <- matrix(0, num_iters,5)
  for (j in 1:5){
    lag <- lag_list[j]
    
    if (i==1){
      effect_vector <- rep(0, num_steps)
    }else {
      effect_vector <- rep(0, num_steps)
      effect_vector[lag+1] = 0.03
    }
    cat("\n","lag=",lag,"\n")
    out <- run_exps(num_individuals,num_steps,lag,effect_vector)
    mat1[,j] <- out[[1]]
    mat2[,j] <- out[[2]]
    mat3[,j] <- out[[3]]
    mat4[,j] <- out[[4]]
    cat("\n")
    cat("rejection_rate_1:",colMeans(mat1))
    cat(" rejection_rate_2:",colMeans(mat2))
    cat(" rejection_rate_3:",colMeans(mat3))
    cat(" rejection_rate_4:",colMeans(mat4))
  }
  save(mat1,mat2,mat3,mat4, file= paste0("result_swd/syn_1_lags_",i, ".RData"))
}




#vary effect size, fix sample size, number of steps and effect lag
num_individuals <- 300
num_steps <- 8
lag <- 2
effect_vector <- rep(0, num_steps)
size_list<- c(0.01, 0.02, 0.03, 0.04,0.05)

for (i in 2:2){
  
  mat1 <- matrix(0, num_iters,5)
  mat2 <- matrix(0, num_iters,5)
  mat3 <- matrix(0, num_iters,5)
  mat4 <- matrix(0, num_iters,5)
  for (j in 1:5){
    if (i==2){
      effect_vector[lag+1] = size_list[j]
    }
    cat("\n","effect_size=",size_list[j],"\n")
    out <- run_exps(num_individuals,num_steps,lag,effect_vector)
    mat1[,j] <- out[[1]]
    mat2[,j] <- out[[2]]
    mat3[,j] <- out[[3]]
    mat4[,j] <- out[[4]]
    cat("\n")
    cat("rejection_rate_1:",colMeans(mat1))
    cat(" rejection_rate_2:",colMeans(mat2))
    cat(" rejection_rate_3:",colMeans(mat3))
    cat(" rejection_rate_4:",colMeans(mat4))
  }
  save(mat1,mat2,mat3,mat4, file= paste0("result_swd/syn_1_sizes_",i, ".RData"))
}

  









####plot the results####



library(tidyverse)

plot_fig <- function(link, x_name, x_numbers, to_do,title){
  
  load(link)
  data_ <- data.frame(Numbers = x_numbers,
                      Bonferroni = colMeans(mat1),
                      Fisher = colMeans(mat3),
                      Ours = colMeans(mat4))
  
  if(to_do == "Type_1") {
    p <- ggplot(data_, aes(x=Numbers))+
      geom_hline(yintercept = 0.10)+
      scale_y_continuous(breaks=c(0, 0.10, 0.30, 0.50), limits = c(0.0,0.50))+
      geom_line(aes(y = Bonferroni,color = "Bonferroni"),size=2) +
      geom_line(aes(y = Fisher, color =  "MCRTs+F"),size=2) +
      geom_line(aes(y = Ours,color =  "MCRTs+Z"),size=2) +
      geom_point(aes(y = Bonferroni,color = "Bonferroni"),size=4)+
      geom_point(aes(y = Fisher,color =  "MCRTs+F"),size=4)+
      geom_point(aes(y = Ours,color = "MCRTs+Z"),size=4)+
      #facet_grid(. ~ title)+
      scale_color_manual(name = NULL, values = c("Bonferroni" = "#999999", "MCRTs+F" = "#009E73", "MCRTs+Z" = "#D55E00"))+
      xlab(x_name)+
      ylab("Type I error rate")
    p <- p + theme_bw(base_line_size =0.5,base_rect_size = 1.0)
    p <-p + theme(
      legend.position="top",
      #legend.key.size = unit(0.72, 'cm'), #change legend key size
      #legend.key.height = unit(0.73, 'cm'), #change legend key height
      #legend.key.width = unit(1.6, 'cm'), #change legend key width
      #legend.background = element_rect(color = "black"),
      #legend.position = c(.79,.80),
      text = element_text(size=32, face="bold"),
      #lot.title = element_text(color="black", size=22, face="bold",hjust = 0.5),
      axis.title.x = element_text(color="black", size=32, face="bold",vjust = -2),
      axis.title.y = element_text(color="black", size=32, face="bold",vjust = 5),
      plot.margin=unit(c(1.0,1.0,1.0,1.0),"cm"))
  } else if (to_do == "Power"){
    
    p <- ggplot(data_, aes(x=Numbers))+
      geom_line(aes(y = Bonferroni,color = "Bonferroni"),size=2) +
      geom_line(aes(y = Fisher, color =  "MCRTs+F"),size=2) +
      geom_line(aes(y = Ours,color =  "MCRTs+Z"),size=2) +
      geom_point(aes(y = Bonferroni,color = "Bonferroni"),size=4)+
      geom_point(aes(y = Fisher,color =  "MCRTs+F"),size=4)+
      geom_point(aes(y = Ours,color = "MCRTs+Z"),size=4)+
      #facet_grid(. ~ title)+
      scale_color_manual(name = NULL, values = c("Bonferroni" = "#999999", "MCRTs+F" = "#009E73", "MCRTs+Z" = "#D55E00"))+
      xlab(x_name)+
      ylab("Power")
    p <- p + theme_bw(base_line_size =0.5,base_rect_size = 1.0)
    p <-p + theme(
      legend.position="top",
      #legend.key.size = unit(0.72, 'cm'), #change legend key size
      #legend.key.height = unit(0.73, 'cm'), #change legend key height
      #legend.key.width = unit(1.6, 'cm'), #change legend key width
      #legend.background = element_rect(color = "black"),
      #legend.position = c(.79,0.20),
      text = element_text(size=32, face="bold"),
      #lot.title = element_text(color="black", size=22, face="bold",hjust = 0.5),
      axis.title.x = element_text(color="black", size=32, face="bold",vjust = -2),
      axis.title.y = element_text(color="black", size=32, face="bold",vjust = 5),
      plot.margin=unit(c(1.0,1.0,1.0,1.0),"cm"))
  } else{
    stop("Error")
  }
  ggsave(p,  width=6,height=6, dpi = 300, filename =paste("result_swd/",x_name,"_",to_do,".png"))
  return(p)
}



tital = "Time steps = 8, Time lag = 2"

p <- plot_fig(link="result_swd/syn_1_samples_1.RData",  x_name="Number of units", x_numbers=c(100,200,300,400,500),to_do = "Type_1",title)
p <- plot_fig(link="result_swd/syn_1_samples_2.RData",  x_name="Number of units", x_numbers=c(100,200,300,400,500),to_do = "Power",title)

tital = "Number of units = 750, Time lag = 2"

p <- plot_fig(link="result_swd/syn_1_steps_1.RData",  x_name="Number of time steps", x_numbers=c(4,6,8,10,12),to_do = "Type_1",title)
p <- plot_fig(link="result_swd/syn_1_steps_2.RData",  x_name="Number of time steps", x_numbers=c(4,6,8,10,12),to_do = "Power",title)

tital = "Time steps = 8, Number of units = 750"

p <- plot_fig(link="result_swd/syn_1_lags_1.RData",  x_name="Time lag", x_numbers=c(0,1,2,3,4),to_do = "Type_1",title)
p <- plot_fig(link="result_swd/syn_1_lags_2.RData",  x_name="Time lag", x_numbers=c(0,1,2,3,4),to_do = "Power",title)


tital = "Time steps = 8, Number of units = 750"

p <- plot_fig(link="result_swd/syn_1_sizes_2.RData",  x_name="Effect size", x_numbers=c(0.01,0.02,0.03,0.04,0.05),to_do = "Power",title)



