#code to produce Fig 5




library(lme4)
library(plm) 
source('crt.R')

#####Simulation to compare CRTs with linear models#####

# Setup
num_individuals <- 200
num_steps <- 8
effect_vector <- rep(0, num_steps)
effect_vector[1] = 0.1
effect_vector[2] = 0.3
effect_vector[3] = 0.6
effect_vector[4] = 0.4
effect_vector[5] = 0.2
num_iters <- 1000


# CRTs
list_crts = list()
list_crts_f = list()
tt = 1
for (index in 1:3){
  cat("\n")
  cat("index=",index)
  
  if (index ==1){
    weight_list <- c(0,0.1)
  }else{
    weight_list <- c(0.1)
  }

  for (weight in weight_list){
    
    cat("\n")
    cat("weight=",weight)
    cat("\n")

    mat1 <- matrix(0, num_iters,5)
    mat2 <- matrix(0, num_iters,5)
    mat3 <- matrix(0, num_iters,5)
    mat4 <- matrix(0, num_iters,5)
    check1 <- 0
    check2 <- 0
    check3 <- 0
    check4 <- 0
    check5 <- 0
    check6 <- 0


    for (seed in 1:num_iters){


      data_xyz <- synthetic(num_individuals,num_steps,effect_vector,seed,weight,index,normalize=TRUE,include_x=TRUE) 
      
      for (lag in 0:4){
        
        
        data_list <- data_swd(data_xyz,lag,sample_split=TRUE)
        
        p_list <- crt_swd(data_list, type='ci',n_beta=500,const=1)
        ci <- global_ci(data_list,p_list,method='mcrts', const=1)
        ci_ <- global_ci(data_list,p_list,method='mcrts+', const=1)
        
        i <- lag + 1

        mat1[seed,i] <- (effect_vector[i]>=ci[1])*(effect_vector[i]<=ci[2])*1.0
        mat2[seed,i] <- ci[2] - ci[1]
        
        mat3[seed,i] <- (effect_vector[i]>=ci_[1])*(effect_vector[i]<=ci_[2])*1.0
        mat4[seed,i] <- ci_[2] - ci_[1]
        
      }
      check2 <- check2 + 1
      check1 <- check1 + mat1[seed,]
      check3 <- check3 + mat2[seed,]
      
      cat("\n")
      cat("seed=",seed)
      cat("\n")
      cat("coverage=", check1/check2)
      cat("\n")
      cat("length=", check3/check2)
      cat("\n")  
      
      
      check5 <- check5 + 1
      check4 <- check4 + mat3[seed,]
      check6 <- check6 + mat4[seed,]
      
      cat("coverage=", check4/check5)
      cat("\n")
      cat("length=", check6/check5)
      cat("\n")  

    }

  
    cat("coverage_mcrts:",colMeans(mat1))
    cat("\n")
    cat("CI length_mcrts:",colMeans(mat2))
    cat("\n")
    
    cat("coverage_mcrts+:",colMeans(mat3))
    cat("\n")
    cat("CI length_mcrts+:",colMeans(mat4))
    cat("\n")
    
    list_crts[[tt]] <- c(index,weight,colMeans(mat1),colMeans(mat2))
    list_crts_f[[tt]] <- c(index,weight,colMeans(mat3),colMeans(mat4))
    tt <- tt+1
  }
}

save(list_crts, file="result_swd/syn_2_crts.RData")
save(list_crts_f, file="result_swd/syn_2_crts_f.RData")














# Mixed effects models
list_mixed = list()


tt <- 1
for (index in 1:3){
  
  
  cat("\n")
  cat("index=",index)
  
  if (index ==1){
    weight_list <- c(0,0.1)
  }else{
    weight_list <- c(0.1)
  }
  
  
  for (weight in weight_list){
    
    cat("\n")
    cat("weight=",weight)
    cat("\n")
    
    mat1 <- matrix(0, num_iters,5)
    mat2 <- matrix(0, num_iters,5)
    for (seed in 1:num_iters){
      #cat("\n")
      #cat("seed=",seed)
      data_xyz <- synthetic(num_individuals,num_steps,effect_vector,seed,weight,index,normalize=FALSE,include_x=TRUE)
      #data_xzy <- subset(data_xyz,t!=num_steps)
      
      model = lmer(y ~  (1 | id) +  x + t + z_0 + z_1 + z_2+ z_3 + z_4+z_5+z_6+z_7,  data=data_xyz, REML=FALSE)
      result<- confint(model, method = "Wald",level = 0.9)
      #print(result)
      
      lower<- result[6:10,1]
      upper<- result[6:10,2]
      
      for (i in 1:5){
        
        mat1[seed,i] <- (effect_vector[i]>=lower[i])*(effect_vector[i]<=upper[i])*1.0
        mat2[seed,i] <- upper[i]- lower[i]
      }
    }
    
    cat("coverage:",colMeans(mat1))
    cat("\n")
    cat("CI length:",colMeans(mat2))
    cat("\n")
    
    list_mixed[[tt]] <- c(index,weight,colMeans(mat1),colMeans(mat2))
    tt <- tt+1
  }
  
}


save(list_mixed, file="result_swd/syn_2_mixed.RData")










####plot the results####



library(tidyverse)

#load("results/syn_2_crts_f.RData")

#load("results/syn_2_mixed.RData")
#load("results/syn_2_crts.RData")



plot_fig <- function(dat_m,dat_crt,dat_crt_f,to_do,index,weight,i){
  
  data_ <- data.frame(lag = c(0,1,2,3,4),
                      mixed = dat_m,
                      crtsf = dat_crt_f,
                      crts = dat_crt )
  if (index == 1){
    q <- "Quadratic"
  } else if (index==2){
    q <- "Exponential"
  }  else{
    q <- "Tanh"  
  }
  
  if (i==1){
    data_$title <- paste("No Interation",sep = "")
  }else{
    data_$title <- paste(q, " Interation",sep = "")
  }
  
  
  if(to_do == "coverage") {
    p <- ggplot(data_, aes(x=lag))+
      geom_hline(aes(yintercept = 0.90))+
      xlim(-0.1, 4.1)+  
      scale_y_continuous(breaks=c(0.00, 0.30, 0.60, 0.90),limits = c(0.0,1.0))+
      geom_line(aes(y = mixed,color = "Mixed"),size=2) +
      geom_point(aes(y = mixed,color = "Mixed"),size=4)+
      geom_line(aes(y = crtsf, color =  "MCRTs+F"),size=2) +
      geom_point(aes(y = crtsf,color =  "MCRTs+F"),size=4)+
      geom_line(aes(y = crts, color =  "MCRTs+Z"),size=2) +
      geom_point(aes(y = crts,color =  "MCRTs+Z"),size=4)+
      facet_grid(. ~ title)+
      scale_color_manual(name = NULL, values = c("Mixed" = "#999999", "MCRTs+F" = "#56B4E9", "MCRTs+Z" = "#E69F00"))+
      xlab("Time lag") +
      ylab("Coverage rate")
    
    p <- p + theme_bw(base_line_size =0.5,base_rect_size = 1.0)
    p <-p + theme( legend.position="top",
                   text = element_text(size=32, face="bold"),
                   #lot.title = element_text(color="black", size=22, face="bold",hjust = 0.5),
                   axis.title.x = element_text(color="black", size=32, face="bold",vjust = -2),
                   axis.title.y = element_text(color="black", size=32, face="bold",vjust = 5),
                   plot.margin=unit(c(1.0,1.0,1.0,1.0),"cm"))
  } else if (to_do == "length"){
    p <- ggplot(data_, aes(x=lag))+
      xlim(-0.1, 4.1)+  
      #scale_y_continuous(breaks=c(0.05, 0.10, 0.15, 0.20), limits = c(0.05,0.20))+
      geom_line(aes(y = mixed,color = "Mixed"),size=2) +
      geom_point(aes(y = mixed,color = "Mixed"),size=4)+
      geom_line(aes(y = crtsf, color =  "MCRTs+F"),size=2) +
      geom_point(aes(y = crtsf,color =  "MCRTs+F"),size=4)+
      geom_line(aes(y = crts, color =  "MCRTs+Z"),size=2) +
      geom_point(aes(y = crts,color =  "MCRTs+Z"),size=4)+
      facet_grid(. ~ title)+
      scale_color_manual(name = NULL, values = c("Mixed" = "#999999", "MCRTs+F" = "#56B4E9", "MCRTs+Z" = "#E69F00"))+
      xlab("Time lag") +
      ylab("Length of CIs")
    p <- p + theme_bw(base_line_size =0.5,base_rect_size = 1.0)
    p <-p + theme( legend.position="top",
                   text = element_text(size=32, face="bold"),
                   axis.title.x = element_text(color="black", size=32, face="bold",vjust = -2),
                   axis.title.y = element_text(color="black", size=32, face="bold",vjust = 5),
                   plot.margin=unit(c(1.0,1.0,1.0,1.0),"cm"))
  } else{
    stop("Error")
  }

  ggsave(p, width=6,height=6, dpi = 300, filename =paste("result_swd/index_",index,"_weight_",weight,"_",to_do,".png"))
  return(p)
}




to_do="coverage"
d <-length(list_mixed)

for (i in 1:d){

  index<- list_mixed[[i]][1]
  weight<-list_mixed[[i]][2]
  dat_m <- list_mixed[[i]][3:7]
  dat_crt_f <-list_crts[[i]][3:7]
  dat_crt <-list_crts_f[[i]][3:7]
  p<- plot_fig(dat_m,dat_crt,dat_crt_f,to_do,index,weight,i)
}

to_do="length"

for (i in 1:d){

  index<- list_mixed[[i]][1]
  weight<-list_mixed[[i]][2]
  dat_m <- list_mixed[[i]][8:12]
  dat_crt_f <-list_crts[[i]][8:12]
  dat_crt <-list_crts_f[[i]][8:12]
  p<- plot_fig(dat_m,dat_crt,dat_crt_f,to_do,index,weight,i)
}






