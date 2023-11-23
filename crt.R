# Load required libraries
library(metap)
library(mutoss)
library(poolr)
library(reshape)
library(igraph)
library(adagio)


remove_singleton <- function(set,G_){
  
  # remove the units with 0 neighbors
  
  set[rowSums(G_[set,])>1.001]
}


pval = function(y, tol=1e-12) {
  
  # Compute the p-value given the test statistics under different assignments
  
  size = length(y)
  M = sum(abs(y[2:size]-y[1])<tol)
  p_val = (sum(y[2:size]>y[1]) + runif(1)*(1+M)) / size
  p_val
}


count_perm <- function(treated_,control_){
  
  #Count the number of permutations
  
  mu_treated <- treated_$y
  mu_treated <- mu_treated[!is.na(mu_treated)]
  mu_control <- control_$y
  mu_control <- mu_control[!is.na(mu_control)]
  
  num1 <- length(mu_treated)
  num2 <- length(mu_control)
  num <- num1 + num2
  group <- combn(1:num, num1)
  NC <- NCOL(group)
  return(NC)
}


example_space = function(num_units) {
  
  # Generate synthetic data from mixtures of Gaussian distribution
  
  X_1 = c(rnorm(num_units * 0.5, sd = 0.1) + 0.5, rnorm(num_units * 0.3, sd = 0.075) + 0.3, rnorm(num_units * 0.2, sd = 0.075) +  0.4)
  X_2 = c(rnorm(num_units * 0.5, sd = 0.1) + 0.5, rnorm(num_units *  0.3, sd = 0.075) + 0.6, rnorm(num_units * 0.2, sd = 0.075) +  0.2)
  
  X_1 = X_1*100
  X_2 = X_2*100
  thepoints = cbind(X_1, X_2)
  Dmat = array(0, c(nrow(thepoints), nrow(thepoints)))
  for (ii in 1:nrow(thepoints)) {
    for (jj in 1:nrow(thepoints)) {
      Dmat[ii, jj] = sqrt(sum((thepoints[ii, ] - thepoints[jj,])^2))
    }
  }
  return(list(D = Dmat, thepoints = thepoints))
}



setup <- function(seed,N,prob,radius_gap,radius_list,signal_list,signal=1){
  
  # Setup the synthetic interference experiments
  
  out = NULL
  while(is.null(out)){
    
    #Skip a few imbalanced runs that have no control units
    tryCatch(
      {
        withTimeout({
          set.seed(seed)
          graph <- example_space(N)
          
          # Generate the distance matrix
          C <- graph$D
          
          # Randomize the assignment
          Z = sample(c(rep(1, round(N*prob) ), rep(0, N-round(N*prob))))
          
          # Generate treated or control outcomes
          Y =  Z*4 + rnorm(N)
          
          # Add the spillover effects
          for (j in 1:length(radius_list)){
            radius = radius_list[j]
            if (signal==1){
              Y = Y + ((C>(radius-radius_gap))* (C<=radius))%*%Z*signal_list[j] 
            }
          }
          

          randoms = list()
          # Run all the CRTs
          set_sum = c()
          for (j in 1:length(radius_list)){
            
            radius = radius_list[j]
            
            # Create the interference matrix
            G0 <- (C<=radius)*1.0 
            
            # Remove the units that are randomized previously
            set_left = setdiff(1:N,set_sum)
            
            # Find the cover
            sol <- setcover(G0[set_left,set_left])
            set = set_left[sol$set]
            
            # Remove the covering points which only covers themselves
            set_ <- remove_singleton(set,G0)
            randoms[[j]] = set_
            set_sum = c(set_sum,set_)
            
          }
          
          
          
          G = list()
          J = list()
          I = list()
          K = list()
          
          
          for (j in 1:length(radius_list)){
            
            radius = radius_list[j]
            

            
            # Define the focal units
            G[[j]] = (C>(radius-radius_gap))*(C<=radius) * 1; diag(G[[j]]) = 0
            
            
            set_sum = c()
            for (j_1 in j:length(radius_list)){
              set_sum = c(set_sum,randoms[[j_1]])
            }
            
            set_nn = which(colSums(G[[j]][randoms[[j]],])>0.001 & colSums((C<=radius-radius_gap)[set_sum,])==0)
            
            set_nn = setdiff(set_nn, set_sum)
            #save the setup
            J[[j]] <- set_nn 
            I[[j]] <- randoms[[j]]
            K[[j]] <- set_sum
            
          }
          
          out = list(G = G,  Y = Y, Z = Z, J = J, I=I, K= K)
          
        },timeout=6000)}, error=function(error_message) {
          message("next seed.")
        }
    )
    
    #adjust random seed
    seed = seed +1234
  }
  

  return(out)
}


subdesign_matrix <- function(Z,num_assignments,II){
  
  # Permute part of the assignment matrix 
  
  # The first copy is not permuted
  Z_status <- c(Z)
  
  # Generate permutations
  perm_list <- 2:num_assignments
  set.seed(NULL)
  for (seed2 in perm_list){
    Z_seed <- Z
    Z_seed[II] <- sample(Z[II])
    Z_status <- cbind(Z_status,Z_seed) 
  }
  set.seed(NULL)
  
  return(Z_status)
}


run_seq_test <- function(Z_,Y_,I, J, G, alpha,sign=1.0,threshold=NULL){
  
  # Compute the randomization p-value
  
  minimum_size <- 1/alpha
  p_value <- 1.0
  size <- 0
  Y_s <- c()
  
  # Extract the relevant part of the inference matrix
  G_JI = G[J,I]
  num_treated <- sum(Z_[I,1]==1)
  num_control <- sum(Z_[I,1]==0)
  b = choose(length(I),num_treated)
  Z_ <- data.frame(Z_)
  if(b>=minimum_size){
    
    for (t in 1:ncol(data.frame(Z_))){
      
      #Compute the difference-in-means statistics
      
      Z_t = Z_[,t]
      status = (G_JI%*%Z_t[I]>0.001)*1.0
      test <- mean(Y_[J[status==1]]) - mean(Y_[J[status==0]])
      test <- test*sign
      
      if (is.na(test)==TRUE){
        test = 1000000000
      }
      
      Y_s <- c(Y_s,test)
    } 
    
    p_value <- pval(Y_s)
    
  }
  
  # Compute the randomization distribution for the forest example
  #if (is.null(threshold)==FALSE){
  #  png(paste0("result_interference/forest_hist_", threshold, ".png", sep=""), width=8, height=8, units="in", res=1200)
  #  par(mar = c(5, 5, 3, 4.5))
  #  hist(Y_s, col="lightsalmon1", xlab='Test statistics',main=paste("P = ",round(p_value,3)), cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2,cex=2)
  #  abline(v=Y_s[1],
  #         col="tomato3",
  #         lty=2,
  #         lwd=5.0)
  #  dev.off()
  #}
  return(p_value)
  
}


run_seq_ci <- function(Z_, Y_,I, J, G, alpha,threshold,sign=1.0){
  
  # Compute the randomization confidence intervals
  
  minimum_size <- 1/alpha
  p_value <- 1.0
  size <- 0
  
  # Extract the relevant part of the inference matrix
  G_JI = G[J,I]
  num_treated <- sum(Z_[I,1]==1)
  num_control <- sum(Z_[I,1]==0)
  
  p_list = c()
  
  # Domain for grid search
  beta_list = seq(-100,100,0.1)
  Z_ <- data.frame(Z_)
  if(choose(length(I),num_treated)>=minimum_size){
    
    Z_1 = Z_[,1]

    for (beta in beta_list){
      
      
      # Adjust the outcomes
      add_effect = G_JI%*%Z_1[I]*beta 
      Y_beta = Y_
      Y_beta[J] = Y_beta[J] - add_effect
      Y_s <- c()
      
      # Loop over the permutations
      for (t in 1:ncol(data.frame(Z_))){
        Z_t = Z_[,t]
        status = (G_JI%*%Z_t[I]>0.001)*1.0
        test <- mean(Y_beta[J[status==1]]) - mean(Y_beta[J[status==0]])
        test <- test*sign
        
        if (is.nan(test)) {
          test = 1000000000
        }
        
        Y_s <- c(Y_s,test)
      }
      p_value <- min(max(mean(Y_s >= Y_s[1]),0.0000000001),0.9999999999)
      p_list <- c(p_list,p_value)
    } 
  }
  
  # Invert the p-values into confidence intervals
  u = beta_list[which.min(abs(p_list-alpha/2))]
  l = beta_list[which.min(abs(p_list-1+alpha/2))]

  return(c(l,u))
}


data_swd <- function(data,lag,exact=FALSE,sample_split=TRUE){
  
  num_steps <- max(data$t)
  data_list <-list()
  if (sample_split == FALSE){
    for (j in 1:(num_steps-lag-1)){
      treated_ <- subset(data,s==j & t == (j+lag))
      control_ <- subset(data,s>(j+lag) & t == (j+lag))
      data_list[[j]] <- list(treated_,control_)
    }
  }else if(sample_split == TRUE){
    
    num_groups <- min(lag+1,num_steps-lag-1)
    
    if (num_groups<=0){
      stop("Error: Please enter a time lag < number of time steps - 1")
    }
    
    group_list = list()
    for (t in 1:num_groups){
      step <- t
      group <- step 
      while ((step + lag + 1) <= num_steps){
        step <- step + lag +1 
        group <- c(group,step)
      }
      
    
      group_list[[length(group_list)+1]] <- group

    }
    
    for (group in group_list){

      group_size <- length(group)
      for (j in 1:(group_size-1)){
        treated_ <- subset(data,s==group[j] & t == (group[j]+lag))
        control_ <- subset(data,s %in% group[(j+1):group_size] & t == (group[j]+lag))
        
        if (exact==FALSE){
          data_list[[length(data_list)+1]] <- list(treated_,control_)
        }else if (exact==TRUE){
          num_perm <- count_perm(treated_,control_)
          if (num_perm> 1/alpha){
            data_list[[length(data_list)+1]] <- list(treated_,control_)
          }
        } else {
          stop("Error: R stops due to an incorrect setup")
        }
      }
    }  
  } else{
    stop("Error: R stops due to an incorrect setup")
  }
  return(data_list)
}
































































crt_swd <- function(data_list,type,beta=0,exact=FALSE,b=1000,n_beta=1000,const=3,print=FALSE){
  
  # Run a sequence of CRTs for stepped-wedge design
 
  d <- length(data_list)
  
  if (print==TRUE){
    cat("\n")
    cat("Number of tests:", d)
    cat("\n")
  }

  if (type=='p_values'){
    p_list <-NULL
    for (j in 1:d){
      data_i <- data_list[[j]]
      if (print==TRUE){
        cat("Permutation between steps", sort(unique(data_i[[1]]$s)),"and", sort(unique(data_i[[2]]$s)))
        cat("\n")
      }
      p_value <- crt_p(data_i, beta, exact, b)
      p_list <- c(p_list,p_value)
    }
    

  }else if (type=='ci'){
    p_matrix_1 <- matrix(0, d, n_beta)
    p_matrix_2 <- matrix(0, d, n_beta)
    for (j in 1:d){
      data_i <- data_list[[j]]
      if (print==TRUE){
        cat("Permutation between steps", sort(unique(data_i[[1]]$s)),"and", sort(unique(data_i[[2]]$s)))
        cat("\n")
      }
      out <- crt_ci(data_i, const, n_beta, exact, b)
      p_matrix_1[j,] <- out[[1]]
      p_matrix_2[j,] <- out[[2]]
    }
    
    p_list <- list(p_matrix_1,p_matrix_2)
  
  } else {
    stop("Error: R stops due to an incorrect setup")
  }

 return(p_list)  
}



global_test <- function(data_list,p_list,method,alpha = 0.1){
  
  # P-value combiners in stepped-wedge design
  
  if (length(p_list)==1){
    
    reject <- (p_list <= alpha)
    
  }else{
    
    if (method=='holms'){
      result<- holm(p_list, alpha,silent=TRUE)
      p_value <- min(p_list)
      reject<- p_value <= alpha/length(p_list) 
    } else if(method=='mcrts'){
      p_value<- sumlog(p_list,log.p = FALSE)['p'][[1]]
      reject<- p_value <= alpha
    } else if(method=='mcrts+'){
      w <- optimal_weight(data_list)
      T_value <- sum(w * qnorm(p_list,mean=0,sd=1))
      p_value <- pnorm(T_value, 0, 1)
      reject<- p_value <= alpha
    } else{
      stop("Error: R stops due to an incorrect setup")
    }
    
  }
  return(reject)
}



global_ci <- function(data_list,p_list,method,alpha=0.1,const=3){
  
  # Invert combined p-values into confidence intervals for lagged effects
  
  p_matrix_1 <- p_list[[1]]
  p_matrix_2 <- p_list[[2]]
  significance <- 0
  n_beta <- NCOL(p_matrix_1)
  alpha_0 <- alpha/2+0.05

  while (significance <(1-alpha) & alpha_0>= 0.00){
    list1 <- NULL
    list1_p <- NULL
    list2 <- NULL
    list2_p <- NULL
    
    
    if (NROW(p_matrix_1) >1){
      if (method=='mcrts+'){
        w<- optimal_weight(data_list)
        for (beta in 1:n_beta){
          add <- 2*const*beta/n_beta -const
          
          T_1 <- sum(w * qnorm(p_matrix_1[,beta],mean=0,sd=1))
          p <- pnorm(T_1, 0, 1)
          if (p >=alpha_0){
            list1 <- c(list1,add)
            list1_p <- c(list1_p, p)
          }
          
          T_2 <- sum(w * qnorm(p_matrix_2[,beta],mean=0,sd=1))
          p <- pnorm(T_2, 0, 1)
          if (p >=alpha_0){
            list2 <- c(list2,add)
            list2_p <- c(list2_p,p)
          }
        }
      } else if(method=='mcrts'){
        
        for (beta in 1:n_beta){
          add <- 2*const*beta/n_beta -const
          p <- sumlog(p_matrix_1[,beta],log.p = FALSE)['p'][[1]] 
          if (p >=alpha_0){
            list1 <- c(list1,add)
            list1_p <- c(list1_p,p)
          }
          
          p <- sumlog(p_matrix_2[,beta],log.p = FALSE)['p'][[1]] 
          if (p >=alpha_0){
            list2 <- c(list2,add)
            list2_p <- c(list2_p,p)
          }
        }
      } else {
        stop("Error: R stops due to an incorrect setup")
      }
    }else if(NROW(p_matrix_1) ==1){
      for (beta in 1:n_beta){
        add <- 2*const*beta/n_beta -const
        p <- p_matrix_1[1,beta]
        if (p >=alpha_0){
          list1 <- c(list1,add)
          list1_p <- c(list1_p,p)
          #cat("beta_max_u:", add,res1)
          #cat("\n")
        }
        p <- p_matrix_2[1,beta]
        if (p >=alpha_0){
          list2 <- c(list2,add)
          list2_p <- c(list2_p,p)
          #cat("beta_min_l:", add,res2)
          #cat("\n")
        }
      } 
    }else {
      stop("Error: R stops due to an incorrect setup")
    }
    
    u <- max(list1)
    l <- min(list2)
    
    ci_levels <- c(list2_p[list2==l],1-list1_p[list1==u])
    ci<-c(l,u)
    
    significance <- ci_levels[2]-ci_levels[1]
    alpha_0 <- alpha_0 - 0.0025
  }
  
  
  #cat("Significance Level:", significance)
  
  return(ci)
}



crt_p <- function(data_i, beta=0, exact=FALSE,b=1000){
  
  # Compute the p-value for a CRT
  
  mu_treated <- data_i[[1]]$y
  mu_treated <- mu_treated[!is.na(mu_treated)]
  mu_control <- data_i[[2]]$y
  mu_control <- mu_control[!is.na(mu_control)]
  
  num1 <- length(mu_treated)
  num2 <- length(mu_control)
  num <- num1 + num2
  
  rr <- perm_matrix(num1,num2, exact, b)
  x <- mu_treated - beta
  y <- mu_control
  obs_t <- mean(x)-mean(y)
  xy <- c(x,y)
  counter_t <- rr %*% xy
  p_value <- pval(c(obs_t,counter_t))
  #mean(counter_t >= obs_t)
  p_value <- min(max(p_value,0.0000000001),0.9999999999)
  return(p_value)
}



crt_ci <- function(data_i,const,n_beta,exact,b){
  
  # Compute the randomization confidence interval for lagged effects
  
  #cat("Permutation between steps", unique(data_i[[1]]$s),"and", unique(data_i[[2]]$s))
  #cat("\n")
  mu_treated <- data_i[[1]]$y
  mu_treated <-mu_treated[!is.na(mu_treated)]
  mu_control <- data_i[[2]]$y
  mu_control <-mu_control[!is.na(mu_control)]
  
  vec1 <- NULL
  vec2 <- NULL
  num1 <- length(mu_treated)
  num2 <- length(mu_control)

  rr <- perm_matrix(num1,num2, exact, b)
  
  for (beta in 1:n_beta){
    add <- 2*const*beta/n_beta -const
    
    y_1 <- mu_treated - add
    y_2 <- mu_control
    y_12 <- c(y_1,y_2)
   
    obs_t <- mean(y_1)-mean(y_2)
    counter_t <- rr %*% y_12
    
    #hist(counter_t)
    #abline(v = obs_t, col = "blue", lwd = 2)
    p_value_1 <- mean(counter_t <= obs_t) 
    p_value_2 <- mean(counter_t >= obs_t) 
    vec1 <- c(vec1, min(max(p_value_1,0.0000000001),0.9999999999)) 
    vec2 <- c(vec2, min(max(p_value_2,0.0000000001),0.9999999999))
    
  }
  return(list(vec1,vec2))
}




perm_matrix <- function(num1,num2, exact,b){
  
  #Run a (MC or exact) permutation test
  
  num <- num1 + num2
  if (exact == TRUE){
    group <- combn(1:num, num1)
    NC <- NCOL(group)
    #cat("Number of permutations:", NC)
    cat("\n")
    rr <- matrix(-1/num2, nrow=NC, ncol=num)
    for (j in 1:NC){
      rr[j,group[,j]] <- 1/num
    }
  }else if (exact == FALSE){
    var <- c(rep(1/num1, num1), rep(-1/num2, num2))
    rr <- matrix(0, nrow=b, ncol=num)
    
    for (j in 1:b){

      rr[j,] <- sample(var, size = num, replace=FALSE)
    } 

  }else{
    stop("Error: R stops due to an incorrect setup")
  }
  return(rr)
}



optimal_weight <- function(data_list){
  
  # Compute the weights in the Z-score p-value combiner
  
  r_list = NULL
  m_total <- 0
  for (data in data_list){
    dat_1 <- data[[1]]$y
    dat_1 <- dat_1[!is.na(dat_1)]
    dat_0 <- data[[2]]$y
    dat_0 <- dat_0[!is.na(dat_0)]

    m_1 <- NROW(dat_1)  
    m_0 <- NROW(dat_0)  
    m_ <- m_0 + m_1
    m_total <- m_total + m_
     
  }
  
  for (data in data_list){
    dat_1 <- data[[1]]$y
    dat_1 <- dat_1[!is.na(dat_1)]
    dat_0 <- data[[2]]$y
    dat_0 <- dat_0[!is.na(dat_0)]
    
    v_1 <- var(dat_1)
    v_0 <- var(dat_0)
    
    m_1 <- NROW(dat_1)  
    m_0 <- NROW(dat_0)  
    r_ <-  1/(v_1*m_total/m_0 + v_0*m_total/m_1)
    r_list <- c(r_list,r_)
  }
  w <- sqrt(r_list/sum(r_list))
  return(w)
}




synthetic <- function(num_individuals,num_steps,effect_vector,seed,weight=0,index=1,normalize=TRUE,include_x=TRUE){
  
  
  #Synthetic simulations for stepped-wedge design
  
  
  #set.seed(NULL) 
  data_matrix <- matrix(0,nrow = num_individuals,ncol = num_steps + 1)
  
  x <- rnorm(num_individuals, 0, 0.25)
  
  random_effect<-rnorm(num_individuals, 0, 0.25)
  
  #set.seed(seed)
  subsample <- round(num_individuals/num_steps)
  design <- c(rep(subsample, num_steps-1),num_individuals-subsample*(num_steps-1))
  s <- NULL
  for (t in 1:num_steps){
    s <- c(s,rep(t,design[t]))
  }
  s <- sample(s, size = num_individuals, replace=FALSE)
  
 
  #beta <- rnorm(dim, 0, 1)
  #beta_t <- 0.1

  
  for (t in 0:num_steps){
    
    if (include_x==TRUE){
      factor_t <- 0.5*(x  + t)  
    }else{
      factor_t<- 0.5*t
    }
    if (index==1){
      add <- 1*factor_t*factor_t 
    }else if (index==2){
      add <- 2*exp(factor_t/2)
    }else if (index==3){
      add <- 5*tanh(factor_t)
    }else{
      stop("Error: R stops due to an incorrect setup")
    }
    
    noise <- rnorm(num_individuals, 0, 0.1)
    
    data_matrix[,t+1] <- (1- weight)*factor_t + weight*add+ random_effect+noise 
    
    for (i in 1:num_individuals){
      gap <- t-s[i]
      if (gap>=0){
        data_matrix[i,t+1] <- data_matrix[i,t+1] + effect_vector[gap+1]
      }
    }
  }
  
  
  x_data <- data.frame(x)
  x_data$id <- seq.int(num_individuals)
  
  s_data <- data.frame(s)
  s_data$id <- seq.int(num_individuals)
  
  mydata <- data.frame(data_matrix)
  colnames(mydata) <- c(0:num_steps)
  mydata$id <- seq.int(num_individuals)
  
  mdata <- melt(mydata, id=c("id"))
  colnames(mdata)[2] <- "t"
  colnames(mdata)[3] <- "y"
  mdata <- mdata[order(mdata$id),]
  
  total <- merge(s_data,x_data,by="id")
  final <- merge(mdata,total,by="id")
  
  final$t <- as.numeric(final$t) -1 
  lag <- final$t - final$s
  num <- NROW(final)
  z_matrix <-  matrix(0,nrow = num,ncol = num_steps)
  for (j in 1:num){
    if (lag[j]>=0){
      z_matrix[j,lag[j]+1] <- 1
    }
  }
  z_data <- data.frame(z_matrix)
  for(i in 0:(NCOL(z_data)-1)){                   
    colnames(z_data)[i+1] <- paste0("z_", i)
  }
  data_xyz <- cbind(final,z_data) 
  
  if (normalize ==TRUE){
    for (i in 1:num_individuals){
      rows <- data_xyz$id
      counter <- subset(data_xyz,id ==i & t==0)
      data_xyz[rows==i,]$y <- data_xyz[rows==i,]$y - counter$y
    }  
  }
  
  return(data_xyz)
}


