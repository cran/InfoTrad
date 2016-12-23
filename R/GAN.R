GAN<- function(data, likelihood=c("LK", "EHO")){ 
  
  #Convert the data input to data.frame type
  data = as.data.frame(data)
  
  # User can choose to run LK or EHO whereas default is "LK"
  choice <- match.arg(likelihood)
  
  # Name the first column as buy, second column as sell
  colnames(data)<-c('Buy', 'Sell')
  
  imbalance = data[,'Buy'] - data[,'Sell']
  data["Imbalance"] = imbalance
  
  # Complete linkage clustering
  clusters = hclust(dist(data[, "Imbalance"]))
  clusterCut = cutree(clusters, 3)
  data["Cluster"] = clusterCut
  
  t = table(clusterCut, data$Cluster)
  
  # Calculate the means of clusters
  f_clus = t['1',]
  f_clus = t['1',][which(f_clus > 0)]
  f_mean = mean(strtoi(names(f_clus)))
  
  s_clus = t['2',]
  s_clus = t['2',][which(s_clus > 0)]
  s_mean = mean(strtoi(names(s_clus)))
  
  t_clus = t['3',]
  t_clus = t['3',][which(t_clus > 0)]
  t_mean = mean(strtoi(names(t_clus)))
  
  # Assign clustering
  means = c(f_mean, s_mean, t_mean)
  good = which(means==max(means))
  bad = which(means==min(means))
  no = which((means > min(means)) & (means < max(means)))
  
  # Calculate the mean of Buy and Sell Orders with respect to their clusters 
  bad_buy_mean = mean(data[data$Cluster==bad,"Buy"])
  good_buy_mean =mean(data[data$Cluster==good,"Buy"])
  no_buy_mean= mean(data[data$Cluster==no,"Buy"])
  
  bad_sell_mean = mean(data[data$Cluster==bad,"Sell"])
  good_sell_mean =mean(data[data$Cluster==good,"Sell"])
  no_sell_mean= mean(data[data$Cluster==no,"Sell"])
  
  # Count the number of bad,good, no news
  n_b = length(which(data$Cluster == bad))
  n_g = length(which(data$Cluster == good))
  n_n = length(which(data$Cluster == no))
  
  # Calculate the weight of each cluster
  w_b = n_b/(n_b+n_g+n_n)
  w_g = n_g/(n_b+n_g+n_n)
  w_n = n_n/(n_b+n_g+n_n)
  
  # Calculating the parameter estimates
  alpha_h = w_b+w_g
  delta_h = w_b/alpha_h
  e_b = (w_b/(w_b+w_n))*(bad_buy_mean) + (w_n/(w_b+w_n))*(no_buy_mean)
  e_s = (w_g/(w_g+w_n))*(good_sell_mean) + (w_n/(w_g+w_n))*(no_sell_mean)
  mu_b = good_buy_mean - e_b
  mu_s = bad_sell_mean - e_s
  mu = (w_g/(w_b+w_g))*mu_b + (w_b/(w_b+w_g))*mu_s
  
  par0 = c(alpha_h,delta_h,mu,e_b, e_s)
  if(choice=="LK"){
    nLL<-LK(data)    # Negative log likelihood in Lin and KE form
  }
  if(choice=="EHO"){
    nLL<-EHO(data)   # Negative log likelihood in Easleyetal form
  }
  
  run<-try(model<-optim(par0, nLL, gr = NULL, method = c("Nelder-Mead"), hessian = FALSE)) # Nelder Mean
  if (class(run)=="try-error"){
    print("Error in optimization procedure")
  }
  
  output<-matrix(nrow=1,ncol=7)
  
  output[1]=model$par[1] # alpha_hat
  output[2]=model$par[2] # delta_hat
  output[3]=model$par[3] # mu_hat
  output[4]=model$par[4] # eb_hat
  output[5]=model$par[5] # es_hat
  output[6]=model$value*(-1)
  
  #Calculates the probability of informed trading 
  output[7]=(output[1]*output[3])/((output[1]*output[3])+output[4]+output[5])
  
  # Naming the columns and rows
  colnames(output)<-c("alpha","delta","mu", "epsilon_b","epsilon_s","LikVal","PIN")
  rownames(output)<-c("Output")
  return (output)
}