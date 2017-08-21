EA_class <- function(a,d,m,eb,es,lk,pin) {
  value <- list(alpha=a, delta=d, mu=m, epsilon_b=eb, epsilon_s=es, LKval=lk, PIN=pin)
  attr(value, "class") <- "EA_class"
  value
}

print<-function(obj){
  UseMethod("print",obj)
}

# S3 Method for print
print.EA_class <- function(obj) {
  cat("Alpha:", obj$alpha, "\n")
  cat("Delta:", obj$delta, "\n")
  cat("Mu:", obj$mu, "\n")
  cat("Epsilon_b:", obj$epsilon_b, "\n")
  cat("Epsilon_s:", obj$epsilon_s, "\n")
  cat("Likelihood Value:", obj$LKval, "\n")
  cat("PIN:", obj$PIN, "\n")
}

EA<- function(data, likelihood=c("LK", "EHO")){ 
  data = as.data.frame(data)
  
  choice <- match.arg(likelihood)
  
  colnames(data)<-c('Buy', 'Sell')
  
  # Imbalance with sign
  sign_imbalance = data[,'Buy'] - data[,'Sell']
  data["sign_Imbalance"] = sign_imbalance
  
  # Absolute Imbalance
  abs_imbalance = abs(data[,'Buy'] - data[,'Sell'])
  data["abs_Imbalance"] = abs_imbalance
  
  # Complete linkage clustering
  clusters = hclust(dist(data[, "abs_Imbalance"]))
  # Divide to Event-No Event Clusters
  clusterCut = cutree(clusters, 2)
  data["Cluster"] = clusterCut
  
  t = table(clusterCut, data$Cluster)
  
  # Calculate the means of abs_imbalance for each clusters
  f_mean = mean(data[which(data[,5]==1),][,4])
  
  s_mean = mean(data[which(data[,5]==2),][,4])
  
  # Assign clustering
  means = c(f_mean, s_mean)
  event = which(means==max(means))
  no_event = which(means==min(means))
  
  # Discard no event rows from data
  event_data = data[data$Cluster==event,]
  # Cluster Event cluster into good-bad event with sign_imbalance
  clusters_E = hclust(dist(event_data[, "sign_Imbalance"]))
  clusterCut_E = cutree(clusters_E, 2, h=0)
  event_data["BadGood_Cluster"] = clusterCut_E
  
  # Calculate the mean of Buy and Sell Orders with respect to their clusters 
  bad_buy_mean = mean(event_data[event_data$BadGood_Cluster==1,"Buy"])
  good_buy_mean =mean(event_data[event_data$BadGood_Cluster==2,"Buy"])
  no_buy_mean = mean(data[data$Cluster==no_event,"Buy"])
  bad_sell_mean = mean(event_data[event_data$BadGood_Cluster==1,"Sell"])
  good_sell_mean =mean(event_data[event_data$BadGood_Cluster==2,"Sell"])
  no_sell_mean = mean(data[data$Cluster==no_event,"Sell"])
  
  # Count the number of bad,good, no news
  n_b = length(which(data$Cluster == no_event))
  n_g = length(which(event_data$BadGood_Cluster == 2))
  n_n = length(which(event_data$BadGood_Cluster == 1))
  
  # Weight
  w_b = n_b/(n_b+n_g+n_n)
  w_g = n_g/(n_b+n_g+n_n)
  w_n = n_n/(n_b+n_g+n_n)
  
  # Calculating the parameter estimates
  alpha_h = w_b+w_g
  delta_h = w_b/alpha_h
  e_b = (w_b/(w_b+w_n))*(bad_buy_mean) + (w_n/(w_b+w_n))*(no_buy_mean)
  e_s = (w_g/(w_g+w_n))*(good_sell_mean) + (w_n/(w_g+w_n))*(no_sell_mean)
  
  mu_b = max(good_buy_mean-e_b,0)
  mu_s = max(bad_sell_mean-e_s,0)
  mu = (w_g/(w_g+w_b))*(mu_b) + (w_b/(w_g+w_b))*(mu_s)
  
  
  mu = (w_g/(w_b+w_g))*mu_b + (w_b/(w_b+w_g))*mu_s
  
  par0 = c(alpha_h,delta_h,mu,e_b, e_s)
  if(choice=="LK"){
    nLL<-LK(data)    # Negative log likelihood in Lin and KE form
  }
  if(choice=="EHO"){
    nLL<-EHO(data)   # Negative log likelihood in Easleyetal form
  }
  
  #The boundary constraints
  low=c(0,0,0,0,0)
  up=c(1,1,Inf,Inf,Inf)
  
  model<-neldermead(par0,nLL,lower=low,upper=up)
  
  output<-matrix(nrow=1,ncol=7)
  
  output[1]=model$par[1] # alpha_hat
  output[2]=model$par[2] # delta_hat
  output[3]=model$par[3] # mu_hat
  output[4]=model$par[4] # eb_hat
  output[5]=model$par[5] # es_hat
  output[6]=model$value*(-1)
  
  #Calculates the probability of informed trading 
  output[7]=(output[1]*output[3])/((output[1]*output[3])+output[4]+output[5])
  ret=EA_class(output[1],output[2],output[3],output[4],output[5],output[6],output[7])
  
  return(ret)
}