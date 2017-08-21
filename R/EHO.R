EHO<- function(data, fixed=c(FALSE,FALSE,FALSE,FALSE,FALSE)){
  params<-fixed
  function(p){
    
  params[!fixed]<-p
  #Assign the number of buy-sell orders
  B<-data[,1]
  S<-data[,2]
  
  #Trading days
  trad_days <- length(B);
  
  #Initialize parameters values
  alpha <- params[1];  #alpha 
  delta  <- params[2]; #delta
  mu <- params[3]; #mu
  epsb <- params[4]; #e_b
  epss <- params[5]; #e_s
  
  #Initialize
  LK_I=0;
  
  for (j in 1:trad_days){
    buy_s     <- B[j];
    sell_s    <- S[j];
    
    #Compute values of interest for the log-likelihood function
    M<- min(buy_s,sell_s)+max(buy_s,sell_s)/2
    Xs<- epss/(mu+epss)
    Xb<- epsb/(mu+epsb)
      
    #Split the log-likelihood in two parts and compute the relevant terms
    part1<- buy_s*log(mu+epsb)+sell_s*log(mu+epss)+M*log(Xb)+M*log(Xs)-(epsb+epss)
    part2<- log(alpha*delta*exp(-mu)*Xb^(buy_s-M)*Xs^(-M)+alpha*(1-delta)*exp(-mu)*Xb^(-M)*Xs^(sell_s-M)+(1-alpha)*Xb^(buy_s-M)*Xs^(sell_s-M))
    
    if (is.nan(part2)){
      part2=0
    }
    LK_I <- LK_I + (part1+part2)
  
  }   
  return(-LK_I)
  }
}
  