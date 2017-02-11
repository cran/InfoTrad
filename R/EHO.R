## This is a comment

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
  LK_I <- c(0);
  
  for (j in 1:trad_days){
    buy_s     <- B[j];
    sell_s    <- S[j];
    
    #Compute values of interest for the log-likelihood function
    M<- min(buy_s,sell_s)+max(buy_s,sell_s)/2
    Xs<- epss/(mu+epss)
    Xb<- epsb/(mu+epsb)
      
    #Split the log-likelihood in two parts and compute the relevant terms
    part1<- -(epsb+epss)+M*log(Xb+Xs)+buy_s*log(mu+epsb)+sell_s*log(mu+epss)
    part2<- log(alpha*(1-delta)*exp(-mu)*Xs^(sell_s-M)*Xb^(-M)+alpha*delta*exp(-mu)*Xb^(buy_s-M)*Xs^(-M)+(1-alpha)*Xs^(sell_s-M)*Xb^(buy_s-M))
    
    LK_I <- LK_I + (part1+part2)
  
  }   
  
  if ((epsb>=0) && (epss>=0) && (mu>=0) && (alpha>=0) && (delta>=0) &&  (alpha<=1) && (delta<=1)){
    LK_out <- -LK_I;
  }else{
    LK_out <- -Inf;
  }
  
  return(LK_out)
  }
}
  