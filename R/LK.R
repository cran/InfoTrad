LK <- function(data, fixed=c(FALSE,FALSE,FALSE,FALSE,FALSE)){
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
    e1<- (-1)*(mu+buy_s*log(1+(mu/epsb)))
    e2<- (-1)*(mu+sell_s*log(1+(mu/epss)))
    e3<- (-1)*(buy_s*log(1+(mu/epsb))+sell_s*log(1+(mu/epss)))
    emax<- max(e1,e2,e3)
    
    #Split the log-likelihood in two parts and compute the relevant terms
    part1<- buy_s*log(mu+epsb)+sell_s*log(mu+epss)-(epsb+epss)+emax
    part2<- log(alpha*delta*exp(e1-emax)+alpha*(1-delta)*exp(e2-emax)+(1-alpha)*exp(e3-emax))
    
    LK_I <- LK_I + (part1+part2)
  
  }   
  return(-LK_I)
  }
}
  