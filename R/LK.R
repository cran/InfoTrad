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
  mu <- params[3];     #mu
  epsb <- params[4];   #epsilon_b
  epss <- params[5];   #epsilon_s
  
  # From now on LK stands for likelihood function.
  #Initialize
  LK_I <- c(0);
  
  for (j in 1:trad_days){
    buy_s     <- B[j];
    sell_s    <- S[j];
    
    #Compute values of interest for the log-likelihood function
    e1<- -mu-sell_s*log(1+(mu/epss))
    e2<- -mu-buy_s*log(1+(mu/epsb))
    e3<- -buy_s*log(1+(mu/epsb))-sell_s*log(1+(mu/epss))
    emax<- max(e1,e2,e3)
    
    #Split the log-likelihood in two parts and compute the relevant terms
    part1<- -epsb-epss+buy_s*log(mu+epsb)+sell_s*log(mu+epss)+emax
    part2<- log(alpha*(1-delta)*exp(e1-emax)+alpha*delta*exp(e2-emax)+(1-alpha)*exp(e3-emax))
    
    LK_I <- LK_I + (part1+part2)
  
  }   
  
  #Defining the constraints
  if ((epsb>=0) && (epss>=0) && (mu>=0) && (alpha>=0) && (delta>=0) &&  (alpha<=1) && (delta<=1)){
    LK_out <- -LK_I;
  }else{
    LK_out <- -Inf;
  }
  
  return(LK_out)
  }
}
  