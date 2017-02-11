YZ<-function(data, likelihood=c("LK","EHO")){
  
  # User can choose to run LK or EHO whereas default is "LK"
  choice <- match.arg(likelihood)
  
  #Separating buy values from sell values  
  Bt<-data[,1]
  St<-data[,2]

  #Replace NA's with zero to allow the existence of zero order days.
  Bt[is.na(Bt)]=0
  St[is.na(St)]=0
  
  #Defining variables for the mean buy and sell orders and the max amount of order
  mean_b=mean(Bt)
  mean_s=mean(St)
  maxim=max(max(Bt),max(St))

  #Setting the input parameter panel
  alpha_in<-c(0.1, 0.3, 0.5, 0.7, 0.9)
  delta_in<-c(0.1, 0.3, 0.5, 0.7, 0.9)
  gamma_in<-c(0.1, 0.3, 0.5, 0.7, 0.9)

  #Creating a matrix for estimates driven through all combinations of initial values
  #Output matrix contains estimates of alpha, delta, mu, eb, es and the likelihood value for all initial combinations.
  output=matrix(nrow=125,ncol=6)
  
  #A dummy is created for moving in out matrix. 
  dum=1
  for (i in (1:5)){
    for( j in (1:5)){
      for(k in (1:5)){
      
        alpha=alpha_in[i]
      
        delta=delta_in[j]
      
        eb=gamma_in[k]*mean_b
      
        mu=(mean_b-eb)/(alpha*(1-delta))
      
        es=mean_s-(alpha*delta*mu)
        
        #YZ indicates that es cannot be negative
        #In addition this part forces to provide initial values for mu that are less maximum value for order.
        if (es<=0 | mu>maxim){
          dum=dum+1
          next
        }
      
        par0<-c(alpha,delta,mu,eb,es)
        if(choice=="LK"){
          nLL<-LK(data) #Negative log likelihood in Lin and KE form
        }
        if(choice=="EHO"){
          nLL<-EHO(data) #Negative log likelihood in Easleyetal form
        }
  
        ##Optim package###
          run<-try(model<-optim(par0, nLL, gr = NULL, method = c("Nelder-Mead"), hessian = FALSE))#Nelder Mean
          if (class(run)=="try-error"){
            dum=dum+1
            next
          }
          output[dum,1]=model$par[1]# alpha_hat
          output[dum,2]=model$par[2]# delta_hat
          output[dum,3]=model$par[3]# mu_hat
          output[dum,4]=model$par[4]# eb_hat
          output[dum,5]=model$par[5]# es_hat
          output[dum,6]=model$value*(-1)

        dum=dum+1
      }
    }
  }
  
  #Sorts the output matrix such that the last term will provide the estimates that maximize LK
  output_s<-output[order(output[,6],na.last=FALSE,decreasing=FALSE),]
  
  # Give column names
  colnames(output_s) = c('Alpha','Delta','Mu','Eb','Es','LKVal')
  
  # Remove rows with Alpha=0
  output_s = output_s[output_s[,"Alpha"]!=0,]
  
  #Final output file 
  fin_output<-matrix(nrow=1,ncol=7)
  #Collects the parameter estimates
  fin_output[1]=output_s[nrow(output_s),1] #alpha_h
  fin_output[2]=output_s[nrow(output_s),2] #delta_h
  fin_output[3]=output_s[nrow(output_s),3] #mu_h
  fin_output[4]=output_s[nrow(output_s),4] #eb_h
  fin_output[5]=output_s[nrow(output_s),5] #es_h
  fin_output[6]=output_s[nrow(output_s),6] #Fval
  
  #Calculates the probability of informed trading 
  fin_output[7]=(fin_output[1]*fin_output[3])/((fin_output[1]*fin_output[3])+fin_output[3]+fin_output[4])
  
  # Naming the columns and rows
  colnames(fin_output)<-c("alpha", "delta", "mu", "epsilon_b", "epsilon_s", "Likval", "PIN")
  rownames(fin_output)<-c("Output")
  return(fin_output)
}

