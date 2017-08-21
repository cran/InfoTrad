\name{EHO}
\alias{EHO}
\title{Likelihood factorization of Easley et. al. (2010) - EHO Factorization}
\description{The function calculates the likelihood factorization of Easley et. al. (2010) and computes paramaters for estimation of PIN value.}
\usage{EHO(data, fixed = c(FALSE, FALSE, FALSE, FALSE, FALSE))}
\arguments{
  \item{data}{Data frame with 2 variables}
  \item{fixed}{Initial values for parameters in the following order: alpha, delta, mu, epsilon_b, epsilon_s}
}
\details{In order to use EHO's return in optimization functions, please \strong{omit} second argument. With this way, EHO will return a function instead of a value. Moreover, argument for data must be a data frame with 2 columns that contain numbers. Not any other type.}
\section{Warning}{This function does not handle NA values. Therefore the datasets should not contain any missing values.}
\value{
  \item{LK_out}{Returns an optim() object including parameter estimates for the likelihood factorization of Easley et. al. (2010)}
}
\references{
Easley, D., Hvidkjaer, S., & O'Hara, M. Factoring information into returns. Journal of Financial and Quantitative Analysis, 45(2):293-309,2010.
}
\author{Duygu Celik and Murat Tinic}

\examples{
  # Sample Data
  #   Buy Sell 
  #1  350  382  
  #2  250  500  
  #3  500  463  
  #4  552  550  
  #5  163  200  
  #6  345  323  
  #7  847  456  
  #8  923  342  
  #9  123  578  
  #10 349  455 
  
  Buy<-c(350,250,500,552,163,345,847,923,123,349)
  Sell<-c(382,500,463,550,200,323,456,342,578,455)
  data=cbind(Buy,Sell)

  # Initial parameter values
  # par0 = (alpha, delta, mu, epsilon_b, epsilon_s)
  par0 = c(0.5,0.5,300,400,500)

  # Call EHO function
  EHO_out = EHO(data)
  model = optim(par0, EHO_out, gr = NULL, method = c("Nelder-Mead"), hessian = FALSE)

  ## Parameter Estimates
  model$par[1] # Estimate for alpha
  # [1] 0.9111102
  model$par[2] # Estimate for delta
  #[1] 0.0001231429
  model$par[3] # Estimate for mu
  # [1] 417.1497
  model$par[4] # Estimate for eb
  # [1] 336.075
  model$par[5] # Estimate for es
  # [1] 466.2539
  
  ## Estimate for PIN
  (model$par[1]*model$par[3])/((model$par[1]*model$par[3])+model$par[4]+model$par[5])
  # [1] 0.3214394
  ####

}

