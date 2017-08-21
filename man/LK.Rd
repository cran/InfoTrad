\name{LK}
\alias{LK}
\title{Likelihood factorization of Lin and Ke (2011) - LK Factorization}
\description{The function calculates the likelihood factorization of Lin and Ke (2011) and computes paramaters for estimation of PIN value.}
\usage{LK(data, fixed = c(FALSE, FALSE, FALSE, FALSE, FALSE))}
\arguments{
  \item{data}{Data frame with 2 variables}
  \item{fixed}{Initial values for parameters in the following order: alpha, delta, mu, epsilon_b, epsilon_s}
}

\details{In order to use LK's return in optimization functions, please \strong{omit} second argument. With this way, LK will return a function instead of a value. Moreover, argument for data must be a data frame with 2 columns that contain numbers. Not any other type.}

\section{Warning }{This function does not handle NA values. Therefore the datasets should not contain any missing value}
\value{
  \item{LK_out}{Returns an optim() object including parameter estimates for the likelihood factorization of Lin and Ke (2011)}
}
\references{
Lin, H.W.W. and Ke, W.C. A computing bias in estimating the probability of informed trading. Journal of Financial Markets, 14(4), pp.625-640, 2011.
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
  
  # Call LK function
  LK_out = LK(data)
  model = optim(par0, LK_out, gr = NULL, method = c("Nelder-Mead"), hessian = FALSE)
  
  ## Parameter Estimates
  model$par[1] # Estimate for alpha
  # [1] 0.480277
  model$par[2] # Estimate for delta
  # [1] 0.830850
  model$par[3] # Estimate for mu
  # [1] 315.259805
  model$par[4] # Estimate for eb
  # [1] 296.862318
  model$par[5] # Estimate for es
  # [1] 434.3046
  
  ## Estimate for PIN 
  (model$par[1]*model$par[3])/((model$par[1]*model$par[3])+model$par[4]+model$par[5])
  # [1] 0.178391
  ####
}

