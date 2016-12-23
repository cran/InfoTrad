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
  #    Buy    Sell
  # 1  350    382
  # 2  150    500
  # 3  500    463
  # 4  552    550
  # 5  163    200
  
  Buy=c(350,150,500,552,163)
  Sell=c(382,500,463,550,200)
  data=cbind(Buy,Sell)
  
  # Initial parameter values
  # par0 = (alpha, delta, mu, epsilon_b, epsilon_s)
  par0 = c(0.5,0.5,300,400,500)
  
  # Call LK function
  LK_out = LK(data)
  model = optim(par0, LK_out, gr = NULL, method = c("Nelder-Mead"), hessian = FALSE)
  
  ## The structure of the model output ##
  model
  
  # $par
  # [1] 6.339028e-01 3.358580e-06 4.302504e+02 1.348861e+02 4.343046e+02
  # $value
  # [1] -19004.07
  
  # $counts
  # function gradient 
  # 392       NA 
  
  # $convergence
  # [1] 0
  
  # $message
  # NULL
  ##
  
  ## Parameter Estimates
  model$par[1] # Estimate for alpha
  # [1] 0.6339028
  model$par[2] # Estimate for delta
  # [1] 3.35858e-06
  model$par[3] # Estimate for mu
  # [1] 430.2504
  model$par[4] # Estimate for eb
  # [1] 134.8861
  model$par[5] # Estimate for es
  # [1] 434.3046
  
  ## Estimate for PIN 
  (model$par[1]*model$par[3])/((model$par[1]*model$par[3])+model$par[4]+model$par[5])
  # [1] 0.3239435
  ####
}

