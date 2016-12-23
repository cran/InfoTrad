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

# Call EHO function
LK_out = EHO(data)
model = optim(par0, LK_out, gr = NULL, method = c("Nelder-Mead"), hessian = FALSE)

# The structure of the model output

model
# $par
# [1] 3.041745e-01 7.608298e-03 4.129369e+02 3.032937e+02 5.261066e+02

# $value
# [1] -22755.78

# $counts
# function gradient 
# 502       NA 

# $convergence
# [1] 1

# $message
# NULL
##

## Parameter Estimates
model$par[1] # Estimate for alpha
# [1] 0.3041745
model$par[2] # Estimate for delta
#[1] 0.007608298
model$par[3] # Estimate for mu
# [1] 412.9369
model$par[4] # Estimate for eb
# [1] 303.2937
model$par[5] # Estimate for es
# [1] 526.1066

## Estimate for PIN
(model$par[1]*model$par[3])/((model$par[1]*model$par[3])+model$par[4]+model$par[5])
# [1] 0.1315227
####

}

