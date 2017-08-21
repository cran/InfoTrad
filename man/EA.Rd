\name{EA}
\alias{EA}
\alias{print.EA_class}
\alias{EA_class}
\alias{print}
\title{Ersan and Alici (2016) Cluster analysis with the altered steps.}
\description{It estimates PIN using Ersan and Alici (2016) modified clustering algorithm.}
\usage{EA(data, likelihood = c("LK", "EHO"))
\method{print}{EA_class}(obj)}
\arguments{
  \item{data}{Data frame with 2 variables}
	\item{likelihood}{Character strings for likelihood algorithm. Default is "LK".}
	\item{obj}{object variable}
}
\details{Argument for data must be a data frame with 2 columns that only contain numbers. Not any other type. You do not have to give names to the columns. We will assign first one as "Buy" and second as "Sell", therefore you should put order numbers with respect to this order.}

\section{Warning}{This function does not handle NA values. Therefore the datasets should not contain any missing value. This function do not conduct the iterative estimation procedure proposed in the same paper.}

\value{
  Returns a list of parameter estimates (output)
  \item{alpha}{A Number}
  \item{delta}{A Number}
  \item{mu}{A Number}
  \item{eb}{A Number}
  \item{es}{A Number}
  \item{LikVal}{A Number}
  \item{PIN}{A Number}
}

\references{
Ersan, Oguz, and Asli Alici . "An unbiased computation methodology for estimating the probability of informed trading (PIN)." Journal of International Financial Markets, Institutions and Money 43 (2016): 74-94.
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
  
  Buy=c(350,250,500,552,163,345,847,923,123,349)
  Sell=c(382,500,463,550,200,323,456,342,578,455)
  data=cbind(Buy,Sell)
  
  # Parameter estimates using the LK factorization of Lin and Ke (2011) 
  # with the modified clustering algorithm of Ersan and Alici (2016).
  # Default factorization is set to be "LK"
  
  result=EA(data)
  print(result)
  
  # Alpha: 0.9511418 
  # Delta: 0.2694005 
  # Mu: 76.7224 
  # Epsilon_b: 493.7045 
  # Epsilon_s: 377.4877 
  # Likelihood Value: 43973.71 
  # PIN: 0.07728924 
  
  
  # Parameter estimates using the EHO factorization of Easley et. al. (2010) 
  # with the modified clustering algorithm of Ersan and Alici (2016).
  
  result=EA(data,likelihood="EHO")
  print(result)
  
  # Alpha: 0.9511418 
  # Delta: 0.2694005 
  # Mu: 76.7224 
  # Epsilon_b: 493.7045 
  # Epsilon_s: 377.4877 
  # Likelihood Value: 43973.71 
  # PIN: 0.07728924 
}

