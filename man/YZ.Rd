\name{YZ}
\alias{YZ}
\title{
Yan and Zhang (2012) Grid-Search based PIN Estimates
}
\description{It estimates PIN using Yan and Zhang (2012) algorithm.}
\usage{
YZ(data, likelihood = c("LK", "EHO"))
}
\arguments{
  \item{data}{Data frame with 2 variables}
	\item{likelihood}{Character strings for likelihood algorithm. Default is "LK".}
}
\details{Argument for data must be a data frame with 2 columns that only contain numbers. Not any other type. You do not have to give names to the columns. We will assign first one as "Buy" and second as "Sell", therefore you should put order numbers with respect to this order.}

\section{Warning}{This function does not handle NA values. Therefore the datasets should not contain any missing value}

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
Y. Yan and S. Zhang. An improved estimation method and empirical properties of the probability of informed trading. Journal of Banking & Finance, 36(2):454-467, 2012.
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
  
  # Parameter estimates using the LK factorization of Lin and Ke (2011) 
  # with the algorithm of Yan and Zhang (2012).
  # Default factorization is set to be "LK"
  
  YZ(data)
  #          alpha        delta      mu   epsilon_b epsilon_s  LikVal      PIN
  # Output 0.5876696 0.0005934669 309.688  156.3483  420.2923  19035.72 0.239897
  
  # Parameter estimates using the EHO factorization of Easley et. al. (2010) 
  # with the algorithm of Yan and Zhang (2012).
  
  YZ(data,likelihood="EHO")
  #            alpha     delta     mu    epsilon_b epsilon_s   LikVal      PIN
  # Output 0.4156183 0.5079664 355.9734  268.2996   261.455  23091.94 0.2183094
}
