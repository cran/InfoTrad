\name{GAN}
\alias{GAN}
\alias{print.GAN_class}
\alias{GAN_class}
\title{GAN et al.(2015) Clustering based PIN Estimates}
\description{It estimates PIN using hierarchical agglomertaive clustering.}
\usage{GAN(data, likelihood = c("LK", "EHO"))
\method{print}{GAN_class}(obj)}
\arguments{
	\item{data}{Data frame with 2 variables}
	\item{likelihood}{Character strings for likelihood algorithm. Default is "LK".}
	\item{obj}{object variable}
}
\details{Argument for data must be a data frame with 2 columns that only contain numbers. Not any other type. You do not have to give names to the columns. We will assign first one as "Buy" and second as "Sell", therefore you should put order numbers with respect to this order. This package uses the hclust() function of Mullner (2013) to cluster the data at default settings.}

\section{Warning}{This function does not handle NA values. Therefore, the dataset should not contain any missing values.}

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

D. Mullner. fastcluster: Fast hierarchical, agglomerative clustering routines for r and python. Journal of Statistical Software, 53(9):1-18, 2013.

Gan, Q., Wei, W. C., & Johnstone, D. A faster estimation method for the probability of informed trading using hierarchical agglomerative clustering. Quantitative Finance, 15(11), 1805-1821, 2015.}
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
  data<-cbind(Buy,Sell)

  # Parameter estimates using the LK factorization of Lin and Ke (2011) 
  # with the algorithm of Gan et. al. (2015).
  # Default factorization is set to be "LK"
  
  result=GAN(data)
  print(result)
  
  # Alpha: 0.3999998 
  # Delta: 0 
  # Mu: 442.1667 
  # Epsilon_b: 263.3333 
  # Epsilon_s: 424.9 
  # Likelihood Value: 44371.84 
  # PIN: 0.2044464 
  
  # Parameter estimates using the EHO factorization of Easley et. al. (2010) 
  # with the algorithm of Gan et. al. (2015)
  
  result=GAN(data, likelihood="EHO")
  print(result)
  
  # Alpha: 0.3230001 
  # Delta: 0.4780001 
  # Mu: 481.3526 
  # Epsilon_b: 356.6359 
  # Epsilon_s: 313.136 
  # Likelihood Value: Inf 
  # PIN: 0.1884001 


}

