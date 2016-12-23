\name{GAN}
\alias{GAN}
\title{GAN et al.(2015) Clustering based PIN Estimates}
\description{It estimates PIN using hierarchical agglomertaive clustering.}
\usage{GAN(data, likelihood = c("LK", "EHO"))}
\arguments{
	\item{data}{Data frame with 2 variables}
	\item{likelihood}{Character strings for likelihood algorithm. Default is "LK".}
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
  # with the algorithm of Gan et. al. (2015).
  # Default factorization is set to be "LK"
  
GAN(data)
#             alpha        delta      mu  epsilon_b epsilon_s   LikVal      PIN
#  Output 0.5896898 4.489097e-06 259.886  159.1055  394.6942 19023.96 0.2167481
  
  # Parameter estimates using the EHO factorization of Easley et. al. (2010) 
  # with the algorithm of Gan et. al. (2015)
  
GAN(data, likelihood="EHO")
#             alpha       delta       mu   epsilon_b epsilon_s  LikVal       PIN
#  Output 0.9343141 0.001986236 209.3302  153.7496  429.5933  22688.5 0.2510904

}

