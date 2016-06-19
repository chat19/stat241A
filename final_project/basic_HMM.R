# This script is part of the final project for
# STAT 241A: Statistical Learning Theory.
# The objective is to implement the hidden markov
# model for stock returns. Volatility is treated
# as hidden states, and the observations are the
# stock returns.
#
# Author: Chun Yu Hong (Johnny) (jcyhong@berkeley.edu)

#---------------------------------------------
# Load the necessary packages.
library("quantmod")
library("reshape2")
library("ggplot2")
library(RColorBrewer)

#---------------------------------------------

# Remark: For simplicity, sigmas refers to sigma squared.

AlphaNorm <- function(y, A, initProb, sigmas){
  # Implement the alpha part of the alpha-beta algorithm.
  #
  # Args:
  #   y: the observations
  #   A: the transition matrix
  #   initProb: the initial probability vector
  #   sigmas: the variance parameters
  # Returns:
  #   alpha: the list of alphas
  
  numStates <- length(initProb)
  TT <- length(y)
  alpha <- list()
  alpha1 <- rep(NA, numStates)
  for(i in 1:numStates){
    alpha1[i] <- dnorm(y[[1]], sd = sqrt(sigmas[i])) * initProb[i]
  }
  alpha[[1]] <- alpha1/sum(alpha1)
  for(i in 2:TT){
    alphaTemp <- rep(NA, numStates)
    for(j in 1:numStates){
      a <- rep(NA, numStates)
      for(jj in 1:numStates){
        a[jj] <- alpha[[i-1]][jj]*A[jj, j]*dnorm(y[i], sd = sqrt(sigmas[j]))
      }
      alphaTemp[j] <- sum(a)
    }
    alpha[[i]] <- alphaTemp/sum(alphaTemp)
  }
  alpha
}

# backward <- function(y, A, initProb, sigmas){
#   numStates <- length(initProb)
#   TT <- length(y)
#   beta <- list()
#   beta1 <- rep(1, numStates)
#   beta[[TT]] <- beta1
#   for(i in (TT-1):1){
#     temp <- rep(NA, numStates)
#     for(qT in 1:numStates){
#       temp2 <- 0
#       for(qT1 in 1:numStates){
#         temp2 <- temp2 + beta[[i + 1]][qT1] * 
#           dnorm(y[i + 1], sd = sqrt(sigmas[qT1])) * A[qT,qT1] 
#       }
#       temp[qT] <- temp2
#     }
#     beta[[i]] <- temp
#   }
#   beta
# }


GammaInfer <- function(alpha, A){
  # Implement the gamma part of the alpha-beta algorithm.
  #
  # Args:
  #   alpha: the list of alphas
  #   A: the transition matrix
  # Returns:
  #   gamma: a lust of gammas
  
  numStates <- nrow(A)
  TT <- length(alpha)
  gamma <- list()
  gamma[[TT]] <- alpha[[TT]]
  for(i in (TT - 1):1){
    temp <- rep(NA, numStates)
    for(qT in 1:numStates){
      finalTerms <- rep(NA, numStates)
      for(qT1 in 1:numStates){
        # Compute all the terms in the denominator.
        temp2 <- rep(NA, numStates)
        for(qTemp in 1:numStates){
          temp2[qTemp] <- alpha[[i]][qTemp] * A[qTemp, qT1]
        }
        # Normalize.
        temp2 <- temp2/sum(temp2)
        finalTerms[qT1] <- temp2[qT]* gamma[[i + 1]][qT1]
      }
      temp[qT] <- sum(finalTerms)
    }
    gamma[[i]] <- temp
  }
  gamma
}

XiInfer <- function(alpha, gamma, A, y, sigmas){
  # Implement the xi part of the alpha-beta algorithm.
  #
  # Args:
  #   y: the observations
  #   gamma: a list of gammas
  #   A: the transition matrix
  #   y: a vector of observations
  #   sigmas: the variance parameters
  # Returns:
  #   xi: a list of xis
  
  numStates = nrow(A)
  TT <- length(alpha)
  # xi is a list of four-dimensional arrays.
  xi <- list()
  for(i in 1:(TT-1)){
    temp <- matrix(0, nrow = numStates, ncol = numStates)
    for(qT in 1:numStates){
      for(qT1 in 1:numStates){
        temp[qT, qT1] <- alpha[[i]][qT]*dnorm(y[i+1], sd = sqrt(sigmas[qT1]))*
          gamma[[i + 1]][qT1]* A[qT, qT1] / alpha[[i+1]][qT1]
      }
    }
    xi[[i]] <- temp/sum(temp)
  }
  xi
}

loglikHMM <- function(y, A, sigmas, initProb, gamma, xi){
  # Compute the log-likelihood of the HMM model.
  #
  # Args:
  #   y: the observations
  #   A: the transition matrix
  #   sigmas: the variance parameters
  #   initProb: the initial probability vector
  #   gamma: a list of gammas
  #   xi: a list of xis
  # Returns:
  #   loglik: the log-likelihood of the HMM model
  
  firstTerm <- sum(initProb * log(initProb))
  
  secondTerm <- sum(Reduce("+", xi) * log(A))
  
  temp <- lapply(1:length(gamma), function(i){
    gamma[[i]] * (y[i])^2
  })
  Er <- Reduce("+", temp)
  thirdTerm <- sum(-1/(2 * sigmas) * Er)
  
  En <- Reduce("+", gamma)
  fourthTerm <- sum((-1/2 * log(2 * pi) - 1/2 * log(sigmas)) * En)
  
  loglik <- firstTerm + secondTerm + thirdTerm + fourthTerm
  loglik
}



EM_HMM <- function(y, params, maxIt = 1000, eps = 10^(-20)){
  # Run the EM algorithm for the hidden Markov model.
  #
  # Args:
  #   y: a vector of observations
  #   params: a list of parameters of the HMM
  #   maxIt: the maximum iterations of EM
  #   eps: the tolerance
  # Returns
  #   a list containing the parameter estimates and the log-likelihood
  
  ACurr<- params[[1]]
  initProbCurr <- params[[2]]
  sigmasCurr <- params[[3]]
  numStates <- length(initProbCurr)
  loglikelihood <- c()
  for(tt in 1:maxIt){
    AOld <- ACurr
    sigmasOld <- sigmasCurr
    initProbOld <- initProbCurr
    # HMM inference.
    alpha <- AlphaNorm(y, ACurr, initProbCurr , sigmasCurr )
    gamma <- GammaInfer(alpha, ACurr)
    xi <- XiInfer(alpha, gamma, ACurr, y, sigmasCurr )
    
    # E step
    En <- Reduce("+", gamma)
    
    Em <- Reduce("+", xi)
    
    temp2 <- lapply(1:length(gamma), function(i){
      gamma[[i]] * (y[i])^2
    })
    Er <- Reduce("+", temp2)
    
    # M step 
    ACurr <- Em/rowSums(Em)
    sigmasCurr  <- Er/En
    initProbCurr  <- gamma[[1]]
    loglikelihood <- c(loglikelihood, 
                       loglik_HMM(y, ACurr, sigmasCurr, initProbCurr, gamma, xi))
    if(sum((AOld - ACurr)^2) <= eps^2 & sum((sigmasOld - sigmasCurr)^2) <= eps^2){
      break
    }
  }
  #return(list(ACurr, initProbCurr, sigmasCurr))
  return(list(ACurr, initProbCurr, sigmasCurr, loglikelihood))
}

ComputeReturns <- function(prices){
  # Compute the stock returns based on prices.
  #
  # Args:
  #   prices: a vector of prices
  # Returns:
  #   a vector of returns, with length one less than prices
  prices <- as.numeric(prices)
  diff(prices)/prices[-length(prices)]
}

ParamInit <- function(numStates){
  # Returns an initial guess of the parameters.
  #
  # Args:
  #   numStates: the number of hidden states
  # Returns:
  #   a list of parameters
  A <- matrix(1/numStates, nrow = numStates, ncol = numStates)
  initProb <- rep(1/numStates, numStates)
  sigmas <- seq(0.1, 0.9, 0.8/(numStates - 1))
  list(A, initProb, sigmas )
}

ComputeReturns <- function(df){
  # Compute returns based on a data frame, in which the
  # sixth column contains prices.
  #
  # Args:
  #   df: a data frame in which the
  #       sixth column contains prices
  # Returns:
  #   a vector of returns
  prices <- df[,6]
  prices <- as.numeric(prices)
  diff(prices)/prices[-length(prices)]
}

# Viterbi algorithm.
FindMaxConfig <- function(y, params){
  # Find the configuration with the highest likelihood
  #
  # Args:
  #   y: a vector of observations
  #   params: a vector of parameters
  # Returns:
  #   a list containing the messages and the max-likelihood
  #   configuration
  TT <- length(y)
  A <- params[[1]]
  initProb <- params[[2]]
  sigmas <- params[[3]]
  numStates <- length(initProb)
  
  messages <- list()
  config <- list()
  
  messageQT1 <- rep(NA, numStates)
  configQT1 <- rep(NA, numStates)
  for(qT1 in 1:numStates){
    temp <- rep(NA, numStates)
    for(qT in 1:numStates){
      temp[qT] <- A[qT1, qT] * dnorm(y[TT], sd = sqrt(sigmas[qT]))
    }
    messageQT1[qT1] <- max(temp)
    configQT1[qT1] <- which.max(temp)
  }
  messages[[TT]] <- messageQT1
  config[[TT]] <- configQT1
  
  for(i in (TT - 1):2){
    messageqI1 <- rep(NA, numStates)
    config_qI1 <- rep(NA, numStates)
    for(qI1 in 1:numStates){
      temp <- rep(NA, numStates)
      for(q_i in 1:numStates){
        temp[q_i] <- A[qI1, q_i] * dnorm(y[i], sd = sqrt(sigmas[q_i])) *
          messages[[i + 1]][q_i]
      }
      messageqI1[qI1] <- max(temp)
      config_qI1[qI1] <- which.max(temp)
    }
    messages[[i]] <- messageqI1
    config[[i]] <- config_qI1
  }
  
  lik <- rep(NA, numStates)
  q1State <- rep(NA, numStates)
  allLik <- sapply(1:numStates, function(q1){
    initProb[q1] * dnorm(y[1], sd = sqrt(sigmas[q1])) * messages[[2]][q1]
  })
  messages[[1]] <- max(allLik)
  config[[1]] <- which.max(allLik)
  
  return(list(messages, config))
}

GetSeq <- function(config){
  # Extract the sequence of max-likelihood states.
  #
  # Args:
  #   config: the configuation of states
  # Returns:
  #   sequ: the sequence of max-likelihood states
  sequ <- config[[2]][[1]]
  # Trace the path.
  for(i in 2:length(config[[2]])){
    sequ[i] <- config[[2]][[i]][sequ[i - 1]]
  }
  sequ
}

ExtractInfo <- function(ticker, interval = 10, centered = T, numRm = 0){
  # Extract the stock data from yahoo based on ticker.
  #
  # Args:
  #   ticker: a character vector 
  #   interval: the time interval between consecutive prices
  #   centered: a boolean indicating whether the prices are centered
  #   numRm: number of observations to be removed
  #
  # Returns:
  #   info: a list of (centered) prices and the data frame
  temp <- getSymbols(ticker,src="yahoo", env = NULL)
  temp <- temp[seq(1, nrow(temp), by = interval),]
  if(numRm > 0){
    n = length(temp)
    temp <- temp[-((n - numRm + 1):n), ]
  }
  info <- list()
  info[[1]] <- ComputeReturns(temp)
  if(centered){
    #avg <- prod(1 + info[[1]]) ^ (1/length(info[[1]])) - 1
    avg <- mean(info[[1]])
    info[[1]] <- info[[1]] - avg
  }
  info[[2]] <- temp
  return(info)
}

InfoAnalysis <- function(info, numStates, ticker){
  # Create plots of the stock, including results
  # from the HMM.
  #
  # Args:
  #   info: a list of info
  #   numStates: an integer indicating the number of hidden states
  #   ticker: a character vector
  # Returns:
  #   a histogram of centered returns
  #   an acf plot of the returns
  #   a plot of the predicted prices
  #   a plot of the volatility stages
  training <- info[[1]]
  resultsTraining <- EM_HMM(training, ParamInit(numStates), 150)
  config <- FindMaxConfig(training, resultsTraining)
  info[[2]] <- cbind(info[[2]], c(NA, GetSeq(config)))
  plot(GetSeq(config), type = "l")
  cols = rainbow(numStates)
  
  par(mfrow = c(2,2), oma = c(0, 0, 2, 0))
  hist(info[[1]], main = "Centered Returns", xlab = "Centered Returns")
  acf(info[[1]], main = "ACF plot of Centered Returns")
  plot(info[[2]][,6], main = "Adjusted Prices", ylab = "Adjusted Prices")
  plot(info[[2]][,7], main = "Volatility", ylab = "Volatility Stage")
  mtext(ticker, outer = TRUE, cex = 1.5)
}


LikelihoodAnalysis <- function(returns, numStatesVec){
  # Plot the log-likelihood against the number of hidden states.
  # 
  # Args:
  #   returns: a vector of returns
  #   numStatesVec: a vector of numbers of hidden states
  # Returns:
  #   logliks: a vector of log-likelihoods
  
  logliks <- sapply(numStatesVec, function(numStates){
    results <- EM_HMM(returns, ParamInit(numStates), 150)
    max(results[[4]])
  })
  plot(numStatesVec, logliks, main = "Log likelihood for different\nnumber of states")
  return(logliks)
}

FindBICs <- function(logliks, returns, numStatesVec){
  # Compute the BICs based on the number of states and log-likelihoods.
  #
  # Args:
  #   logliks: a vector of log-likelihoods
  #   returns: a vector of returns
  #   numStatesVec: a vector of numbers of hidden states
  # Returns:
  #   a plot of the BICs against the number of states in the HMM
  #   a vector of BICs
  num_obs <- length(returns)
  numOfParams <- sapply(numStatesVec, function(i){
    numInitProb <- i - 1
    numA <- i * (i - 1)
    numSigmas <- i
    numInitProb + numA + numSigmas
  }) 
  BICs <- -2 * logliks + numOfParams * log(num_obs)
  plot(numStatesVec, BICs, main = "BICs for different number of states",
       type = "b", xlab = "Number of States")
  BICs
}


Simulation <- function(numRm, mu, centeredReturns){
  # Run simulations based on HMM.
  # 
  # Args:
  #   numRm: number of observations removed from the observations
  #   mu: the mean parameter
  #   centeredReturns: a vector of centered returns
  # Returns:
  #   a vector of simulated returns
  training <- centeredReturns
  resultsTraining <- EM_HMM(training, ParamInit(2), 150)
  A <- resultsTraining[[1]]
  sigmas <- resultsTraining[[3]]
  config <- FindMaxConfig(training, resultsTraining)
  maxLikSeq <- GetSeq(config)
  curr <- maxLikSeq[length(maxLikSeq)]
  predictedSeq <- rep(NA, numRm)
  for(i in 1:numRm){
    if(runif(1) <= A[curr, 1]){
      predictedSeq[i] <- 1
    }else{
      predictedSeq[i] <- 2
    }
    curr <- predictedSeq[i] 
  }
  mu + rnorm(numRm, sd = sqrt(sigmas[predictedSeq]))
}



info <- ExtractInfo("AMZN")
InfoAnalysis(info, 2, "AMZN")
dev.off()
resultsTraining <- EM_HMM(info[[1]], ParamInit(2), 150)
logliks <- LikelihoodAnalysis(info[[1]], 2:8)
BICs <- FindBICs(logliks, info[[1]], 2:8)
plot(2:8, BICs, main = "BICs for different number of states (AMZN)",
     type = "b", xlab = "Number of States")

info2 <- ExtractInfo("PLNR")
InfoAnalysis(info2, 3, "PLNR")
dev.off()
resultsTraining2 <- EM_HMM(info2[[1]], ParamInit(3), 150)
logliks2 <- LikelihoodAnalysis(info2[[1]], 2:8)
BICs2 <- FindBICs(logliks2, info2[[1]], 2:8)

png("BICs")
plot(2:8, BICs, main = "BICs for different number of states",
     type = "b", xlab = "Number of States", ylab = "BIC",
     ylim = c(min(c(BICs, BICs2)), max(c(BICs, BICs2))))
lines(2:8, BICs2, type = "b", lty = 2, col = "red")
points(2, BICs[1], pch = 16)
points(3, BICs2[2], pch = 16, col = "red")
legend("topleft", legend = c("AMZN", "PLNR"), lty = c(1,2), col = c("black", "red"))
dev.off()


GetVolatility <- function(ticker){
  # Extract the estimated volatility based on the HMM.
  #
  # Args:
  #   ticker: a character vector
  # Returns:
  #   a vector of volatility stages
  info <- ExtractInfo(ticker)
  logliks <- LikelihoodAnalysis(info[[1]], 2:6)
  BICs <- FindBICs(logliks, info[[1]], 2:6)
  numVS <- which.min(BICs) + 1
  results <- EM_HMM(info[[1]], ParamInit(numVS), 150)
  results[[3]]
}

# Mega-cap stocks.
tickerList <- c("AMZN", "AAPL", "MSFT", "XOM", "T",
                 "PLNR", "HSC", "TREX", "ANF", "ALK")
volatilityAll <- lapply(tickerList, GetVolatility)
sqrtVar <- lapply(volatilityAll, sqrt)

# Small-cap stocks.
tickerList2 <- c("JNJ", "WFC", "GE", "CVX", "WMT",
                 "UBCP", "CIZN", "KFFB", "DOM", "CSPI")
volatilityAll2 <- lapply(tickerList2, GetVolatility)
sqrtVar2 <- lapply(volatilityAll2, sqrt)

sqrtVar20 <- c(sqrtVar, sqrtVar2)

png("sigsvol.png")
plot(NULL, xlim=c(1,3), ylim=c(0,0.55), ylab="sigma", xlab="Volatility Stage (BIC-optimal)",
     xaxt = "n", main = "Sigmas vs Volatility stages",
     cex.main=1.5, cex.lab=1.5,cex.axis=1.5)
cols <- c(rep(c("red", "blue"), each = 5), rep(c("red", "blue"), each = 5))
sapply(1:20, function(i){
  vols <- sqrtVar20[[i]]
  lines(1:length(vols), vols, col = cols[i], type = "b")
})
axis(1, at = 1:3)
legend("topleft", c("Mega cap", "Small cap"), col = c("red", "blue"), lty = c(1,1))
dev.off()

# Study Amazon (AMZN) more carefully.
trainingInfo <- ExtractInfo("AMZN", numRm = 10)
mu <- mean(trainingInfo[[1]])
centeredReturns <- trainingInfo[[1]]
forecast1 <- Simulation(10, mu, centeredReturns)
forecast2 <- Simulation(10, mu, centeredReturns)
forecast3 <- Simulation(10, mu, centeredReturns)

actual <- info[[1]][215:224]

dframe <- data.frame(forecast1, forecast2, forecast3,
                     actual, date = index(info[[2]])[216:225])
dframeLong <- melt(dframe, id="date")  # convert to long format

ggplot(data=dframeLong,
       aes(x=date, y=value, colour=variable)) +
  geom_line() + ylab("Returns") + xlab("Date") + 
  ggtitle("AMZN Predicted Returns vs Actual Returns")