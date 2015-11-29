library("quantmod")
library(RColorBrewer)

# Remark: sigmas refers to sigma squared.

alpha_norm <- function(y, A, init_prob, sigmas){
  num_states <- length(init_prob)
  TT <- length(y)
  alpha <- list()
  alpha1 <- rep(NA, num_states)
  for(i in 1:num_states){
    alpha1[i] <- dnorm(y[[1]], sd = sqrt(sigmas[i])) * init_prob[i]
  }
  alpha[[1]] <- alpha1/sum(alpha1)
  for(i in 2:TT){
    alpha_temp <- rep(NA, num_states)
    for(j in 1:num_states){
      a <- rep(NA, num_states)
      for(jj in 1:num_states){
        a[jj] <- alpha[[i-1]][jj]*A[jj, j]*dnorm(y[i], sd = sqrt(sigmas[j]))
      }
      alpha_temp[j] <- sum(a)
    }
    alpha[[i]] <- alpha_temp/sum(alpha_temp)
  }
  alpha
}

# backward <- function(y, A, init_prob, sigmas){
#   num_states <- length(init_prob)
#   TT <- length(y)
#   beta <- list()
#   beta1 <- rep(1, num_states)
#   beta[[TT]] <- beta1
#   for(i in (TT-1):1){
#     temp <- rep(NA, num_states)
#     for(q_t in 1:num_states){
#       temp2 <- 0
#       for(q_t1 in 1:num_states){
#         temp2 <- temp2 + beta[[i + 1]][q_t1] * 
#           dnorm(y[i + 1], sd = sqrt(sigmas[q_t1])) * A[q_t,q_t1] 
#       }
#       temp[q_t] <- temp2
#     }
#     beta[[i]] <- temp
#   }
#   beta
# }


gamma_infer <- function(alpha, A){
  num_states <- nrow(A)
  TT <- length(alpha)
  gamma <- list()
  gamma[[TT]] <- alpha[[TT]]
  for(i in (TT - 1):1){
    temp <- rep(NA, num_states)
    for(q_t in 1:num_states){
      final_terms <- rep(NA, num_states)
      for(q_t1 in 1:num_states){
        # Compute all the terms in the denominator.
        temp2 <- rep(NA, num_states)
        for(q_temp in 1:num_states){
          temp2[q_temp] <- alpha[[i]][q_temp] * A[q_temp, q_t1]
        }
        # Normalize.
        temp2 <- temp2/sum(temp2)
        final_terms[q_t1] <- temp2[q_t]* gamma[[i + 1]][q_t1]
      }
      temp[q_t] <- sum(final_terms)
    }
    gamma[[i]] <- temp
  }
  gamma
}

xi_infer <- function(alpha, gamma, A, y, sigmas){
  num_states = nrow(A)
  TT <- length(alpha)
  # xi is a list of four-dimensional arrays.
  xi <- list()
  for(i in 1:(TT-1)){
    temp <- matrix(0, nrow = num_states, ncol = num_states)
    for(q_t in 1:num_states){
      for(q_t1 in 1:num_states){
        temp[q_t, q_t1] <- alpha[[i]][q_t]*dnorm(y[i+1], sd = sqrt(sigmas[q_t1]))*
          gamma[[i + 1]][q_t1]* A[q_t, q_t1] / alpha[[i+1]][q_t1]
      }
    }
    xi[[i]] <- temp/sum(temp)
  }
  xi
}

loglik_HMM <- function(y, A, sigmas, init_prob, gamma, xi){
  first_term <- sum(init_prob * log(init_prob))
  second_term <- sum(Reduce("+", xi) * log(A))
  temp <- lapply(1:length(gamma), function(i){
    gamma[[i]] * (y[i])^2
  })
  Er <- Reduce("+", temp)
  third_term <- sum(-1/(2 * sigmas) * Er)
  En <- Reduce("+", gamma)
  fourth_term <- sum((-1/2 * log(2 * pi) - 1/2 * log(sigmas)) * En)
  loglik <- first_term + second_term + third_term + fourth_term
  loglik
}



EM_HMM <- function(y, params, max_it = 1000, eps = 10^(-20)){
  A_curr <- params[[1]]
  init_prob_curr <- params[[2]]
  sigmas_curr <- params[[3]]
  num_states <- length(init_prob_curr)
  loglikelihood <- c()
  for(tt in 1:max_it){
    A_old <- A_curr
    sigmas_old <- sigmas_curr
    init_prob_old <- init_prob_curr
    # HMM inference.
    alpha <- alpha_norm(y, A_curr , init_prob_curr , sigmas_curr )
    gamma <- gamma_infer(alpha, A_curr )
    xi <- xi_infer(alpha, gamma, A_curr, y, sigmas_curr )
    
    # E step
    En <- Reduce("+", gamma)
    
    Em <- Reduce("+", xi)
    
    temp2 <- lapply(1:length(gamma), function(i){
      gamma[[i]] * (y[i])^2
    })
    Er <- Reduce("+", temp2)
    
    # M step 
    A_curr  <- Em/rowSums(Em)
    sigmas_curr  <- Er/En
    init_prob_curr  <- gamma[[1]]
    loglikelihood <- c(loglikelihood, 
                       loglik_HMM(y, A_curr, sigmas_curr, init_prob_curr, gamma, xi))
    if(sum((A_old - A_curr)^2) <= eps^2 & sum((sigmas_old - sigmas_curr)^2) <= eps^2){
      break
    }
  }
  #return(list(A_curr, init_prob_curr, sigmas_curr))
  return(list(A_curr, init_prob_curr, sigmas_curr, loglikelihood))
}

compute_returns <- function(prices){
  prices <- as.numeric(prices)
  diff(prices)/prices[-length(prices)]
}

param_init <- function(num_states){
  A <- matrix(1/num_states, nrow = num_states, ncol = num_states)
  init_prob <- rep(1/num_states, num_states)
  sigmas <- seq(0.1, 0.9, 0.8/(num_states - 1))
  list(A, init_prob, sigmas )
}

compute_returns <- function(df){
  prices <- df[,6]
  prices <- as.numeric(prices)
  diff(prices)/prices[-length(prices)]
}

# Viterbi algorithm.
find_max_config <- function(y, params){
  TT <- length(y)
  A <- params[[1]]
  init_prob <- params[[2]]
  sigmas <- params[[3]]
  num_states <- length(init_prob)
  
  messages <- list()
  config <- list()
  
  message_q_T_1 <- rep(NA, num_states)
  config_q_T_1 <- rep(NA, num_states)
  for(q_T_1 in 1:num_states){
    temp <- rep(NA, num_states)
    for(q_T in 1:num_states){
      temp[q_T] <- A[q_T_1, q_T] * dnorm(y[TT], sd = sqrt(sigmas[q_T]))
    }
    message_q_T_1[q_T_1] <- max(temp)
    config_q_T_1[q_T_1] <- which.max(temp)
  }
  messages[[TT]] <- message_q_T_1
  config[[TT]] <- config_q_T_1
  
  for(i in (TT - 1):2){
    message_q_i_1 <- rep(NA, num_states)
    config_q_i_1 <- rep(NA, num_states)
    for(q_i_1 in 1:num_states){
      temp <- rep(NA, num_states)
      for(q_i in 1:num_states){
        temp[q_i] <- A[q_i_1, q_i] * dnorm(y[i], sd = sqrt(sigmas[q_i])) *
          messages[[i + 1]][q_i]
      }
      message_q_i_1[q_i_1] <- max(temp)
      config_q_i_1[q_i_1] <- which.max(temp)
    }
    messages[[i]] <- message_q_i_1
    config[[i]] <- config_q_i_1
  }
  
  lik <- rep(NA, num_states)
  q_1_state <- rep(NA, num_states)
  all_lik <- sapply(1:num_states, function(q_1){
    init_prob[q_1] * dnorm(y[1], sd = sqrt(sigmas[q_1])) * messages[[2]][q_1]
  })
  messages[[1]] <- max(all_lik)
  config[[1]] <- which.max(all_lik)
  
  return(list(messages, config))
}

# Extract the sequence of max-likelihood states.
get_seq <- function(config){
  sequ <- config[[2]][[1]]
  # Trace the path.
  for(i in 2:length(config[[2]])){
    sequ[i] <- config[[2]][[i]][sequ[i - 1]]
  }
  sequ
}

extract_info <- function(ticker, interval = 10, centered = T){
  temp <- getSymbols(ticker,src="yahoo", env = NULL)
  temp <- temp[seq(1, nrow(temp), by = interval),]
  info <- list()
  info[[1]] <- compute_returns(temp)
  if(centered){
    #avg <- prod(1 + info[[1]]) ^ (1/length(info[[1]])) - 1
    avg <- mean(info[[1]])
    info[[1]] <- info[[1]] - avg
  }
  info[[2]] <- temp
  return(info)
}

info_analysis <- function(info, num_states, ticker){
  training <- info[[1]]
  results_training <- EM_HMM(training, param_init(num_states), 150)
  config <- find_max_config(training, results_training)
  info[[2]] <- cbind(info[[2]], c(NA, get_seq(config)))
  plot(get_seq(config), type = "l")
  cols = rainbow(num_states)
  
  par(mfrow = c(2,2), oma = c(0, 0, 2, 0))
  hist(info[[1]], main = "Centered Returns", xlab = "Centered Returns")
  acf(info[[1]], main = "ACF plot of Centered Returns")
  plot(info[[2]][,6], main = "Adjusted Prices", ylab = "Adjusted Prices")
  plot(info[[2]][,7], main = "Volatility", ylab = "Volatility Stage")
  mtext(ticker, outer = TRUE, cex = 1.5)
}


likelihood_analysis <- function(returns, num_states_vec){
  logliks <- sapply(num_states_vec, function(num_states){
    results <- EM_HMM(returns, param_init(num_states), 150)
    max(results[[4]])
  })
  plot(num_states_vec, logliks, main = "Log likelihood for different\nnumber of states")
  return(logliks)
}

find_BICs <- function(logliks, returns, num_states_vec){
  num_obs <- length(returns)
  num_of_params <- sapply(num_states_vec, function(i){
    num_init_prob <- i - 1
    num_A <- i * (i - 1)
    num_sigmas <- i
    num_init_prob + num_A + num_sigmas
  }) 
  BICs <- -2 * logliks + num_of_params * log(num_obs)
  plot(num_states_vec, BICs, main = "BICs for different number of states",
       type = "b", xlab = "Number of States")
  BICs
}


info <- extract_info("AMZN")
info_analysis(info, 2, "AMZN")
dev.off()
results_training <- EM_HMM(info[[1]], param_init(2), 150)
logliks <- likelihood_analysis(info[[1]], 2:8)
BICs <- find_BICs(logliks, info[[1]], 2:8)
plot(2:8, BICs, main = "BICs for different number of states (AMZN)",
     type = "b", xlab = "Number of States")

info2 <- extract_info("PLNR")
info_analysis(info2, 3, "PLNR")
dev.off()
results_training2 <- EM_HMM(info2[[1]], param_init(3), 150)
logliks2 <- likelihood_analysis(info2[[1]], 2:8)
BICs2 <- find_BICs(logliks2, info2[[1]], 2:8)

png("BICs")
plot(2:8, BICs, main = "BICs for different number of states",
     type = "b", xlab = "Number of States", ylab = "BIC",
     ylim = c(min(c(BICs, BICs2)), max(c(BICs, BICs2))))
lines(2:8, BICs2, type = "b", lty = 2, col = "red")
points(2, BICs[1], pch = 16)
points(3, BICs2[2], pch = 16, col = "red")
legend("topleft", legend = c("AMZN", "PLNR"), lty = c(1,2), col = c("black", "red"))
dev.off()


get_volatility <- function(ticker){
  info <- extract_info(ticker)
  logliks <- likelihood_analysis(info[[1]], 2:6)
  BICs <- find_BICs(logliks, info[[1]], 2:6)
  num_VS <- which.min(BICs) + 1
  results <- EM_HMM(info[[1]], param_init(num_VS), 150)
  results[[3]]
}

ticker_list <- c("AMZN", "AAPL", "MSFT", "XOM", "T",
                 "PLNR", "HSC", "TREX", "ANF", "ALK")
volatility_all <- lapply(ticker_list, get_volatility)
sqrt_var <- lapply(volatility_all, sqrt)

ticker_list2 <- c("JNJ", "WFC", "GE", "CVX", "WMT",
                 "UBCP", "CIZN", "KFFB", "DOM", "CSPI")
volatility_all2 <- lapply(ticker_list2, get_volatility)
sqrt_var2 <- lapply(volatility_all2, sqrt)

sqrt_var_20 <- c(sqrt_var, sqrt_var2)

plot(NULL, xlim=c(1,3), ylim=c(0,0.55), ylab="sigma", xlab="Volatility Stage (BIC-optimal)",
     xaxt = "n", main = "Sigmas vs Volatility stages")
cols <- c(rep(c("red", "blue"), each = 5), rep(c("red", "blue"), each = 5))
sapply(1:20, function(i){
  vols <- sqrt_var_20[[i]]
  lines(1:length(vols), vols, col = cols[i], type = "b")
})
axis(1, at = 1:3)
legend("topleft", c("Mega cap", "Small cap"), col = c("red", "blue"), lty = c(1,1))

