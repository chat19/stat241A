library("quantmod")

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

EM_HMM <- function(y, params, max_it = 1000, eps = 10^(-20)){
  A_curr <- params[[1]]
  init_prob_curr <- params[[2]]
  sigmas_curr <- params[[3]]
  num_states <- length(init_prob_curr)
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
    #loglikelihood <- c(loglikelihood, loglik(dt,A_curr,B,f_curr,init_prob, Em,En,Er, gamma))
    if(sum((A_old - A_curr)^2) <= eps^2 & sum((sigmas_old - sigmas_curr)^2) <= eps^2){
      break
    }
  }
  return(list(A_curr, init_prob_curr, sigmas_curr))
  #return(list(A_curr, init_prob_curr, sigmas_curr, loglikelihood))
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

getSymbols("VTSMX",src="yahoo")
training <- compute_returns(VTSMX[seq(1, nrow(VTSMX), by = 30),])
training <- training - mean(training)
results_training <- EM_HMM(training, param_init(2), 150)
results_training






plot(density(training))

training <- rnorm(666, sd = 0.1)
