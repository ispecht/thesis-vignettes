## Test of MCMC with overlapping states

# Probability at (0,0), (1,0), (1,1), (0,1)
p <- c(0.1,0.2,0.3,0.4)

# Likelihood
L <- function(state){
  p[state] + p[ifelse(state + 1 == 5, 1, state)]
}

states <- 1 # 1 = bottom, 2 = right, 3 = top, 4 = left

# Number of iterations
N_iters <- 1e5

for (i in 2:N_iters) {
  # If bottom or top, propose left and right states with equal probabilities
  if(states[i-1] %in% c(1,3)){
    prop <- ifelse(runif(1) < 1/2, 2, 4)
  }else{
    prop <- ifelse(runif(1) < 1/2, 1, 3)
  }

  # M-H criterion, symmetric proposal
  if(runif(1) < L(prop)/L(states[i-1])){
    states[i] <- prop
  }else{
    states[i] <- states[i-1]
  }
}

# True probabilities for each state
true_probs <- sapply(1:4, L)
true_probs <- true_probs / sum(true_probs)
print(true_probs)

# Empirical probabilities
print(table(states)/N_iters)

