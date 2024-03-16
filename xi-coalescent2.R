## Simulation for Example 3.6
library(parallel)

# Function that returns TRUE if the labeled Galton-Watson process contains an isomorphic subtree to (R,g) such that A(L) = B(phi(L)); FALSE otherwise.
simulate <- function(z, N){
  # Draw the number of nodes in generation i, n[i]
  n <- c()
  for (i in 1:4) {
    if(i == 1){
      n[i] <- rpois(1, 2)
    }else{
      # Sum of k i.i.d. Pois(2) r.v.s is Pois(2k)
      n[i] <- rpois(1, 2*n[i-1])
    }
  }
  # Whether generation 4 includes the label \kappa(r_{41}) = 4 is Bernoulli with probability 1 - (9/10)^n[4]
  if(runif(1) < 1 - ((N-1)/N)^n[4]){
    return(T)
  }else{
    return(F)
  }
}

res <- mclapply(1:1e6, simulate, N = 10, mc.cores = 12)
mean(unlist(res))

## Repeat experiment for N = 100
res <- mclapply(1:1e6, simulate, N = 100, mc.cores = 12)
mean(unlist(res))

# Function that returns TRUE if the labeled Galton-Watson process contains a UNIQUE isomorphic subtree to (R,g) such that A(L) = B(phi(L)); FALSE otherwise.
simulate_unique <- function(z, N){
  # Draw the number of nodes in generation i, n[i]
  n <- c()
  # Ancestor function
  h <- list()
  # Labels
  labs <- list()
  for (i in 1:4) {
    if(i == 1){
      n[i] <- rpois(1, 2)
      h[[i]] <- rep(0, n[i])
    }else{
      h[[i]] <- integer(0)
      if(n[i-1] > 0){
        for (j in 1:n[i-1]) {
          h[[i]] <- c(h[[i]], rep(j, rpois(1, 2)))
        }
      }
      n[i] <- length(h[[i]])
    }
    labs[[i]] <- sample(1:N, n[i], replace = T)
  }
  #labs

  # Get the ancestry of person 1
  if(sum(1 == labs[[4]]) != 1){
    return(FALSE)
  }else{
    anc <- which(labs[[4]] == 1)
    for (i in 4:2) {
      anc <- c(h[[i]][anc[1]], anc)
    }
  }
  all_labs <- unlist(labs)
  is_unique <- c()
  for (i in 1:4) {
    if(sum(labs[[i]][anc[i]] == all_labs) == 1){
      is_unique[i] <- T
    }else{
      is_unique[i] <- F
    }
  }
  if(all(is_unique)){
    return(T)
  }else{
    return(F)
  }
}

N <- 10

res_unique <- mclapply(1:1e6, simulate_unique, N = N, mc.cores = 12)
mean(unlist(res_unique))

# Poisson(2) MGF
pois_mgf <- function(t){
  exp(2 * (exp(t) - 1))
}

# q1 through q4
get_qs <- function(N){
  q <- rep(NA, 4)
  q[4] <- 1 - 4/N
  for (i in 3:1) {
    q[i] <- (1 - 4/N) * pois_mgf(log(q[i+1]))
  }
  q
}
qs <- get_qs(N)

# Degree sequence, starting at root thru 3rd generation (4th generation irrelevant)
ds <- rep(1, 4)

# Coalescent
coalescent <- function(N){
  choose(N-1, 3) * factorial(3) * (1/N)^sum(ds) * prod((exp((get_qs(N) - 1)*2) * 2^ds) / factorial(ds))
}
coalescent(N)

