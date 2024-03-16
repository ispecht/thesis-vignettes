## Code for Example B.1
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
