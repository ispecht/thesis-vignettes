# MIT License
#
# Copyright (c) 2023 Ivan Specht
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

## Simulation for Example 3.6
library(parallel)

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

# Degree sequence, starting at root thru 3rd generation (4th generation irrelevant)
ds <- rep(1, 4)

# Coalescent
coalescent <- function(N){
  qs <- get_qs(N)
  choose(N-1, 3) * factorial(3) * (1/N)^sum(ds) * prod((exp((get_qs(N) - 1)*2) * 2^ds) / factorial(ds))
}



# Simulation for N = 10
N <- 10
res_unique <- mclapply(1:1e6, simulate_unique, N = N, mc.cores = 12)
mean(unlist(res_unique))
coalescent(N)

# Simulation for N = 100
N <- 100
res_unique <- mclapply(1:1e6, simulate_unique, N = N, mc.cores = 12)
mean(unlist(res_unique))
coalescent(N)

