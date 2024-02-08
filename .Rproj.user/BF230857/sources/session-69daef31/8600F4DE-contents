## Simulation for Example 3.2
library(parallel)

# Helper function that draws the x offspring of a node i
offspring <- function(i, x){
  return(sample(setdiff(0:10, i), x))
}

# Do the offspring of a node contain nodes 2 and 3?
success <- function(i,x){
  if(all(
    c(2,3) %in% offspring(i,x)
  )){
    return(T)
  }else{
    return(F)
  }
}

# Function that returns TRUE if the finite-population Galton-Watson process contains an isomorphic subtree to (R,g) such that kappa(L) = Lambda(phi(L)); FALSE otherwise.
simulate <- function(z){
  # Draw the number of nodes in generation 1
  x <- rpois(1, 1)
  # If 0, return FALSE
  if(x == 0){
    return(FALSE)
  }else{
    # Nodes in s1
    s1 <- offspring(0, x)
    # For each node in s1, draw the number of offspring
    xs <- rpois(x, 1)
    # Children of a node in s1 contains 2 and 3?
    if(any(
      mapply(success, s1, xs)
    )){
      return(T)
    }else{
      return(F)
    }
  }
}

res <- mclapply(1:1e6, simulate, mc.cores = 12)
mean(unlist(res))
