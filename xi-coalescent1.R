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
