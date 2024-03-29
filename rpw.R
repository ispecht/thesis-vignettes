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

### Compare empirical and approximate theoretical distrubutions for the RPW model
library(parallel)
library(extraDistr)
library(ggplot2)
library(tikzDevice)
library(cowplot)

## Simulate the RPW model, starting with 0 white balls and 1 black ball
# Assume p = p_W = 3p_B
# An efficient way of simulating this: select the steps at which the mutations happen as Geom(3p) draws
# At these steps, choose a ball. If white, return black. If black, return white with probability 1/3
rpw <- function(x, p, n){
  # Number of white balls
  w <- 0
  # Number of total balls
  tot <- 1
  while (tot < n) {
    # How many balls added to the urn before the transition?
    before <- rgeom(1, p)
    if(tot + before >= n){
      before <- n - tot
      w <- w + rbinom(1, before, rbeta(1, w, tot-w))
      break
    }
    if(w!=0){
      w <- w + rbinom(1, before, rbeta(1, w, tot-w))
    }
    tot <- tot + before
    # Did we draw a white ball?
    if(runif(1) < w/tot){
      # Add a black ball
      tot <- tot + 1
    }else{
      # Add a white ball with probability p_B = p/3
      if(runif(1) < 1/3){
        w <- w+1
      }
      tot <- tot + 1
    }
  }
  print(w)
  return(w/n)
}

## Theoretical CDF
t_cdf <- function(x, p, k){
  if(x == 0){
    (1-p/3)^k
  }else{
    # Integrate theoretical PDF
    (1-p/3)^k + (1-(1-p/3)^k)*(-1+(1-x)^(1+k)+x+k *x)/(k *x)
  }
}




# How many model runs?
n_iters <- 1e6

# Function to plot empirical and theoretical density
get_plot <- function(n = 1e6, p = 1e-6){

  # Set bin width
  bw <- sqrt(p)

  # Set max plot range
  max_range = 50*bw

  # Bin right edges
  bins <- seq(bw, max_range, bw)

  # Simulate RPW model
  sims <- unlist(mclapply(1:n_iters, rpw, p=p, n =n, mc.cores = 12))

  # Get proportion of iterations in each bin (bins being right edges)
  e_cdf <- unlist(mclapply(bins, function(x){mean(sims < x)}, mc.cores = 12))

  # Theoretical CDF
  # Take k = sqrt(1/p)
  t_vals <- sapply(bins, t_cdf, p=p, k = sqrt(1/p))

  df <- data.frame(x = rep(bins, 2), y = c(t_vals, e_cdf), CDF = c(rep("Theoretical", length(bins)), rep("Empirical", length(bins))))

  p <- ggplot(df, aes(x = x, y=y, group = CDF)) +
    geom_line(aes(color = CDF)) +
    ylab("Cumulative Density") +
    xlab("Proportion of White Balls") +
    scale_color_manual(values = c("black", "red")) +
    theme_minimal() +
    theme(legend.position = "top")

  return(p)
}

ns <- c(1e5, 1e6, 1e7)
ps <- c(1e-7, 1e-6, 1e-5)

plots <- list()
for (i in 1:length(ps)) {
  plots <- c(plots, list(get_plot(p = ps[i])))
  print(i)
}
for (i in 1:length(ns)) {
  plots <- c(plots, list(get_plot(n = ns[i])))
  print(i)
}

all_plots <- plot_grid(plots[[1]], plots[[4]], plots[[2]], plots[[5]], plots[[3]], plots[[6]], ncol = 2, labels = c("A","D","B","E","C","F"))

ggsave("./figs/rpw.pdf", width = 6.5, height = 7.5)
ggsave("./figs/rpw.png", width = 6.5, height = 7.5)

