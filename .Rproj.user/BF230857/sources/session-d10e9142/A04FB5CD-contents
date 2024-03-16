library(ggplot2)

t <- c(5,5,5,5,5, 4,4, 3, 2, 1)
thetas <- c(1,2,3,4,5, 1.5, 3.5, 4.25, 2.875, 2.875)

h <- c(6,6,7,7,8,9,8,9,10,NA)
n <- length(h)

Lineage <- c("Variant B", "Variant B", "Variant A", "Variant A", "Variant A")

# vertical segments
xs <- c()
ystart <- c()
yend <- c()
for (i in 1:n) {
  kids <- which(h == i)
  if(length(kids) > 0){
    xs <- c(xs, t[i])
    ystart <- c(ystart, min(thetas[kids]))
    yend <- c(yend, max(thetas[kids]))
  }
}

phy <- ggplot() +
  geom_segment(mapping = aes(x = t[h], xend = t, y = thetas, yend = thetas), linewidth = 0.5) +
  geom_segment(mapping = aes(x = xs, xend = xs, y = ystart, yend = yend), linewidth = 0.5) +
  geom_point(mapping= aes(x = t[1:5], y = thetas[1:5], color = Lineage), size = 2) +
  geom_point(mapping = aes(x = 3.25, y = 1.5), shape = 4, size = 4, color = "red") +
  scale_color_manual(values = c("blue", "red"), name = "Lineage") +
  xlab("Evolutionary Time") +
  scale_y_continuous(breaks = NULL) +
  theme_minimal() +
  theme(axis.title.y = element_blank(), legend.position = 'top')

phy

ggsave(file = "./figs/example_phylogeny.pdf", width = 5, height = 3)
ggsave(file = "./figs/example_phylogeny.png", width = 5, height = 3)
