#plot_yeast_erc_old_and_new.R was used to compare tables of ERC values (compare to Nathan Clarks ERC values online)
#Brandon Jerngian

library("ggplot2")
library("scales")
library("grid")
library("gridExtra")
library("ggExtra")
library("reshape2")
# setwd("C:/Users/brand/Desktop/streamlined_ERC_pfam/ERC/")
print("loading All")
# yeast_both = read.csv("yeast_points_all.csv")
# print("loading Old")
# yeast_old = read.csv("yeast_points_old.csv")
print("loading New")
# yeast_new = read.csv("yeast_points_new_fixed.csv")

print("Plotting Histogram")
# p <- ggplot(yeast_both[sample(nrow(yeast_both), 10000),], aes(x = old_data, y = new_data)) +
# p <- ggplot(yeast_both[complete.cases(yeast_both), ], aes(x = old_data, y = new_data)) +
#   geom_point(shape = 1, color = "black") +
#   labs(x = "Nathan's ERC Values", y = "Our ERC Values", title = "ERC Correlation for Yeast Genes") +
#   theme(plot.title = element_text(size = 22, face="bold")) +
#   scale_x_continuous(breaks = pretty_breaks(n = 10)) +
#   scale_y_continuous(breaks = pretty_breaks(n = 10))
# 
# marg = ggMarginal(p, type = "histogram", binwidth = 0.05)
# 
# ggsave("erc_marg_web2.png", plot = marg, width = 10, height = 10, units = 'in')

# print("Plotting Old Scatter")
# h1 <- ggplot(yeast_old, aes(old_data)) +
#   geom_histogram(binwidth = 0.04) +
#   labs(x = "Nathan's Yeast ERC Values", y = "Count", title = "Nathan's Yeast Gene ERC Values") +
#   theme(plot.title = element_text(size = 14, face="bold")) +
#   coord_cartesian(xlim = c(-1, 1)) +
#   scale_x_continuous(breaks = pretty_breaks(n = 10)) +
#   scale_y_continuous(breaks = pretty_breaks(n = 10))
# 
small <- list()
for (ii in 1:500000){
  if(!is.na(yeast_new$new_data[ii])){
    if(yeast_new$new_data[ii] < 0.1){
    small <- c(small, yeast_new$new_data[ii])
    print( ii )
    }
  }
}
small <- melt(small$value)

print("Plotting New Scatter")
# h2 <- ggplot(, aes(data.frame(yeast_new[seq(500000),]))) +
h2 <- ggplot(, aes(small)) +
  geom_histogram(binwidth = 0.005) +
  labs(x = "Yeast Domain P Values\n\nWith 13 Species:\nS cerevisiae, K lactis, A gossypii,
       L elongisporus, D hansenii, C glabrata, C albicans,
       C dubliniensis, C lusitaniae, C tropicalis, C guilliermondii, K thermotolerans, S stipitis", y = "Count", title = "Yeast Domain P Values") +
  theme(plot.title = element_text(size = 14, face="bold")) 
  # coord_cartesian(xlim = c(0, 0.1)) +
  # scale_x_continuous(breaks = pretty_breaks(n = 10)) +
  # scale_y_continuous(breaks = pretty_breaks(n = 10))

h2
# 
# 
# gr = grid.arrange(h1, h2)
# 
ggsave("erc_hist_all_fixed.png", plot = h2)
