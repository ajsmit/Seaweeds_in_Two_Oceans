### Spatial analysis/MEM functions

# Made by AJ Smit (2015)

source("theme/themes.R")

ggplotMEM <- function(obj) {
	thelen <- length(obj$obs)
	themax <- which.max(obj$obs)
	thecol <- rep("b", thelen)
	thecol[themax] <- "a"
	obj.df <- data.frame(x = 1:thelen,
	                     y = obj$obs,
	                     qu95 = apply(obj$sim, 2, quantile, 0.95),
	                     cols = thecol)
	ggplot(data = obj.df, aes(x = x, y = y)) +
	geom_bar(stat = "identity", aes(col = cols, fill = cols), show.legend = FALSE) +
	# geom_point(aes(y = qu95), shape = 4, col = "black", size = 0.8) +
	geom_line(aes(y = qu95), col = "black", linetype = "solid", size = 0.5) +
	scale_x_continuous(breaks = seq(1, thelen, 1)) +
	xlab("MEMs") + ylab(expression(R^2)) +
	theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.6),
	      axis.text = element_text(size = 8),
	      plot.title = element_text(hjust = 0, size = 10))
}

## ggplot2 function for the South African Coast:
ggmap <- function() {
	sa_lats <- c(-35.5, -26); sa_lons <- c(14, 34)
	load("/Users/ajsmit/Dropbox/repos/tempssa_v3.0/graph/sa_shore.Rdata")
	sa_shore <- fortify(sa_shore)
	sa_plot <- ggplot() + theme_grey() +
	geom_polygon(data = sa_shore, aes(x = X, y = Y, group = PID),
	             show.legend = FALSE, fill = "#F9FAEC", col = "#F9FAEC", size = 0.2) +
	coord_map(xlim = sa_lons, ylim = sa_lats) +
	  coord_fixed(ratio = 1, expand = TRUE) +
	  scale_x_continuous(breaks = seq(15, 35, 5),
	                     labels = scales::unit_format("°E", sep = "")) +
	  scale_y_continuous(breaks = seq(-35, -25, 5),
	                     labels = c("35°S", "30°S", "25°S")) +
	  xlab(NULL) + ylab(NULL) +
	  theme(strip.background = element_blank(),
	        strip.text = element_text(size = 10, hjust = 0))
	return(sa_plot)
}

## Function to plot the background map in base graphics:
plotmap <- function() {
	par(mar = c(5.1, 4.1, 4.1, 2.1))
    plot(shore, type = "l", lwd = 0.1, add = FALSE, col = "#dfd7d2",
         border = "#dfd7d2")
    }
