  # Set up the fonts:
  # arial <- c("~/Dropbox/R/fonts/Arial.afm",
  #            "~/Dropbox/R/fonts/Arial_Bold.afm",
  #            "~/Dropbox/R/fonts/Arial_Italic.afm",
  #            "~/Dropbox/R/fonts/Arial_Bold_Italic.afm")

  require(grid); require(ggplot2)
  # require(showtext)

  theme_set(theme_bw())
  bw_update <- theme_bw() +
  theme(plot.background = element_blank(),
        panel.background = element_rect(colour = NA, fill = NA),
        panel.border = element_rect(colour = "black",
                                    fill = NA, size = 0.8),
        panel.grid.minor = element_line(colour = NA),
        panel.grid.major = element_line(colour = "black",
                                        size = 0.2,
                                        linetype = "dotted"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        axis.ticks = element_line(size = 0.5),
        axis.ticks.length = unit(2, "mm"),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.key = element_blank(),
        legend.key.height = unit(.22, "cm"),
        #legend.key.width = unit(0.30, "cm"),
        legend.background = element_blank(),
        plot.title = element_text(size = 10, hjust = 0),
        strip.background = element_rect(colour = "white", fill = "white"),
        strip.text = element_text(size = 8))
  theme_set(bw_update) # Set theme