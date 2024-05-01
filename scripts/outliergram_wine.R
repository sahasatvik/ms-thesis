#!/usr/bin/env Rscript

library(cowplot)
library(dplyr)
library(ggplot2)

source("outliergram.R")
load("../data/Winedata.rda")

idx  <- (Winedata$ppm < 5.62) & (Winedata$ppm > 5.37)
wave <- Winedata$ppm[idx]
wine <- Winedata$spectra[, idx]
colnames(wine) <- wave

og <- outliergram(wine)

orange <- "#e99352"
purple <- "purple"
red <- "red"

theme.data <- list(
        scale_color_manual(values = c(orange, purple, red)),
        labs(
                x = "Wavelength",
                y = "Response",
                title = "Data"
        ),
        theme(
                legend.position = "none",
                plot.background = element_rect(fill = "transparent", color = NA),
                panel.background = element_rect(fill = "#f7f7f7", colour = NA),
                text = element_text(family = "Fira Sans"),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                # axis.ticks.x = element_blank(),
                # axis.ticks.y = element_blank(),
                # axis.text.x = element_blank(),
                # axis.text.y = element_blank(),
        )
)

wine.df <- as.data.frame(wine)
wine.df$id <- 1:nrow(wine)
wine.df <- melt(wine.df, id.var = "id")
wine.df$variable <- wave[wine.df$variable]
wine.df <- merge(wine.df, og$df[, c("id", "outlier")])

p.data <- ggplot(wine.df, aes(x = variable, y = value / 1e5, group = as.factor(id))) +
        geom_line(data = filter(wine.df, outlier == "F"),   linewidth = 0.5, color = "grey") +
        geom_line(data = filter(wine.df, outlier == "d"),   linewidth = 0.2, color = purple, alpha = 0.6) +
        geom_line(data = filter(wine.df, outlier == "mbd"), linewidth = 0.2, color = red, alpha = 0.8) +
        scale_x_reverse() +
        theme.data

png("../images/outliergram_wine.png", width = 8, height = 4, units = "in", res = 160)

plot_grid(p.data, og$plot)

dev.off()
