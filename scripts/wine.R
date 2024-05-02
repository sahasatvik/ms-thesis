#!/usr/bin/env Rscript

library(reshape2)
library(cowplot)
library(dplyr)
library(ggplot2)
library(ggrepel)

load("../data/Winedata.rda")

idx  <- (Winedata$ppm < 5.62) & (Winedata$ppm > 5.37)
wave <- Winedata$ppm[idx]
wine <- Winedata$spectra[, idx]
colnames(wine) <- wave

orange <- "#e99352"
purple <- "purple"
red <- "red"
blue <- "blue"
green <- "seagreen"

theme.data <- list(
        scale_color_manual(values = c(orange, purple, red)),
        labs(
                x = "Wavelength",
                y = "Response",
                title = "NMR spectra of wine samples"
        ),
        theme(
                legend.position = "none",
                plot.background = element_rect(fill = "transparent", color = NA),
                panel.background = element_rect(fill = "#f7f7f7", colour = NA),
                text = element_text(family = "Fira Sans"),
                # panel.grid.major = element_blank(),
                # panel.grid.minor = element_blank(),
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

p.data <- ggplot(wine.df, aes(x = variable, y = value / 1e5, group = as.factor(id))) +
        geom_line(linewidth = 0.2, color = "grey") +
        geom_line(data = filter(wine.df, id == 1), linewidth = 0.2, color = purple) +
        geom_line(data = filter(wine.df, id == 37), linewidth = 0.2, color = green) +
        geom_line(data = filter(wine.df, id %in% c(2, 3)), linewidth = 0.2, color = red) +
        geom_line(data = filter(wine.df, id %in% c(23, 35)), linewidth = 0.2, color = blue) +
        geom_text_repel(data = filter(wine.df, (id == 1) & (abs(variable - 5.5906) < 2e-4)), aes(label = as.character(id)), direction = "y", nudge_y = 0.2, color = purple) +
        geom_text_repel(data = filter(wine.df, (id == 37) & (abs(variable - 5.403) < 2e-4)), aes(label = as.character(id)), direction = "x", nudge_x = -0.0005, color = green) +
        geom_text_repel(data = filter(wine.df, (id == 2) & (abs(variable - 5.388) < 2e-4)), aes(label = as.character(id)), direction = "y", nudge_y = -0.1, color = red) +
        geom_text_repel(data = filter(wine.df, (id == 3) & (abs(variable - 5.415) < 2e-4)), aes(label = as.character(id)), direction = "y", nudge_y = -0.1, color = red) +
        geom_text_repel(data = filter(wine.df, (id == 23) & (abs(variable - 5.53) < 2e-4)), aes(label = as.character(id)), direction = "y", nudge_y = -0.1, color = blue) +
        geom_text_repel(data = filter(wine.df, (id == 35) & (abs(variable - 5.54) < 2e-4)), aes(label = as.character(id)), direction = "y", nudge_y = 0.1, color = blue) +
        scale_x_reverse(expand = expansion(mult = 0.01)) +
        theme.data

cairo_pdf("../images/wine.pdf", width = 8, height = 4)

print(p.data)

dev.off()
