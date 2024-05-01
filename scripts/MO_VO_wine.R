#!/usr/bin/env Rscript

library(reshape2)
library(cowplot)
library(dplyr)
library(ggplot2)
library(ggrepel)

source("MO_VO.R")
load("../data/Winedata.rda")

idx  <- (Winedata$ppm < 5.62) & (Winedata$ppm > 5.37)
wave <- Winedata$ppm[idx]
wine <- Winedata$spectra[, idx]
colnames(wine) <- wave

wine.mv <- MO_VO(
        wine,
        outlyingness = function(y) ((y - median(y)) / mad(y))
)
wine.mv <- mutate(wine.mv, outlier = factor(
        as.numeric(abs(MO) > 1.5) +
        2 * as.numeric(VO > sqrt(0.1)),
        levels = 0:3,
        labels = c("F", "c", "s", "cs")
))

orange <- "#e99352"
purple <- "purple"
red <- "red"
green <- "seagreen"


theme.data <- list(
        scale_color_manual(values = c(orange, red, purple, green)),
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

theme.mv <- list(
        scale_color_manual(values = c(orange, red, purple, green)),
        scale_shape_manual(values = c(3, 4, 4, 4)),
        labs(
                x = "|MO|",
                y = "VO",
                title = "MO-VO diagram"
        ),
        theme(
                legend.position = "none",
                plot.background = element_rect(fill = "transparent", color = NA),
                panel.background = element_rect(fill = "#f7f7f7", colour = NA),
                text = element_text(family = "Fira Sans"),
        )
)

wine.df <- as.data.frame(wine)
wine.df$id <- 1:nrow(wine)
wine.df <- melt(wine.df, id.var = "id")
wine.df$variable <- wave[wine.df$variable]
wine.df <- merge(wine.df, wine.mv)

p.data <- ggplot(wine.df, aes(x = variable, y = value / 1e5, group = as.factor(id))) +
        geom_line(data = filter(wine.df, outlier == "F"),  linewidth = 0.5, color = "grey") +
        geom_line(data = filter(wine.df, outlier == "s"),  linewidth = 0.2, color = purple, alpha = 0.8) +
        geom_line(data = filter(wine.df, outlier == "c"),  linewidth = 0.2, color = red, alpha = 0.8) +
        geom_line(data = filter(wine.df, outlier == "cs"), linewidth = 0.2, color = green, alpha = 0.8) +
        scale_x_reverse(expand = expansion(mult = 0.01)) +
        theme.data

p.mv <- ggplot(wine.mv, aes(x = abs(MO), y = VO, shape = outlier, color = outlier)) +
        geom_point() +
        geom_text_repel(aes(label = ifelse(outlier != "F", as.character(id), "")), direction = "both") +
        scale_y_log10() +
        theme.mv

png("../images/MO_VO_wine.png", width = 8, height = 4, units = "in", res = 160)

plot_grid(p.data, p.mv)

dev.off()
