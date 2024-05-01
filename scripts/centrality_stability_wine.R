#!/usr/bin/env Rscript

library(reshape2)
library(cowplot)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggrepel)

source("centrality_stability.R")
load("../data/Winedata.rda")

idx  <- (Winedata$ppm < 5.62) & (Winedata$ppm > 5.37)
wave <- Winedata$ppm[idx]
wine <- Winedata$spectra[, idx]
colnames(wine) <- wave

cs <- centrality_stability(
        wine,
        outlyingness = function(y) abs((y - median(y)) / mad(y))
)
wine.cs <- cs$df
wine.cs <- mutate(wine.cs, outlier = factor(
        as.numeric(centrality > 0.55) +
        2 * as.numeric(stability > 0.15),
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

theme.cs <- list(
        scale_color_manual(values = c(orange, red, purple, green)),
        scale_shape_manual(values = c(3, 4, 4, 4)),
        xlim(0.1, 0.7),
        labs(
                x = "Centrality deviation",
                y = "Stability deviation",
                title = "Centrality-Stability diagram"
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
wine.df <- merge(wine.df, wine.cs)

p.data <- ggplot(wine.df, aes(x = variable, y = value / 1e5, group = as.factor(id))) +
        geom_line(data = filter(wine.df, outlier == "F"),  linewidth = 0.5, color = "grey") +
        geom_line(data = filter(wine.df, outlier == "s"),  linewidth = 0.2, color = purple, alpha = 0.8) +
        geom_line(data = filter(wine.df, outlier == "c"),  linewidth = 0.2, color = red, alpha = 0.8) +
        geom_line(data = filter(wine.df, outlier == "cs"), linewidth = 0.2, color = green, alpha = 0.8) +
        scale_x_reverse(expand = expansion(mult = 0.01)) +
        theme.data

p.cs <- ggplot(wine.cs, aes(x = centrality, y = stability, shape = outlier, color = outlier)) +
        geom_point() +
        geom_text_repel(aes(label = ifelse(outlier != "F", as.character(id), "")), direction = "both") +
        scale_y_log10() +
        theme.cs

png("../images/centrality_stability_wine.png", width = 8, height = 4, units = "in", res = 160)

plot_grid(p.data, p.cs)

dev.off()



O.df <- as.data.frame(cs$O)
O.df[, "id"] <- as.factor(1:nrow(wine))
O.df <- pivot_longer(O.df, -id)
O.df$name <- as.double(O.df$name)


png("../images/outlyingness_heatmap_wine.png", width = 8, height = 4, units = "in", res = 160)

theme.hm <- list(
        labs(
                x = "Wavelength",
                y = NULL,
                title = "Outlyingness heatmap"
        ),
        theme(
                legend.position = "none",
                plot.background = element_rect(fill = "transparent", color = NA),
                panel.background = element_rect(fill = "#f7f7f7", colour = NA),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                text = element_text(family = "Fira Sans"),
        ),
        # scale_fill_gradient2(high = orange),
        scale_fill_gradientn(colors = c("white", orange, purple, "black")),
        # scale_fill_viridis_c(option = "inferno", direction = -1),
        scale_x_reverse(expand = c(0, 0)),
        scale_y_discrete(guide = guide_axis(n.dodge = 2), expand = c(0, 0))
)

p.hm <- ggplot(O.df, aes(x = name, y = id)) +
        geom_tile(aes(fill = value)) +
        theme.hm

p.hm2 <- ggplot(mutate(O.df, value = ifelse(id == 37, 0, value)), aes(x = name, y = id)) +
        geom_tile(aes(fill = value)) +
        theme.hm +
        ggtitle("(Zeroed curve #37)")

plot_grid(p.hm, p.hm2)

dev.off()
