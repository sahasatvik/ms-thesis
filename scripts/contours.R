#!/usr/bin/env Rscript

library(ddalpha)
library(ggplot2)
library(tibble)
library(dplyr)
library(extrafont)

orange <- "#e99352"
purple <- "purple"

theme.data <- list(
        labs(
                x = "x",
                y = "y",
                color = "depth"
        ),
        theme(
                legend.position = "none",
                plot.background = element_rect(fill = "transparent", color = NA),
                panel.background = element_rect(fill = "#f7f7f7", colour = NA),
                text = element_text(family = "Fira Sans"),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.ticks.x = element_blank(),
                axis.ticks.y = element_blank(),
                axis.text.x = element_blank(),
                axis.text.y = element_blank(),
        )
)


show_contour <- function(
        df,
        grid.n = 50,
        grid.breaks = seq(0, 1, length.out = 16)^2,
        depth = depth.spatial,
        title = ""
) {
        delta.X <- diff(range(df$X)) * 0.1
        delta.Y <- diff(range(df$Y)) * 0.1
        grid.X <- seq(min(df$X) - delta.X, max(df$X) + delta.X, length.out = grid.n)
        grid.Y <- seq(min(df$Y) - delta.Y, max(df$Y) + delta.Y, length.out = grid.n)
        grid.D <- expand.grid(X = grid.X, Y = grid.Y) %>%
                mutate(D = depth(cbind(X, Y), df))


        ggplot() +
                # geom_contour(data = grid.D, aes(x = X, y = Y, z = D, colour = after_stat(level)), linewidth = 1, breaks = grid.breaks) +
                # geom_raster(data = grid.D, aes(x = X, y = Y, fill = log(D))) +
                geom_point(data = df, aes(x = X, y = Y), colour = purple, size = 1, shape = 3) +
                geom_contour(data = grid.D, aes(x = X, y = Y, z = D, colour = after_stat(level)), linewidth = 0.6, breaks = grid.breaks) +
                scale_colour_viridis_c(option = "inferno", direction = -1) +
                xlim(min(grid.X), max(grid.X)) +
                ylim(min(grid.Y), max(grid.Y)) +
                theme.data +
                labs(
                        title = title,
                )
}

