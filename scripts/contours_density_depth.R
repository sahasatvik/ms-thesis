#!/usr/bin/env Rscript

library(ddalpha)
library(ggplot2)
library(tibble)
library(dplyr)
library(extrafont)
library(cowplot)

# source("contours.R")

orange <- "#e99352"
purple <- "purple"

theme.data <- list(
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
        ),
        coord_cartesian(xlim = c(-3, 3), ylim = c(-3, 3), expand = FALSE)
        # scale_x_continuous(expand = expansion(mult = 0.02)),
        # scale_y_continuous(expand = expansion(mult = 0.02))
)



A <- mvrnorm(50, mu = c(0, 0), Sigma = matrix(c(1, 0, 0, 1), 2))
df <- data.frame(X = A[, 1], Y = A[, 2])
# df <- data.frame(X = Anscombe$education, Y = Anscombe$income)

grid.n <- 100
depth <- depth.spatial

# delta.X <- diff(range(df$X)) * 0.1
# delta.Y <- diff(range(df$Y)) * 0.1
# grid.X <- seq(min(df$X) - delta.X, max(df$X) + delta.X, length.out = grid.n)
# grid.Y <- seq(min(df$Y) - delta.Y, max(df$Y) + delta.Y, length.out = grid.n)
grid.X <- seq(-3, 3, length.out = grid.n)
grid.Y <- seq(-3, 3, length.out = grid.n)
grid.D <- expand.grid(X = grid.X, Y = grid.Y) %>%
        mutate(D = depth(cbind(X, Y), df))

grid.breaks <- seq(0, 1, length.out = 16)^2

p.depth <- ggplot() +
        geom_point(data = df, aes(x = X, y = Y), colour = purple, size = 1, shape = 3) +
        # geom_contour(data = grid.D, aes(x = X, y = Y, z = D, colour = after_stat(level)), linewidth = 0.6, breaks = grid.breaks) +
        geom_contour(data = grid.D, aes(x = X, y = Y, z = D, colour = after_stat(level)), linewidth = 0.6) +
        scale_colour_viridis_c(option = "inferno", direction = -1) +
        xlim(min(grid.X), max(grid.X)) +
        ylim(min(grid.Y), max(grid.Y)) +
        theme.data +
        labs(
                x = NULL,
                y = NULL,
                title = "Spatial depth contours",
        )


p.density <- ggplot() +
        geom_point(data = df, aes(x = X, y = Y), colour = purple, size = 1, shape = 3) +
        geom_density_2d(data = df, aes(x = X, y = Y, colour = after_stat(level)), linewidth = 0.6) +
        scale_colour_viridis_c(option = "inferno", direction = -1) +
        xlim(min(grid.X), max(grid.X)) +
        ylim(min(grid.Y), max(grid.Y)) +
        theme.data +
        labs(
                x = NULL,
                y = NULL,
                title = "Density contours",
        )


cairo_pdf("../images/contours_density_depth.pdf", onefile = TRUE, width = 8, height = 4)

plot_grid(p.density, p.depth)

dev.off()
