#!/usr/bin/env Rscript

library(ddalpha)
library(ggplot2)
library(cowplot)
library(tibble)
library(dplyr)
library(extrafont)

orange <- "#e99352"
purple <- "purple"

N <- 150

cairo_pdf("../images/ddplots.pdf", onefile = TRUE, width = 8, height = 4)

theme.data <- list(
        geom_point(size = 1.5, alpha = 1),
        scale_color_manual(values = c(orange, purple)),
        scale_shape_manual(values = c(3, 4)),
        labs(
                x = "x",
                y = "y",
                title = "Data"
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

theme.dd <- list(
        geom_point(size = 1.5, alpha = 1),
        geom_abline(aes(intercept = 0, slope = 1)),
        scale_color_manual(values = c(orange, purple)),
        scale_shape_manual(values = c(3, 4)),
        xlim(0, 1),
        ylim(0, 1),
        labs(
                x = "Depth (Orange)",
                y = "Depth (Purple)",
                title = "DD Plot"
        ),
        theme(
                legend.position = "none",
                plot.background = element_rect(fill = "transparent", color = NA),
                panel.background = element_rect(fill = "#f7f7f7", colour = NA),
                text = element_text(family = "Fira Sans"),
                axis.title.x = element_text(colour = orange),
                axis.title.y = element_text(colour = purple),
        )
)

show_ddplot <- function(A, B) {
        dfA <- tibble(group = "A", x = A[, 1], y = A[, 2])
        dfB <- tibble(group = "B", x = B[, 1], y = B[, 2])
        X <- rbind(dfA, dfB)
        X$group <- factor(X$group, levels = c("A", "B"))

        X <- X %>%
                mutate(depth_A = depth.spatial(cbind(x, y), A)) %>%
                mutate(depth_B = depth.spatial(cbind(x, y), B))

        X.plot <- ggplot(
                X[sample(nrow(X)), ],
                aes(x = x, y = y, colour = group, shape = group)
        ) + theme.data

        X.dd <- ggplot(
                X[sample(nrow(X)), ],
                aes(x = depth_A, y = depth_B, colour = group, shape = group)
        ) + theme.dd

        plot_grid(X.plot, X.dd)
}


# No difference
A <- mvrnorm(N, mu = c(0, 0), Sigma = matrix(c(1, 0, 0, 1), 2))
B <- mvrnorm(N, mu = c(0, 0), Sigma = matrix(c(1, 0, 0, 1), 2))
show_ddplot(A, B)


# Location difference
A <- mvrnorm(N, mu = c(-1.0, 0), Sigma = matrix(c(1, 0, 0, 1), 2))
B <- mvrnorm(N, mu = c( 1.0, 0), Sigma = matrix(c(1, 0, 0, 1), 2))
show_ddplot(A, B)


# Scale difference
A <- mvrnorm(N, mu = c(0, 0), Sigma = matrix(c(1, 0, 0, 1), 2))
B <- mvrnorm(N, mu = c(0, 0), Sigma = matrix(c(4, 0, 0, 4), 2))
show_ddplot(A, B)


# Scale difference (shape)
A <- mvrnorm(N, mu = c(0, 0), Sigma = matrix(c(1, -0.5, -0.5, 1), 2))
B <- mvrnorm(N, mu = c(0, 0), Sigma = matrix(c(1, 0.5, 0.5, 1), 2))
show_ddplot(A, B)


# Location and scale difference
A <- mvrnorm(N, mu = c(-1.0, 0), Sigma = matrix(c(1, 0, 0, 1), 2))
B <- mvrnorm(N, mu = c( 1.0, 0), Sigma = matrix(c(9, 0, 0, 9), 2))
show_ddplot(A, B)


dev.off()
