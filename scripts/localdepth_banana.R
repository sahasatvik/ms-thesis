#!/usr/bin/env Rscript

library(cowplot)

source("contours.R")
source("../../datadepth/localdepth.R")


Z <- mvrnorm(500, mu = c(0, 0), Sigma = matrix(c(1, 0.9, 0.9, 1), 2))
a <- 1
b <- 1
df <- data.frame(
        X = a * Z[, 1],
        Y = Z[, 2] / a - b * ((a * Z[, 1])^2 + a^2)
)

spatial <- show_contour(
        df,
        grid.n = 50,
        depth = depth.spatial,
        title = "Spatial depth"
)

spatial.local <- show_contour(
        df,
        grid.n = 50,
        depth = function(U, X) localdepth(
                as.matrix(U),
                as.matrix(X),
                beta = 0.2,
                depth = depth.spatial
        ),
        title = "Local spatial depth, β = 0.2"
)


cairo_pdf("../images/localdepth_banana.pdf", onefile = TRUE, width = 8, height = 4)

plot_grid(spatial, spatial.local)

dev.off()

