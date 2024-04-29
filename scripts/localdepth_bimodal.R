#!/usr/bin/env Rscript

source("contours.R")
source("../../datadepth/localdepth.R")


cairo_pdf("../images/localdepth_bimodal.pdf", onefile = TRUE, width = 8, height = 4)

Z1 <- mvrnorm(200, mu = c(2, 2), Sigma = matrix(c(1, 0, 0, 1), 2))
Z2 <- mvrnorm(200, mu = c(-2, -2), Sigma = matrix(c(1, 0, 0, 1), 2))
Z  <- rbind(Z1, Z2)
df <- data.frame(
        X = Z[, 1],
        Y = Z[, 2]
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


plot_grid(spatial, spatial.local)

dev.off()

