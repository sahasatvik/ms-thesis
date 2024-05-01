#!/usr/bin/env Rscript

library(carData)
data(Anscombe)

source("contours.R")

df <- data.frame(X = Anscombe$education, Y = Anscombe$income)

# Z <- mvrnorm(100, mu = c(0, 0), Sigma = matrix(c(1, 0.9, 0.9, 1), 2))
# df <- data.frame(X = Z[, 1], Y = Z[, 2])

grid.n <- 300

p.mahalanobis <- show_contour(
        df,
        grid.n = grid.n,
        title = "Mahalanobis",
        depth = depth.Mahalanobis,
)

p.spatial <- show_contour(
        df,
        grid.n = grid.n,
        title = "Spatial (Affine invariant)",
        depth = function(U, X) depth.spatial(U, X, mah.estimate = "moment"),
)

p.halfspace <- show_contour(
        df,
        grid.n = grid.n,
        title = "Halfspace",
        depth = function(U, X) depth.halfspace(U, X, exact = TRUE),
)

p.projection <- show_contour(
        df,
        grid.n = grid.n,
        title = "Projection",
        depth = depth.projection,
)

p.simplicial <- show_contour(
        df,
        grid.n = grid.n,
        title = "Simplicial",
        depth = function(U, X) depth.simplicial(U, X, exact = TRUE),
)

p.oja <- show_contour(
        df,
        grid.n = grid.n,
        title = "Oja (Affine invariant)",
        depth = function(U, X) depth.simplicialVolume(U, X, exact = TRUE, mah.estimate = "moment"),
)


cairo_pdf("../images/contours.pdf", onefile = TRUE, width = 8, height = 12)

plot_grid(
        p.halfspace,
        p.mahalanobis,
        p.spatial,
        p.projection,
        p.simplicial,
        p.oja,
        ncol = 2
)

dev.off()
