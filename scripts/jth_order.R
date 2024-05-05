#!/usr/bin/env Rscript

depth.FM <- function(X, Xref, depth = depth.Mahalanobis) {
        T <- ncol(Xref)
        rowMeans(sapply(1:T, function(i) { depth(matrix(X[, i]), matrix(Xref[, i])) }))
}

depth.FMJ <- function(X, Xref, J = 2, n = 1000, depth = depth.Mahalanobis) {
        if (J == 1) { return(depth.FM(X, Xref, depth = depth)) }
        T <- ncol(Xref)
        S <- replicate(n, sample(1:T, J, replace = FALSE))
        rowMeans(apply(S, 2, function(idx) { depth(X[, idx], Xref[, idx]) }))
}
