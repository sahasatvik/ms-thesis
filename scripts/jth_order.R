#!/usr/bin/env Rscript

depth.FM <- function(X, Xref, depth = depth.Mahalanobis) {
        T <- ncol(Xref)
        rowMeans(sapply(1:T, function(i) { depth(matrix(X[, i]), matrix(Xref[, i])) }))
}

depth.FMJ <- function(X, Xref, J = 2, n = 1000, depth = depth.Mahalanobis, diag = FALSE) {
        if (J == 1) { return(depth.FM(X, Xref, depth = depth)) }
        T <- ncol(Xref)
        S <- replicate(n, sample(1:T, J, replace = diag))
        rowMeans(apply(S, 2, function(idx) { depth(X[, idx], Xref[, idx]) }))
}

depth.FMJ.outliers <- function(
        X, Xref,
        J = 2,
        n = 1000,
        depth.1 = depth.Mahalanobis,
        depth.j = depth.Mahalanobis,
        diag = FALSE,
        F = 1.5
) {
        depths <- matrix(depth.FM(X, Xref, depth = depth.1))
        for (j in 2:J) {
                depths <- cbind(depths, matrix(depth.FMJ(X, Xref, J = j, n = n, depth = depth.j, diag = diag)))
        }
        S <- t(apply(depths, 1, function(r) log(r / r[J])))
        S <- matrix(S[, -J], nrow(X))
        cutoffs <- apply(S, 2, function(s) quantile(s, probs = 0.75) + F * IQR(s))
        outliers <- apply((t(S) > cutoffs) * (1:(J - 1)), 2, function(x) min(x[x != 0]))
        return(list(depths = depths, S = S, cutoffs = cutoffs, outliers = outliers))
}
