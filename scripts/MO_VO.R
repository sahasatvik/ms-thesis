#!/usr/bin/env Rscript

MO_VO <- function(
        X,
        X.ref = X,
        outlyingness = function(x, x.ref) (x - median(x.ref)) / mad(x.ref)
) {
        O <- sapply(1:ncol(X), function(i) outlyingness(X[, i], X.ref[, i]))
        MO <- rowMeans(O)
        VO <- rowMeans((O - MO)^2)
        return(data.frame(
                id = 1:nrow(X),
                MO = MO,
                VO = VO
        ))
}
