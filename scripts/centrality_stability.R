#!/usr/bin/env Rscript

centrality_stability <- function(
        X,
        outlyingness = function(col) ((col - mean(col)) / sd(col))^2
) {
        O <- apply(X, 2, outlyingness)
        a <- rowMeans(1 / (1 + O))
        b <- rowMeans(1 + O)
        return(list(
                O = O,
                df = data.frame(
                        id = 1:nrow(X),
                        centrality = 1 - a,
                        stability  = b - 1 / a
                ))
        )
}
