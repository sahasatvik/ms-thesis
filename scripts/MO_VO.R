#!/usr/bin/env Rscript

MO_VO <- function(
        X,
        outlyingness = function(col) ((col - mean(col)) / sd(col))^2
) {
        O <- apply(X, 2, outlyingness)
        MO <- rowMeans(O)
        VO <- rowMeans((O - MO)^2)
        return(data.frame(
                id = 1:nrow(X),
                MO = MO,
                VO = VO
        ))
}
