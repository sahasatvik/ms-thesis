#!/usr/bin/env Rscript

library(ggplot2)
library(ddalpha)
library(np)

depth <- depth.Mahalanobis
# depth <- depth.spatial
# depth <- depth.halfspace

# source("localdepth.R")

localdepth.reduced <- function(u, X, beta, depth, symm = FALSE) {
        # reflected <- t(apply(X, 1, function(x) { 2 * u - x }))
        reflected <- t(2 * u - t(X))
        augmented <- rbind(X, reflected)
        depths <- depth(X, augmented)
        q <- quantile(depths, prob = 1 - beta)
        idx <- (depths >= q)
        if (symm) return(list(X = rbind(X[idx, ], reflected[idx, ]), idx = idx))
        return(list(X = X[idx, ], idx = idx))
}

localregress <- function(X, Y, points, beta, depth, symm = TRUE) {
        apply(points, 1, function(x) {
                reduced <- localdepth.reduced(x, X, beta = beta, depth = depth, symm = symm)
                # reduced <- reduced[!duplicated(reduced), ]
                depths  <- depth(X[reduced$idx, ], reduced$X)
                weighted.mean(Y[reduced$idx], depths)
        })
}


loocv <- function(X, Y, f) {
        Yhat <- sapply(1:nrow(X), function(i) f(X[-i, ], Y[-i, ], matrix(X[i, ], 1)))
        mean((Y - Yhat)^2, na.rm = TRUE)
}

loocv.apply <- function(X, Y, Z, fz) {
        sapply(Z, function(z) {
                loocv(X, Y, function(X, Y, p) fz(X, Y, p, z))
        })
}

loocv.localregress <- function(X, Y, beta) {
        loocv.apply(X, Y, beta, function(X, Y, p, z)
                localregress(X, Y, p, beta = z, depth = depth, symm = TRUE)
        )
}

loocv.search <- function(X, Y, g, z.lo, z.hi, k = 3) {
        Z <- seq(z.lo, z.hi, length.out = k)
        mse <- g(X, Y, Z)
        return(Z[which.min(mse)])
}


do_regression <- function(
        X, Y,
        X.test, Y.test,
        data, data.test,
        formula,
        title,
        cv.beta.min = 0.10,
        cv.beta.max = 0.20,
        cv.beta.k = 50,
        beta = NULL
) {

        if (is.null(beta)) {
                beta <- loocv.search(X, Y, loocv.localregress, cv.beta.min, cv.beta.max, k = cv.beta.k)
        }
        Yhat.depth <- localregress(X, Y, X.test, beta = beta, depth = depth)
        mse.depth <- mean((Y.test - Yhat.depth)^2)

        cat("\n")
        cat(sprintf("Dataset - %s\n", title))
        cat(sprintf("Local depth based regression:\n"))
        cat(sprintf("    beta       = %f\n", beta))
        cat(sprintf("    MSE        = %f\n", mse.depth))


        bw.lc <- np::npregbw(formula = as.formula(formula), data = data, regtype = "lc")
        bw.ll <- np::npregbw(formula = as.formula(formula), data = data, regtype = "ll")

        fit.lc <- np::npreg(bw.lc, newdata = data.test, y.eval = TRUE)
        fit.ll <- np::npreg(bw.ll, newdata = data.test, y.eval = TRUE)

        mse.lc <- fit.lc$MSE
        mse.ll <- fit.ll$MSE

        cat(sprintf("Local smoothing (locally constant):\n"))
        cat(sprintf("    MSE        = %f\n", mse.lc))
        cat(sprintf("Local smoothing (locally linear):\n"))
        cat(sprintf("    MSE        = %f\n", mse.ll))
        cat("\n")

        print(summary(fit.lc))
        print(summary(fit.ll))

        return(list(beta = beta, Yhat.depth = Yhat.depth))
}


gendata <- function(n) {
        X <- mvrnorm(400, mu = c(0, 0), Sigma = matrix(c(1, 0, 0, 1), 2))
        # X <- matrix(runif(400 * 2, min = -3, max = 3), c(400, 2))
        # Y.real <- apply(X, 1, function(x) sqrt(x[1]^2 + x[2]^2))
        # # Y.real <- apply(X, 1, function(x) x[1] + x[2])
        # Y <- Y.real + rnorm(length(Y.real), sd = 0.5)

        r1 <- runif(0.75 * n, min = 2, max = 4)
        r2 <- sqrt(runif(0.25 * n, min = 0, max = 1.5^2))
        t1 <- runif(0.75 * n, min = 0, max = 2 * pi)
        t2 <- runif(0.25 * n, min = 0, max = 2 * pi)

        X1 <- cbind(r1 * cos(t1), r1 * sin(t1))
        X2 <- cbind(r2 * cos(t2), r2 * sin(t2))

        Y1 <- apply(X1, 1, function(x) 0.5 * sqrt(x[1]^2 + x[2]^2))
        # Y1 <- apply(X1, 1, function(x) abs(x[1] + x[2]^2))
        Y2 <- apply(X2, 1, function(x) 5 - sqrt(x[1]^2 + x[2]^2))
        # Y2 <- apply(X2, 1, function(x) 5 + 3 * abs(x[1] - x[2]) + abs(x[1]^2 - x[2]^2))

        X <- rbind(X1, X2)
        Y.real <- c(Y1, Y2)
        Y <- Y.real + rnorm(length(Y.real), sd = 0.5)

        Y <- as.matrix(Y)
        df.real <- data.frame(x = X[, 1], y = X[, 2], z = Y.real)
        df <- data.frame(x = X[, 1], y = X[, 2], z = Y)

        return(list(X = X, Y = Y, Y.real = Y.real, df = df, df.real = df.real))
}


train <- gendata(250)
test  <- gendata(250)
X <- train$X
Y <- train$Y
df <- train$df
X.test <- test$X
Y.test <- test$Y
df.test <- test$df


scale_limits = c(-0.5, 5.5)
# scale_limits = c(-1, 18)

ggplot(df, aes(x = x, y = y)) +
        geom_point(aes(color = z)) +
        scale_color_gradientn(colours = rainbow(5), limits = scale_limits) +
        ggtitle("Train data")

ggplot(df.test, aes(x = x, y = y)) +
        geom_point(aes(color = z)) +
        scale_color_gradientn(colours = rainbow(5), limits = scale_limits) +
        ggtitle("Test data")

reg <- do_regression(
        X, Y,
        X.test, Y.test,
        df, df.test, "z ~ x + y",
        "Simple",
        # beta = 0.02,
        cv.beta.min = 0.01,
        cv.beta.max = 0.10,
        cv.beta.k = 20,
)

beta <- reg$beta
Yhat.depth <- reg$Yhat.depth

df <- data.frame(x = X.test[, 1], y = X.test[, 2], z = Yhat.depth)
ggplot(df, aes(x = x, y = y)) +
        geom_point(aes(color = z)) +
        scale_color_gradientn(colours = rainbow(5), limits = scale_limits) +
        ggtitle("Depth based regression")

df <- data.frame(x = X.test[, 1], y = X.test[, 2], z = Yhat.depth - Y.test)
ggplot(df, aes(x = x, y = y)) +
        geom_point(aes(color = z)) +
        scale_color_gradient2(low = "red", high = "blue") +
        ggtitle("Residuals $\\hat{y} - y$")

hist(Yhat.depth - Y.test)

# beta <- 0.02

# XX  <- t(matrix(c(0, 0,  1, 0,  0, 2,  -3, 0,  0, -4), 2))
# nbr <- c()
# for (i in 1:nrow(XX)) {
#         x   <- XX[i, ]
#         res <- localdepth_.reduced(x, X, beta = beta, depth = depth, symm = FALSE)$X
#         nbri <- data.frame(x = res[, 1], y = res[, 2], l = i, s = 1)
#         nbri <- rbind(nbri, data.frame(x = x[1], y = x[2], l = 0, s = 2))
#         nbr  <- rbind(nbr, nbri)
# }

# nbr$l <- as.factor(nbr$l)
# nbr$s <- as.factor(nbr$s)

# ggplot(nbr) +
#         geom_point(aes(x = x, y = y), data = df, color = "#cccccc") +
#         geom_point(aes(x = x, y = y, color = l, shape = s)) +
#         scale_shape_manual(values = c(16, 8)) +
#         ggtitle("Neighbourhoods")
