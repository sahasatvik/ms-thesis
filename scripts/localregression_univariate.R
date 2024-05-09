#!/usr/bin/env Rscript

library(ddalpha)
library(carData)
library(cowplot)
library(ggplot2)
library(extrafont)

# depth <- function(X, Y) depth.Mahalanobis(X, Y)
# depth <- function(X, Y) apply(X, 1, function(x) 2 * min(mean(x >= Y), mean(x <= Y)))
depth <- function(X, Y) apply(X, 1, function(x) (1 - abs(mean(sign(x - Y)))))


localdepth.1.neighbours <- function(u, X, beta, depth, symm = FALSE) {
        # reflected <- sapply(X, function(x) { 2 * u - x })
        reflected <- 2 * u - X
        augmented <- c(X, reflected)
        depths <- depth(matrix(X), matrix(augmented))
        q <- quantile(depths, prob = 1 - beta)
        if (symm) return(augmented[depths >= q])
        return(X[depths >= q])
}


localregress <- function(X, Y, points, beta, depth, symm = TRUE) {
        sapply(points, function(x) {
                neighbours <- localdepth.1.neighbours(x, X, beta = beta, depth = depth, symm = symm)
                depths <- depth(matrix(X), matrix(neighbours))
                weighted.mean(Y, depths)
        })
}


loocv <- function(X, Y, f) {
        y <- sapply(1:length(X), function(i) f(X[-i], Y[-i], X[i]))
        mean((Y - y)^2, na.rm = TRUE)
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

loocv.ksmooth <- function(X, Y, h) {
        loocv.apply(X, Y, h, function(X, Y, p, z)
                ksmooth(X, Y, x.points = p, bandwidth = z, kernel = "norm")$y
        )
}

loocv.loess1 <- function(X, Y, h) {
        loocv.apply(X, Y, h, function(X, Y, p, z) {
                lo <- loess(y ~ x, data.frame(x = X, y = Y), span = z, degree = 1)
                predict(lo, data.frame(x = p))
        })
}

loocv.loess2 <- function(X, Y, h) {
        loocv.apply(X, Y, h, function(X, Y, p, z) {
                lo <- loess(y ~ x, data.frame(x = X, y = Y), span = z, degree = 2)
                predict(lo, data.frame(x = p))
        })
}

loocv.search <- function(X, Y, g, z.lo, z.hi, k = 3) {
        Z <- seq(z.lo, z.hi, length.out = k)
        mse <- g(X, Y, Z)
        return(Z[which.min(mse)])
}


orange <- "#e99352"
purple <- "purple"

theme.data <- list(
        # scale_color_manual(values = c(orange, purple)),
        # scale_fill_gradient(low = orange, high = purple, na.value = "transparent"),
        # scale_shape_manual(values = c(3, 4)),
        theme(
                legend.position = "none",
                plot.background = element_rect(fill = "transparent", color = NA),
                panel.background = element_rect(fill = "#f7f7f7", colour = NA),
                text = element_text(family = "Fira Sans"),
                # panel.grid.major = element_blank(),
                # panel.grid.minor = element_blank(),
                # axis.ticks.x = element_blank(),
                # axis.ticks.y = element_blank(),
                # axis.text.x = element_blank(),
                # axis.text.y = element_blank(),
        ),
        scale_x_continuous(expand = expansion(mult = 0.01))
        # scale_y_continuous(expand = expansion(mult = 0.0))
)


do.regression <- function(
        X, Y,
        title,
        xlab = "x",
        ylab = "y",
        n.points = 100,
        cv.beta.min = 0.01,
        cv.beta.max = 0.20,
        cv.beta.k = 50,
        cv.bw.min = 0.1,
        cv.bw.max = 5,
        cv.bw.k = 50,
        cv.span1.min = 0.1,
        cv.span1.max = 1.0,
        cv.span1.k = 50,
        cv.span2.min = 0.1,
        cv.span2.max = 1.0,
        cv.span2.k = 50
) {

        grid <- seq(min(X), max(X), length.out = n.points)

        beta <- loocv.search(
                X, Y,
                loocv.localregress,
                cv.beta.min,
                cv.beta.max,
                k = cv.beta.k
        )

        score.depth <- loocv.localregress(X, Y, beta)
        y.depth <- localregress(X, Y, grid, beta = beta, depth = depth, symm = TRUE)
        Yhat.depth <- localregress(X, Y, X, beta = beta, depth = depth, symm = TRUE)
        mse.depth <- mean((Y - Yhat.depth)^2)

        cat("\n")
        cat(sprintf("Dataset - %s\n", title))
        cat(sprintf("Local depth based regression:\n"))
        cat(sprintf("    beta       = %f\n", beta))
        cat(sprintf("    CV score   = %f\n", score.depth))
        cat(sprintf("    MSE        = %f\n", mse.depth))


        bw <- loocv.search(
                X, Y,
                loocv.ksmooth,
                cv.bw.min,
                cv.bw.max,
                k = cv.bw.k
        )

        score.kernel <- loocv.ksmooth(X, Y, bw)
        y.kernel <- ksmooth(X, Y, x.points = grid, bandwidth = bw, kernel = "norm")$y
        Yhat.kernel <- ksmooth(X, Y, x.points = X, bandwidth = bw, kernel = "norm")$y
        mse.kernel <- mean((Y - Yhat.kernel)^2)

        cat(sprintf("Kernel based regression:\n"))
        cat(sprintf("    bandwidth  = %f\n", bw))
        cat(sprintf("    CV score   = %f\n", score.kernel))
        cat(sprintf("    MSE        = %f\n", mse.kernel))


        span1 <- loocv.search(
                X, Y,
                loocv.loess1,
                cv.span1.min,
                cv.span1.max,
                k = cv.span1.k
        )

        score.loess1 <- loocv.loess1(X, Y, span1)
        lo1 <- loess(y ~ x, data.frame(x = X, y = Y), span = span1, degree = 1)
        y.loess1 <- predict(lo1, data.frame(x = grid))
        Yhat.loess1 <- predict(lo1, data.frame(x = X))
        mse.loess1 <- mean((Y - Yhat.loess1)^2)

        cat(sprintf("Local linear regression:\n"))
        cat(sprintf("    span       = %f\n", span1))
        cat(sprintf("    CV score   = %f\n", score.loess1))
        cat(sprintf("    MSE        = %f\n", mse.loess1))


        span2 <- loocv.search(
                X, Y,
                loocv.loess2,
                cv.span2.min,
                cv.span2.max,
                k = cv.span2.k
        )

        score.loess2 <- loocv.loess2(X, Y, span2)
        lo2 <- loess(y ~ x, data.frame(x = X, y = Y), span = span2, degree = 2)
        y.loess2 <- predict(lo2, data.frame(x = grid))
        Yhat.loess2 <- predict(lo2, data.frame(x = X))
        mse.loess2 <- mean((Y - Yhat.loess2)^2)

        cat(sprintf("Local quadratic regression:\n"))
        cat(sprintf("    span       = %f\n", span2))
        cat(sprintf("    CV score   = %f\n", score.loess2))
        cat(sprintf("    MSE        = %f\n", mse.loess2))
        cat("\n")


        df <- data.frame(
                x = X,
                y = Y
        )

        df.regress <- data.frame(
                x = grid,
                y.depth = y.depth,
                y.kernel = y.kernel,
                y.loess1 = y.loess1,
                y.loess2 = y.loess2
        )


        # plot(X, Y)
        # lines(grid, y.depth,  col = "red")
        # lines(grid, y.kernel, col = "green")
        # lines(grid, y.loess1, col = "blue")
        # lines(grid, y.loess2, col = "cyan")
        # title(title)

        ggplot(df.regress) +
                geom_point(data = df, aes(x = x, y = y), color = purple, shape = 4) +
                geom_line(aes(x = x, y = y.depth),  color = "black") +
                geom_line(aes(x = x, y = y.kernel), color = "red",   linetype = "dashed") +
                geom_line(aes(x = x, y = y.loess1), color = "green", linetype = "dashed") +
                geom_line(aes(x = x, y = y.loess2), color = "blue",  linetype = "dashed") +
                theme.data +
                labs(
                        x = xlab,
                        y = ylab,
                        title = title
                )

}


# p.cars <- do.regression(
#         cars$speed,
#         cars$dist,
#         "Cars",
#         xlab = "Speed",
#         ylab = "Distance",
#         cv.span2.min = 0.2
# )

p.mcycle <- do.regression(
        mcycle$times,
        mcycle$accel,
        "Motorcycle accidents",
        xlab = "Time after impact",
        ylab = "Acceleration"
)

p.prestige <- do.regression(
        carData::Prestige$income,
        carData::Prestige$prestige,
        "Pineo-Porter prestige scores",
        xlab = "Average income",
        ylab = "Prestige",
        cv.beta.min = 0.05,
        cv.beta.max = 0.30,
        cv.bw.min = 1000,
        cv.bw.max = 3000,
        cv.bw.k   = 100
)

# X <- c(rnorm(50, mean = 0, sd = 2), rnorm(50, mean = 10, sd = 1))
# Y <- sin(pi * X / 10) + rnorm(length(X), sd = 0.2)
# do.regression(
#         X,
#         Y,
#         "Sinusoidal (sin(pi x / 10) vs x)"
# )
# x <- seq(min(X), max(X), length.out = 100)
# y <- sin(pi * x / 10)
# lines(x, y, col = "black")


cairo_pdf("../images/localregression_univariate.pdf", onefile = TRUE, width = 8, height = 4)

plot_grid(p.mcycle, p.prestige)

dev.off()



# Dataset - Cars (distance vs speed)
# Local depth based regression:
#     beta       = 0.106939
#     CV score   = 248.263062
#     MSE        = 159.495555
# Kernel based regression:
#     bandwidth  = 4.400000
#     CV score   = 248.691254
#     MSE        = 200.619709
# Local linear regression:
#     span       = 0.320408
#     CV score   = 235.847627
#     MSE        = 175.153019
# Local quadratic regression:
#     span       = 0.395918
#     CV score   = 239.262779
#     MSE        = 169.564191

# There were 50 or more warnings (use warnings() to see the first 50)

# Dataset - Motorcycle (acceleration vs time)
# Local depth based regression:
#     beta       = 0.114694
#     CV score   = 584.884360
#     MSE        = 454.782872
# Kernel based regression:
#     bandwidth  = 2.500000
#     CV score   = 595.958565
#     MSE        = 446.308651
# Local linear regression:
#     span       = 0.210204
#     CV score   = 554.758697
#     MSE        = 472.380389
# Local quadratic regression:
#     span       = 0.338776
#     CV score   = 540.450820
#     MSE        = 464.511804


# Dataset - Prestige (prestige vs income)
# Local depth based regression:
#     beta       = 0.192857
#     CV score   = 139.710280
#     MSE        = 113.255943
# Kernel based regression:
#     bandwidth  = 2454.545455
#     CV score   = 133.699582
#     MSE        = 690.007941
# Local linear regression:
#     span       = 0.981633
#     CV score   = 123.228910
#     MSE        = 119.864050
# Local quadratic regression:
#     span       = 1.000000
#     CV score   = 126.237009
#     MSE        = 118.369160


# Dataset - Sinusoidal (sin(pi x / 10) vs x)
# Local depth based regression:
#     beta       = 0.072041
#     CV score   = 0.044778
#     MSE        = 0.025564
# Kernel based regression:
#     bandwidth  = 1.200000
#     CV score   = 0.041262
#     MSE        = 0.259527
# Local linear regression:
#     span       = 0.504082
#     CV score   = 0.036211
#     MSE        = 0.032993
# Local quadratic regression:
#     span       = 1.000000
#     CV score   = 0.034056
#     MSE        = 0.031603

