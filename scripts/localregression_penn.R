#!/usr/bin/env Rscript


library(ddalpha)
library(tidyr)
library(dplyr)
library(stringr)
library(cowplot)
library(ggplot2)
library(extrafont)
# library(np)

# depth <- depth.Mahalanobis
# depth <- depth.halfspace
depth <- function(x, data) depth.spatial(x, data, mah.estimate = "none")

# source("localdepth.R")

localdepth.neighbours <- function(u, X, beta, depth, symm = FALSE) {
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
        Y.hat <- c()
        idx <- c()
        for (i in 1:nrow(points)) {
                x <- points[i, ]
                neighbours <- localdepth.neighbours(x, X, beta = beta, depth = depth, symm = symm)
                depths <- depth(X[neighbours$idx, ], neighbours$X)
                # weighted.mean(Y, depths)
                Y.hat.x <- colSums(depths * Y[neighbours$idx, ]) / sum(depths)
                Y.hat <- rbind(Y.hat, Y.hat.x)
                idx <- rbind(idx, neighbours$idx)
        }
        rownames(Y.hat) <- NULL
        return(list(Y.hat = Y.hat, idx = idx))
}


# loocv <- function(X, Y, f) {
#         Yhat <- sapply(1:nrow(X), function(i) f(X[-i, ], Y[-i, ], matrix(X[i, ], 1)))
#         mean((Y - Yhat)^2, na.rm = TRUE)
# }

# loocv.apply <- function(X, Y, Z, fz) {
#         sapply(Z, function(z) {
#                 loocv(X, Y, function(X, Y, p) fz(X, Y, p, z))
#         })
# }

# loocv.localregress <- function(X, Y, beta) {
#         loocv.apply(X, Y, beta, function(X, Y, p, z)
#                 localregress(X, Y, p, beta = z, depth = depth, symm = TRUE)
#         )
# }

# loocv.search <- function(X, Y, g, z.lo, z.hi, k = 3) {
#         Z <- seq(z.lo, z.hi, length.out = k)
#         mse <- g(X, Y, Z)
#         return(Z[which.min(mse)])
# }

library(Ecdat)
data(SumHes)

t <- min(SumHes$year):max(SumHes$year)
X <- t(matrix(SumHes$gdp, 26))
Y <- t(matrix(SumHes$sr, 26))
colnames(X) <- t
colnames(Y) <- t

country <- unique(SumHes$country)
test <- c(19, 66, 107, 102, 117)
X.test <- X[test, ]
Y.test <- Y[test, ]

X.train <- X[-test, ]
Y.train <- Y[-test, ]

lr <- localregress(X.train, Y.train, X.test, 0.1, depth)
Y.hat <- lr$Y.hat
idx <- lr$idx

X.train.df <- as.data.frame(X.train)
X.train.df$id <- 1:nrow(X.train)
X.train.df <- pivot_longer(X.train.df, -id, names_to = "t", values_to = "gdp")
X.train.df$t <- as.numeric(X.train.df$t)

Y.train.df <- as.data.frame(Y.train)
Y.train.df$id <- 1:nrow(Y.train)
Y.train.df <- pivot_longer(Y.train.df, -id, names_to = "t", values_to = "sr")
Y.train.df$t <- as.numeric(Y.train.df$t)


orange <- "#e99352"
purple <- "purple"
red <- "red"
blue <- "blue"
green <- "seagreen"

theme.data <- list(
        # scale_color_manual(values = c(orange, purple, red)),
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
)

# pdf("penntable.pdf", width = 12, height = 6)
# par(mfrow = c(1, 2))

rows <- list()
for (i in 1:nrow(X.test)) {
        # matplot(t, t(X), type = "l", lty = 1, col = "grey")
        # matplot(t, t(X[idx[i, ], ]), type = "l", lty = 1, col = "green", add = TRUE)
        # lines(t, X.test[i, ], type = "l", lty = 1, lw = 2, col = "black")
        # title(country[test[i]])
        # print(dim(X[idx[i, ], ]))
        df.i <- data.frame(
                t = t,
                gdp = X.test[i, ],
                sr = Y.test[i, ],
                sr.hat = Y.hat[i, ]
        )

        p.X <- ggplot(X.train.df, aes(x = t, y = gdp)) +
                geom_line(aes(group = id), color = "grey", linewidth = 0.2) +
                geom_line(data = X.train.df %>% filter(id %in% which(idx[i, ])), aes(group = id), color = orange) +
                geom_line(data = df.i, color = "black") +
                theme.data +
                theme(
                        plot.margin = margin(5, 5, 0, 5, unit = "pt"),
                ) +
                labs(
                        x = NULL,
                        y = "GDP",
                        title = str_to_title(country[test[i]])
                )

        p.Y <- ggplot(Y.train.df, aes(x = t)) +
                geom_line(aes(group = id, y = sr), color = "grey", linewidth = 0.2) +
                geom_ribbon(
                        data = Y.train.df %>%
                                filter(id %in% which(idx[i, ])) %>%
                                group_by(t) %>%
                                summarize(
                                        ymin = min(sr),
                                        ymax = max(sr)
                                ),
                        aes(ymin = ymin, ymax = ymax),
                        fill = orange,
                        alpha = 0.20
                ) +
                geom_line(data = Y.train.df %>% filter(id %in% which(idx[i, ])), aes(group = id, y = sr), color = orange) +
                geom_line(data = df.i, aes(y = sr), color = purple) +
                geom_line(data = df.i, aes(y = sr.hat), color = "black") +
                theme.data +
                theme(
                        plot.margin = margin(5, 10, 0, 5, unit = "pt"),
                ) +
                labs(
                        x = NULL,
                        y = "Saving rate"
                )

        if (i == nrow(X.test)) {
                p.X <- p.X + xlab("Year")
                p.Y <- p.Y + xlab("Year")
        }

        # matplot(t, t(Y), type = "l", lty = 1, col = "grey")
        # matplot(t, t(Y[idx[i, ], ]), type = "l", lty = 1, col = "green", add = TRUE)
        # lines(t, Y.hat[i, ], type = "l", lty = 1, lw = 2, col = "black")
        # lines(t, Y[test[i], ], type = "l", lty = 2, lw = 2, col = "red")
        # lines(t, apply(Y[idx[i, ], ], 2, max), type = "l", lty = 2, lw = 1, col = "black")
        # lines(t, apply(Y[idx[i, ], ], 2, min), type = "l", lty = 2, lw = 1, col = "black")

        rows[[i]] <- plot_grid(p.X, p.Y)

}



cairo_pdf("../images/localregression_penntable.pdf", onefile = TRUE, width = 8, height = 11)

do.call(plot_grid, c(rows, ncol = 1))
# plot_grid(plots[[1]], plots[[2]], ncol = 2)

dev.off()
