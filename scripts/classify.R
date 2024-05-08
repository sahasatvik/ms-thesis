#!/usr/bin/env Rscript

library(ddalpha)
library(ggplot2)
library(cowplot)
library(tibble)
library(dplyr)
library(akima)
library(extrafont)

orange <- "#e99352"
purple <- "purple"


theme.data <- list(
        scale_color_manual(values = c(orange, purple)),
        scale_fill_gradient(low = orange, high = purple, na.value = "transparent"),
        scale_shape_manual(values = c(3, 4)),
        labs(
                x = "x",
                y = "y",
                title = "Data"
        ),
        theme(
                legend.position = "none",
                plot.background = element_rect(fill = "transparent", color = NA),
                panel.background = element_rect(fill = "#f7f7f7", colour = NA),
                text = element_text(family = "Fira Sans"),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.ticks.x = element_blank(),
                axis.ticks.y = element_blank(),
                axis.text.x = element_blank(),
                axis.text.y = element_blank(),
        ),
        scale_x_continuous(expand = expansion(mult = 0.0)),
        scale_y_continuous(expand = expansion(mult = 0.0))
)

theme.dd <- list(
        # geom_abline(aes(intercept = 0, slope = 1)),
        scale_color_manual(values = c(orange, purple)),
        scale_fill_gradient(low = orange, high = purple, na.value = "transparent"),
        scale_shape_manual(values = c(3, 4)),
        theme(
                legend.position = "none",
                plot.background = element_rect(fill = "transparent", color = NA),
                panel.background = element_rect(fill = "#f7f7f7", colour = NA),
                text = element_text(family = "Fira Sans"),
                axis.title.x = element_text(colour = orange),
                axis.title.y = element_text(colour = purple),
        )
)

show_classify <- function(
        A, B,
        depth = depth.spatial,
        depth.name = "spatial",
        train.ratio = 0.5,
        grid.n = 200,
        expand = 0.05,
        accuracy.only = FALSE
) {
        dfA <- tibble(group = "A", k = 1, x = A[, 1], y = A[, 2])
        dfB <- tibble(group = "B", k = 2, x = B[, 1], y = B[, 2])
        df <- rbind(dfA, dfB)
        df$group <- factor(df$group, levels = c("A", "B"))

        n <- nrow(df)
        idx.train <- sample(1:n, train.ratio * n)
        df.train <- df[ idx.train, ] %>% arrange(k)
        df.test  <- df[-idx.train, ]

        train.A <- df.train[df.train$group == "A", c("x", "y")]
        train.B <- df.train[df.train$group == "B", c("x", "y")]

        df.test <- df.test %>%
                mutate(depth.A = depth(cbind(x, y), train.A)) %>%
                mutate(depth.B = depth(cbind(x, y), train.B))

        classifier <- ddalpha.train(as.matrix(df.train[, c("x", "y", "k")]), depth = depth.name, separator = "polynomial")
        df.test <- mutate(
                        df.test,
                        k.pred = unlist(ddalpha.classify(classifier, cbind(x, y)))
                )

        if (accuracy.only) { return(mean(df.test$k == df.test$k.pred)) }

        delta.X <- diff(range(df$x)) * expand
        delta.Y <- diff(range(df$y)) * expand
        grid.X <- seq(min(df$x) - delta.X, max(df$x) + delta.X, length.out = grid.n)
        grid.Y <- seq(min(df$y) - delta.Y, max(df$y) + delta.Y, length.out = grid.n)
        grid <- expand.grid(X = grid.X, Y = grid.Y) %>%
                mutate(
                        depth.A = depth(cbind(X, Y), train.A),
                        depth.B = depth(cbind(X, Y), train.B),
                        k.pred = unlist(ddalpha.classify(classifier, cbind(X, Y)))
                )

        # grid.interp <- with(grid, interp(
        #         x = depth.A,
        #         y = depth.B,
        #         z = k.pred,
        #         linear = TRUE,
        #         extrap = FALSE,
        #         xo = seq(min(grid$depth.A), max(grid$depth.A), length.out = grid.n),
        #         yo = seq(min(grid$depth.B), max(grid$depth.B), length.out = grid.n)
        # ))

        # grid.interp <- as.data.frame(interp2xyz(grid.interp))
        # grid.interp <- mutate(
        #         grid.interp,
        #         k.pred = case_when(z <= 0.8 ~ NA, z <= 1.5 ~ 1, z > 2.2 ~ NA, z > 1.5 ~ 2, .default = NA)
        # )

        dd.x <- seq(0, 1, length.out = grid.n)
        dd.degree <- classifier$classifiers[[1]]$degree
        dd.polynomial <- classifier$classifiers[[1]]$polynomial
        dd.y <- rowSums(sapply(1:dd.degree, function(i) dd.x^i * dd.polynomial[i]))
        dd <- data.frame(x = dd.x, y = dd.y)
        print(classifier)

        df.plot <- ggplot(
                        df.train[sample(nrow(df.train)), ]
                ) +
                geom_raster(data = grid, aes(x = X, y = Y, fill = k.pred), alpha = 0.10) +
                geom_contour(data = grid, aes(x = X, y = Y, z = k.pred), breaks = c(1.5), color = "black") +
                geom_point(aes(x = x, y = y, colour = group, shape = group), size = 1.5, alpha = 1) +
                theme.data +
                ggtitle("Training Data")

        if (classifier$patterns[[1]]$name == 1) {
                df.dd <- ggplot(
                                df.test[sample(nrow(df.test)), ]
                        ) +
                        geom_ribbon(data = dd, aes(x = x, ymin = y, ymax = 1), fill = purple, alpha = 0.10) +
                        geom_ribbon(data = dd, aes(x = x, ymax = y, ymin = 0), fill = orange, alpha = 0.10) +
                        geom_line(data = dd, aes(x = x, y = y)) +
                        geom_point(aes(x = depth.A, y = depth.B, colour = group, shape = group), size = 1.5, alpha = 1) +
                        theme.dd +
                        coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
                        labs(
                                x = "Depth (Orange)",
                                y = "Depth (Purple)",
                                title = "Test Data (DD Plot)"
                        )
        } else {
                df.dd <- ggplot(
                                df.test[sample(nrow(df.test)), ]
                        ) +
                        geom_ribbon(data = dd, aes(x = x, ymin = y, ymax = 1), fill = orange, alpha = 0.10) +
                        geom_ribbon(data = dd, aes(x = x, ymax = y, ymin = 0), fill = purple, alpha = 0.10) +
                        geom_line(data = dd, aes(x = x, y = y)) +
                        geom_point(aes(x = depth.B, y = depth.A, colour = group, shape = group), size = 1.5, alpha = 1) +
                        theme.dd +
                        coord_flip(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
                        labs(
                                x = "Depth (Purple)",
                                y = "Depth (Orange)",
                                title = "Test Data (DD Plot)"
                        )
        }

        return(list(
                p.data = df.plot,
                p.dd = df.dd,
                df = df,
                df.train = df.train,
                df.test = df.test,
                grid = grid,
                classifier = classifier
        ))
}


cairo_pdf("../images/classify.pdf", onefile = TRUE, width = 8, height = 4)

N <- 200

A <- mvrnorm(N, mu = c( 1.0, 1.0), Sigma = matrix(c(1, 0, 0, 1), 2))
B <- mvrnorm(N, mu = c(-1.0, 0.0), Sigma = matrix(c(9, -2, -2, 9), 2))
sc <- show_classify(A, B)

plot_grid(sc$p.data, sc$p.dd)
df.test <- sc$df.test
mean(df.test$k == df.test$k.pred)

dev.off()
