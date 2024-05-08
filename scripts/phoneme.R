#!/usr/bin/env Rscript

library(ddalpha)
library(fda.usc)
library(fds)

library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)

source("MO_VO.R")
source("jth_order.R")

test.maxf <- function(
        Xs, f,
        train.ratio = 0.5,
        n.projections = 0,
        sigma = "vexponential"
) {
        k <- length(Xs)
        t <- 1:ncol(Xs[[1]])
        n <- sapply(Xs, nrow)
        n.train <- floor(n * train.ratio)
        idx.train <- sapply(1:k, function(i) sample(n[i], n.train[i]), simplify = FALSE)

        if (n.projections > 0) {
                P <- rproc2fdata(n.projections, t, sigma = sigma, norm = TRUE)$data
                for (i in 1:k) {
                        Xs[[i]] <- Xs[[i]] %*% t(P)
                }
        }

        X.train <- list()
        X.test <- list()
        for (i in 1:k) {
                X <- Xs[[i]]
                X.train[[i]] <- X[ idx.train[[i]], ]
                X.test[[i]]  <- X[-idx.train[[i]], ]
        }

        confusion.matrix <- c()
        for (i in 1:k) {
                crit <- sapply(1:k, function(j) f(X.test[[i]], X.train[[j]]))
                counts <- sapply(1:k, function(j) sum(max.col(crit) == j))
                confusion.matrix <- rbind(confusion.matrix, counts)
        }
        rownames(confusion.matrix) <- 1:k

        return(list(
                accuracy = sum(diag(confusion.matrix)) / sum(confusion.matrix),
                confusion.matrix = confusion.matrix
        ))
}

depth.MO_VO <- function(X, X.ref, depth = depth.Mahalanobis) {
        Y     <- as.matrix(MO_VO(    X, X.ref)[, -1])
        Y.ref <- as.matrix(MO_VO(X.ref, X.ref)[, -1])
        depth(Y, Y.ref)
}

VO <- function(X, X.ref, depth = depth.Mahalanobis) {
        -MO_VO(X, X.ref)[, 3]
}


data(aa)
data(ao)
data(dcl)
data(iy)
data(sh)

phoneme <- list(t(aa$y), t(ao$y), t(dcl$y), t(iy$y), t(sh$y))
phoneme.names <- c("aa", "ao", "dcl", "iy", "sh")

# results <- data.frame(
#         FM   = replicate(20, test.maxf(phoneme, depth.FM)$accuracy),
#         FM.2 = replicate(20, test.maxf(phoneme, depth.FMJ)$accuracy),
#         RM   = replicate(20, test.maxf(phoneme, depth.MO_VO)$accuracy),
#         VOM  = replicate(20, test.maxf(phoneme, VO)$accuracy),
#         # P.H  = replicate(20, test.maxf(phoneme, depth.halfspace, n.projections = 10)$accuracy),
#         P.M  = replicate(20, test.maxf(phoneme, depth.Mahalanobis, n.projections = 10)$accuracy),
#         P.Sp = replicate(20, test.maxf(phoneme, depth.spatial, n.projections = 10)$accuracy)
# )

# write.table(results, "phoneme_results.txt")

results <- read.table("phoneme_results.txt")

mat.to.df <- function(X, group = "") {
        n <- nrow(X)
        df <- as.data.frame(X)
        df$id <- 1:n
        df$group <- group
        # return(df)
        df <- pivot_longer(df, -c(id, group), names_to = "t", values_to = "y")
        df$t <- as.numeric(df$t)
        df
}

phoneme.df <- bind_rows(lapply(1:5, function(i) mat.to.df(phoneme[[i]], group = phoneme.names[i])))
phoneme.df$group <- factor(phoneme.df$group, levels = phoneme.names)

phoneme.medians <- phoneme.df %>%
        group_by(group, t) %>%
        summarise(median = median(y))

theme.data <- list(
        scale_x_continuous(expand = expansion(mult = 0.01)),
        labs(
                x = "Frequency",
                y = "Log periodogram",
                title = "Phoneme periodograms",
                color = ""
        ),
        theme(
                legend.position = "bottom",
                plot.background = element_rect(fill = "transparent", color = NA),
                panel.background = element_rect(fill = "#f7f7f7", colour = NA),
                text = element_text(family = "Fira Sans"),
        )
)

id.show <- sample(max(phoneme.df$id), 100)

p.data <- ggplot(filter(phoneme.df, id %in% id.show)) +
        geom_line(aes(x = t, y = y, group = interaction(id, group), color = group), alpha = 0.05) +
        geom_line(data = phoneme.medians, aes(x = t, y = median, color = group), linewidth = 1) +
        theme.data

results.df <- pivot_longer(results, everything(), names_to = "method", values_to = "accuracy")
results.df$method <- factor(
        results.df$method,
        levels = c("FM", "FM.2", "RM", "VOM", "P.M", "P.Sp"),
        labels = c("FM", "FM²",  "RM", "VOM", "PM",  "PSp")
)


p.results <- ggplot(results.df, aes(x = accuracy, y = method)) +
        geom_boxplot(aes(fill = method)) +
        geom_jitter(size = 0.5, height = 0.1, alpha = 0.7) +
        scale_fill_manual(values=c(rep("#F8766D", 2), rep("#00BA38", 2), rep("#619CFF", 2))) +
        # xlim(NA, 1.0) +
        coord_flip() +
        labs(
                x = "Accuracy",
                y = "Method",
                title = "Classification accuracies"
        ) +
        theme(
                legend.position = "none",
                plot.background = element_rect(fill = "transparent", color = NA),
                panel.background = element_rect(fill = "#f7f7f7", colour = NA),
                text = element_text(family = "Fira Sans")
        )



cairo_pdf("../images/phoneme.pdf", width = 8, height = 4)

plot_grid(p.data, p.results)

dev.off()
