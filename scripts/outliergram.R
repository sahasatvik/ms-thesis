#!/usr/bin/env Rscript

library(ggplot2)
library(ggrepel)
library(extrafont)

MEI <- function(X, Xref) apply(X, 1, function(x) mean(apply(Xref, 1, function(y) mean(x <= y))))
MBD <- function(X, Xref) rowMeans(apply(combn(1:nrow(Xref), 2), 2, function(idx) {
                U <- apply(Xref[idx, ], 2, max)
                L <- apply(Xref[idx, ], 2, min)
                apply(X, 1, function(x) mean((L <= x) & (x <= U)))
}))

outliergram <- function(X, mbd.cut = 0.2) {
        n <- nrow(X)
        mei <- MEI(X, X)
        mbd <- MBD(X, X)
        # mbd <- roahd::MBD(X)

        a0 <- a2 <- -2 / (n * (n - 1))
        a1 <- 2 * (n + 1) / (n - 1)

        parabola <- function(x) (a0 + a1 * x + a2 * n^2 * x^2)

        x <- seq(0, 1, length.out = 100)
        y <- parabola(x)

        d <- parabola(mei) - mbd
        q <- quantile(d, probs = c(0.25, 0.75))
        d.cut <- q[2] + 1.5 * diff(q)

        y.cut <- y - d.cut
        outlier <- factor(
                (d >= d.cut) +
                2 * (mbd < mbd.cut),
                labels = c("F", "d", "mbd")
        )

        df.og <- data.frame(
                id = 1:n,
                x = mei,
                y = mbd,
                d = d,
                outlier = outlier
        )

        df.lines <- data.frame(
                x = x,
                y = y,
                y.cut = y.cut
        )

        orange <- "#e99352"
        purple <- "purple"
        red <- "red"

        theme.og <- list(
                scale_color_manual(values = c(orange, purple, red)),
                scale_shape_manual(values = c(3, 4, 4)),
                xlim(0, 1),
                labs(
                        x = "MEI",
                        y = "MBD",
                        title = "Outliergram"
                ),
                theme(
                        legend.position = "none",
                        plot.background = element_rect(fill = "transparent", color = NA),
                        panel.background = element_rect(fill = "#f7f7f7", colour = NA),
                        text = element_text(family = "Fira Sans"),
                )
        )

        p.og <- ggplot(df.og, aes(x = x, y = y)) +
                geom_ribbon(data = df.lines, aes(ymin = y.cut, ymax = y), fill = orange, alpha = 0.2) +
                geom_line(data = df.lines, color = orange) +
                geom_point(aes(shape = outlier, color = outlier)) +
                geom_text_repel(aes(label = ifelse(outlier != "F", as.character(id), ""), color = outlier), direction = "both") +
                theme.og

        return(list(plot = p.og, df = df.og))
}
