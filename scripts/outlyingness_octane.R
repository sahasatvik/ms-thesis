#!/usr/bin/env Rscript

library(rrcov)
library(reshape2)
library(cowplot)
library(tidyr)
library(dplyr)
library(ggplot2)

source("outliergram.R")
source("centrality_stability.R")
source("MO_VO.R")

data(octane)
wave <- seq(1102, 1552, by = 2)
octane <- octane[, -1]
colnames(octane) <- wave

octane.outliers <- c(25, 26, 36, 37, 38, 39)


og <- outliergram(octane, mbd.cut = 0.1, ids = octane.outliers, x.min = 0.1, x.max = 0.9)

cs <- centrality_stability(
        octane,
        outlyingness = function(y) abs((y - median(y)) / mad(y))
)
octane.cs <- cs$df
octane.cs <- mutate(octane.cs, outlier = factor(
        as.numeric(stability > 1.0),
        levels = 0:1,
        labels = c("F", "s")
))

octane.mv <- MO_VO(
        octane,
        outlyingness = function(y, y.ref) ((y - median(y.ref)) / mad(y.ref))
)
octane.mv <- mutate(octane.mv, outlier = factor(
        as.numeric(VO > 10),
        levels = 0:1,
        labels = c("F", "s")
))


orange <- "#e99352"
purple <- "purple"
red <- "red"

theme.data <- list(
        scale_color_manual(values = c(orange, purple, red)),
        labs(
                x = "Wavelength",
                y = "Response",
                title = "Data"
        ),
        theme(
                legend.position = "none",
                plot.background = element_rect(fill = "transparent", color = NA),
                panel.background = element_rect(fill = "#f7f7f7", colour = NA),
                text = element_text(family = "Fira Sans"),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                # axis.ticks.x = element_blank(),
                # axis.ticks.y = element_blank(),
                # axis.text.x = element_blank(),
                # axis.text.y = element_blank(),
        )
)

theme.cs <- list(
        scale_color_manual(values = c(orange, purple)),
        scale_shape_manual(values = c(3, 4)),
        labs(
                x = "Centrality deviation",
                y = "Stability deviation",
                title = "Centrality-Stability diagram"
        ),
        theme(
                legend.position = "none",
                plot.background = element_rect(fill = "transparent", color = NA),
                panel.background = element_rect(fill = "#f7f7f7", colour = NA),
                text = element_text(family = "Fira Sans"),
        )
)

theme.mv <- list(
        scale_color_manual(values = c(orange, purple)),
        scale_shape_manual(values = c(3, 4, 4, 4)),
        labs(
                x = "|MO|",
                y = "VO",
                title = "MO-VO diagram"
        ),
        theme(
                legend.position = "none",
                plot.background = element_rect(fill = "transparent", color = NA),
                panel.background = element_rect(fill = "#f7f7f7", colour = NA),
                text = element_text(family = "Fira Sans"),
        )
)


octane.df <- as.data.frame(octane)
octane.df$id <- 1:nrow(octane)
octane.df <- melt(octane.df, id.var = "id")
octane.df$variable <- wave[octane.df$variable]
# octane.df <- merge(octane.df, og$df[, c("id", "outlier")])

p.data <- ggplot(octane.df, aes(x = variable, y = value, group = as.factor(id))) +
        geom_line(linewidth = 0.5, color = "grey") +
        geom_line(data = filter(octane.df, id %in% octane.outliers), linewidth = 0.2, color = purple, alpha = 0.6) +
        # geom_line(data = filter(octane.df, outlier == "mbd"), linewidth = 0.2, color = red, alpha = 0.8) +
        scale_x_continuous(expand = expansion(mult = 0.01)) +
        theme.data


p.cs <- ggplot(octane.cs, aes(x = centrality, y = stability, shape = outlier, color = outlier)) +
        geom_point() +
        geom_text_repel(aes(label = ifelse(outlier != "F", as.character(id), "")), direction = "both") +
        scale_y_log10() +
        theme.cs


p.mv <- ggplot(octane.mv, aes(x = abs(MO), y = VO, shape = outlier, color = outlier)) +
        geom_point() +
        geom_text_repel(aes(label = ifelse(outlier != "F", as.character(id), "")), direction = "both") +
        scale_y_log10() +
        theme.mv



cairo_pdf("../images/outlyingness_octane.pdf", width = 8, height = 8)

plot_grid(p.data, og$plot, p.cs, p.mv, ncol = 2)

dev.off()


cairo_pdf("../images/octane.pdf", width = 8, height = 4)

ggplot(octane.df, aes(x = variable, y = value, group = as.factor(id))) +
        geom_line(linewidth = 0.5, color = "grey") +
        geom_line(data = filter(octane.df, id %in% octane.outliers), linewidth = 0.2, color = purple, alpha = 0.6) +
        scale_x_continuous(expand = expansion(mult = 0.01)) +
        scale_color_manual(values = c(orange, purple, red)) +
        labs(
                x = "Wavelength",
                y = "Response",
                title = "NIR spectra of gasoline samples"
        ) +
        theme(
                legend.position = "none",
                plot.background = element_rect(fill = "transparent", color = NA),
                panel.background = element_rect(fill = "#f7f7f7", colour = NA),
                text = element_text(family = "Fira Sans"),
        )

dev.off()


O.df <- as.data.frame(cs$O)
O.df[, "id"] <- as.factor(1:nrow(octane))
O.df <- pivot_longer(O.df, -id)
O.df$name <- as.double(O.df$name)

theme.hm <- list(
        labs(
                x = "Wavelength",
                y = NULL,
                title = "Outlyingness heatmap"
        ),
        theme(
                legend.position = "none",
                plot.background = element_rect(fill = "transparent", color = NA),
                panel.background = element_rect(fill = "#f7f7f7", colour = NA),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                text = element_text(family = "Fira Sans"),
        ),
        # scale_fill_gradient2(high = orange),
        scale_fill_gradientn(colors = c("white", orange, purple, "black")),
        # scale_fill_viridis_c(option = "inferno", direction = -1),
        scale_x_continuous(expand = c(0, 0)),
        scale_y_discrete(guide = guide_axis(n.dodge = 2), expand = c(0, 0))
)

p.hm <- ggplot(O.df, aes(x = name, y = id)) +
        geom_tile(aes(fill = value)) +
        theme.hm

p.hm2 <- ggplot(mutate(O.df, value = ifelse(id %in% octane.outliers, 0, value)), aes(x = name, y = id)) +
        geom_tile(aes(fill = value)) +
        theme.hm +
        ggtitle("(Zeroed curves #25, 26, 36-39)")


cairo_pdf("../images/outlyingness_heatmap_octane.pdf", width = 8, height = 4)

plot_grid(p.hm, p.hm2)

dev.off()
