#!/usr/bin/env Rscript

library(reshape2)
library(cowplot)
library(dplyr)
library(ggplot2)

load("../data/Winedata.rda")

idx  <- (Winedata$ppm < 5.62) & (Winedata$ppm > 5.37)
wave <- Winedata$ppm[idx]
wine <- Winedata$spectra[, idx]
colnames(wine) <- wave
