# data
library("tidyverse")

weights_url <- "https://raw.githubusercontent.com/nonsciencemark/flow-cytometry-pigment-integrals/main/data/weights.csv"

weights <- read.csv(url(weights_url))
