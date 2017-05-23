## Perform and visualize F-tests on the microbiome data
## using four normalization methods. Further, bootsrap
## quartile intervals around the empirical results and
## visualize.

library(edgeR)
library(ecodist)
library(entropy)
library(vegan)
library(plyr); library(dplyr)
library(ggplot2)
source("~/Documents/PurdomGroup/setupBiomeData.R")

## get the raw data
data <- load_mbiome()
abundance <- data$abundance

sID <- abundance$SampleID

## make sure we have sample data for all the data points
sample_data <- data$sample_data %>%
    dplyr::filter(SampleID %in% sID)

## remove the sampleID and data points we don't have
## sample data for
abundance <- abundance %>%
    dplyr::filter(SampleID %in% sample_data$SampleID)
sID <- abundance$SampleID
abundance <- abundance %>%
    dplyr::select(-SampleID)
abundance <- cbind(sample_data, abundance)

## calculate F statistics using the original data
scores <- run_ftests(abundance)

ggplot(scores, aes(x = method, y = F, shape = metric)) +
    geom_point() +
    facet_grid(~SampleType)

## use bootstrap to get quartile intervals of F
## statistics
set.seed(1919441)
prebootData <- abundance

## function to do the bootstrap sampling
sample_fun <- function(.data){
    count <- dim(.data)[1]
    .data %>% sample_n(count, replace = TRUE)
}

## create the bootsrapped data
bootData <- llply(1:5, function(i){
    abd <- prebootData %>%
        dplyr::group_by(Till, Covercrop, SampleType) %>%
        do(sample_fun(.)) %>%
        ungroup()
    run_ftests(abd)
})

## concatenate row-wise and calculate quantile intervals
tmp <- do.call(rbind, bootData)
limits <- tmp %>%
    dplyr::group_by(method, f) %>%
    dplyr::summarise(q2 = quantile(Fscore)[2],
                     q3 = quantile(Fscore)[4])

## join the intervals to the original F statistics
plotData <- left_join(scores, limits, by = c("method", "f"))

## Plot in a couple of ways
ggplot(plotData, aes(x = method, y = Fscore, shape = f, ymax = q3, ymin = q2)) +
    geom_point() +
    geom_errorbar()

ggplot(plotData, aes(x = method, y = Fscore, color = f)) +
    geom_point()

ggplot(tmp %>%
       dplyr::mutate(group = paste0(method, f)),
       aes(x = method, y = Fscore, group = group, color = f)) +
    geom_boxplot()

