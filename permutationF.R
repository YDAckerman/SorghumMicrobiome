## Run permutation tests on the abundance
## F statistic. I don't think I ended up using
## this...

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

set.seed(100)
numPerms <- 100
seeds <- sample(500:1000, numPerms)
permData <- llply(1:numPerms, function(i){NA})

## function to permute the data
tmpfn <- function(.data, seed){
    set.seed(seed)
    j <- sample(1:nrow(.data), nrow(.data),
                replace = FALSE)
    .data$Till <- .data$Till[j]
    .data$Covercrop <- .data$Covercrop[j]
    .data
}

for(i in 1:numPerms){
    ## permute the treatment
    abd <- abundance %>%
        dplyr::group_by(SampleType) %>%
        do(tmpfn(., seed = seeds[i])) %>%
        dplyr::ungroup()
    ## run the tests:
    permData[[i]] <- run_ftests(abd)
}

tmp <- do.call(rbind, permData)

ggplot(tmp, aes(x = method, y = F, shape = metric)) +
    geom_boxplot() +
    facet_wrap(~SampleType)


tmp <- left_join(tmp, scores %>% dplyr::rename(Fscore = F),
                 by = c("method", "metric", "SampleType"))

## incorporate actual scores:
ggplot(tmp, aes(x = method, y = F, shape = metric, color = metric)) +
    geom_boxplot() +
    geom_point(aes(x = method, y = Fscore, shape = metric), color = "red") +
    facet_wrap(~SampleType)

save(tmp, file = "~/Documents/PurdomGroup/permTest.rda")
