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

n_methods <- c("TMM", "RLE", "RAR", "REL", "NON")
dist_methods <- list(dist, bcdist, chi_dist) #entropy
names(dist_methods) <- c("EUC", "BCD", "CHI") # "SHA"
tests <- expand.grid(method = n_methods, f = names(dist_methods))

dsets <- llply(n_methods, function(m){
    ## normalize data using TMM, RLE, and rarification:
    if(m == "RAR"){
        ret <- rarefy(abundance, by = "row", seed = 122623)
    } else if(m == "TMM" | m == "RLE") {
        d <- DGEList(counts = t(abundance))
        d <- calcNormFactors(d, method = m)
        ret <- abundance / d$samples$norm.factors
    } else  if (m == "REL"){
        row_sums <- rowSums(abundance)
        ret <- abundance / row_sums
    } else {
        ret <- abundance
    }
    as.matrix(ret)
})

names(dsets) <- n_methods

## split into the groups we think should be most similar
groups <-  sample_data %>%
    dplyr::select(Till, Covercrop, Batch, SampleType) %>%
    dplyr::distinct()


Ks <- sample_data %>%
    dplyr::group_by(Till, Covercrop) %>%
    dplyr::summarise(K = length(unique(Batch))) %>%
    dplyr::ungroup()


scores <- mdply(tests, function(method, f){
    fnc <- dist_methods[[f]]
    dmat <- as.matrix(fnc(dsets[[method]]))
    
    ## calculate the within group variation
    WGV <- mdply(groups, function(Till, Covercrop, Batch, SampleType){
        ## this will give us the indices of this group's sample
        i <- which(sample_data$Till == Till &
                   sample_data$Covercrop == Covercrop &
                   sample_data$Batch == Batch &
                   sample_data$SampleType == SampleType)
        ## this will give us the number of groups (3, for all cases)
        j <- which(Ks$Till == Till &
                   Ks$Covercrop == Covercrop)
        ## this will give us the total number of samples
        l <- which(sample_data$Till == Till &
                   sample_data$Covercrop == Covercrop &
                   sample_data$SampleType == SampleType)
        group <- as.dist(dmat[i,i])
        g_mean <- mean(group)
        sq_err <- sum((group - g_mean)^2)
        data.frame(WGV = sq_err / (length(l) - Ks$K[j]),
                   GMEAN = g_mean, SIZE = length(i),
                   stringsAsFactors = FALSE)
    }, .inform = TRUE)

    ## calculate overall means
    overall_mean <- mdply(groups, function(Till, Covercrop, Batch, SampleType){
        i <- which(sample_data$Till == Till &
                   sample_data$Covercrop == Covercrop &
                   sample_data$SampleType == SampleType)
        data.frame(OVR_MN = mean(as.dist(dmat[i,i])))
    })
        
    WGV <- left_join(WGV, overall_mean,
                     by = c("Till", "Covercrop", "SampleType"))

    ## calculate the F statistics
    WGV %>%
        group_by(Till, Covercrop, SampleType) %>%
            dplyr::mutate(GV = SIZE * (GMEAN - OVR_MN)^2) %>%
                dplyr::summarise(F = (sum(GV) / (n() - 1)) / sum(WGV))
})

ggplot(scores, aes(x = method, y = F, color = f)) +
    geom_point() +
    facet_wrap(~SampleType)
