## perform a within/between euclidean distance comparison of the data
## using four normalization methods: TMM, RLE, RAR, NON

library(edgeR)
library(vegan)
library(dplyr)
library(plyr)
source("~/Documents/PurdomGroup/setupBiomeData.R")

data <- load_mbiome()
abundance <- data$abundance

## remove the sampleID for now:
sID <- abundance$SampleID
abundance <- abundance[, -which(colnames(abundance) == "SampleID")]

methods <- c("TMM", "RLE", "RAR", "NON")
dsets <- llply(methods, function(m){
    ## normalize data using TMM, RLE, and rarification:
    if(m == "RAR"){
        abdnc_m <- rarefy(abundance, by = "row", seed = 122623)
    } else if(m == "TMM" | m == "RLE") {
        d <- DGEList(counts = t(abundance))
        d <- calcNormFactors(d, method = m)
        abdnc_m <- abundance / d$samples$norm.factors
    } else {
        abdnc_m <- abundance
    }
    ## insert sample meta data:
    abdnc_m$SampleID <- sID
    left_join(abdnc_m, data$sample_data, by = "SampleID")
})

names(dsets) <- methods

## within- and between-group dist:

## within
dist_within <- ldply(methods, function(set){
    data <- dsets[[set]] %>%
        dplyr::filter(!is.na(Replicate) & !is.na(SampleType))
    ddply(data, .(SampleType), function(d){
        i <- which(colnames(d) %in% c("SampleID", "Batch", "Till",
                                      "Covercrop", "Replicate", "SampleType"
                                      ))
        ## calculate euclidean distances
        dmat <- dist(d[, -i])
        data.frame(method = set, max = max(dmat), min = min(dmat),
                   mean = mean(dmat), median = median(dmat))
    })
})

## between
dist_between <- ldply(methods, function(set){
    data <- dsets[[set]] %>%
        dplyr::filter(!is.na(Replicate) & !is.na(SampleType))
    dists <- ddply(data, .(Replicate), function(d){
        ldply(c("Rhizosphere", "Root", "Soil"), function(stype){
            data <- d %>%
                dplyr::filter(SampleType != stype)
            i <- which(colnames(data) %in% c("SampleID", "Batch", "Till",
                                      "Covercrop", "Replicate", "SampleType"
                                             ))
            ## calculate euclidean distances
            dmat <- dist(data[, -i])
            data.frame(method = set, excluding = stype, dist = max(dmat))
        })
    })
    dists %>%
        dplyr::group_by(method, excluding) %>%
        dplyr::summarize(max = max(dist),
                         min = min(dist),
                         mean = mean(dist),
                         median = median(dist))
})

## we want to compare within/between distances, so let's take
## the ratios and plot:
tmp <- left_join(dist_between, dist_within, by = "method")
tmp <- tmp %>%
    dplyr::filter(SampleType != excluding) %>%
    dplyr::mutate(m = mean.y / mean.x)

ggplot(tmp, aes(x = method, y = m, color = SampleType)) +
    geom_point()
