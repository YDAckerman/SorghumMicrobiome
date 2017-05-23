## General, haphazard explorations.

library(biom)
library(plyr)
library(ggplot2)
library(reshape2)

## leave as source for now:
source("~/Documents/PurdomGroup/setupBiomeData.R")

data <- load_mbiome()

tmp_data <- data$abundance
sampdat <- data$sample_data
obsdat <- data$obs_data

i <- which(colnames(tmp_data) == "SampleID")
tmp_data <- tmp_data[, -i]

## summary Stats:
## lib size
lib_size <- rowSums(tmp_data)
hist(lib_size)


## species present
num_species <- rowSums(tmp_data > 0)
hist(num_species)

## heat map
tmp <- as.data.frame(tmp_data > 0)

tmp <- melt(data = tmp_data)
tmp$value <- tmp$value > 0

ggplot(tmp, aes(x = Var1, y = Var2)) +
    geom_tile(aes(fill = value)) +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    xlab("Sample") +
    ylab("OTU")

## the following pca and mds illustrate exactly
## why normailzation methods are needed.

## pca
pca <- prcomp(tmp_data)

ggplot(as.data.frame(pca$rotation), aes(x = PC1, y = PC2)) +
    geom_point()

## mds
d <- dist(tmp_data) # euclidean distances between the rows
fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
qplot(fit$points[,1], fit$points[,2], geom = "point")

## but what if I pca/mds after normalizing by library size:
row_sums <- rowSums(tmp_data)

ndata <- tmp_data / row_sums

npca <- prcomp(ndata)

ggplot(as.data.frame(npca$rotation), aes(x = PC1, y = PC2)) +
    geom_point()

d <- dist(ndata) # euclidean distances between the rows
fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
qplot(fit$points[,1], fit$points[,2], geom = "point")

## what happens if I up-sample:
us_data <- ldply(1:nrow(tmp_data), function(i){
    row <- tmp_data[i,]
    vals <- rep(1:7234, times = row)
    up_sample <- sample(vals,
                        107278 - sum(row),
                        replace = TRUE)
    up_sample <- table(up_sample)
    missing <- !(1:7234 %in% names(up_sample))
    missing_names <- (1:7234)[missing]
    missing_values <- rep(0, length(missing_names))
    names(missing_values) <- missing_names
    up_sample <- c(up_sample, missing_values)
    up_sample <- up_sample[order(as.integer(names(up_sample)))]
    us_row <- row + up_sample
})

d <- dist(us_data) # euclidean distances between the rows
fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
qplot(fit$points[,1], fit$points[,2], geom = "point")

## and if I down sample?

#' Function Name
#'
#' Function Description
#' @param
#' @keywords
#' @export
#' @examples
rarefy <- function(data, by = "row"){
    if(by !%in% c("row", "column")){
        stop("by must be 'row' or 'column")
    }
    if(b == "column"){
        tmp_data <- t(data)
    } else{
        tmp_data <- data
    }
    level <- min(rowSums(tmp_data))
    cNames <- colnames(data)
    nCols <- length(nCols)
    downSamples <- ldply(1:nrow(tmp_data), function(i){
        row <- tmp_data[i,]
        vals <- rep(1:nCols, times = row)
        dwn_sample <- sample(vals,
                             level,
                             replace = TRUE) ## hmmm... replace?
        dwn_sample <- table(dwn_sample)
        missing <- !(1:nCols %in% names(dwn_sample))
        missing_names <- (1:nCols)[missing]
        missing_values <- rep(0, length(missing_names))
        names(missing_values) <- missing_names
        dwn_sample <- c(dwn_sample, missing_values)
        dwn_sample <- dwn_sample[order(as.integer(names(dwn_sample)))]
        dwn_sample
    })
    downSamples
}

ds_data <- rarefy(tmp_data)

d <- dist(ds_data) # euclidean distances between the rows
fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
qplot(fit$points[,1], fit$points[,2], geom = "point")

## I think these rarefaction curves are correct now!
set.seed(35)
i <- sample(1:101, 10)
N <- seq(from = 3000, to = 38000, by = 5000)
times <-  3
rarefaction_curves <- ldply(i, function(sample_i){
    row <- tmp_data[sample_i,]
    ldply(1:times, function(iterate_j){
        ldply(seq(100, sum(row), by = 100),
                             function(rarefy_N){
                                 vals <- rep(1:7234, times = row)
                                 dwn_sample <- sample(vals,
                                                      rarefy_N,
                                                      replace = TRUE)
                                 num_unique <- length(unique(dwn_sample))
                                 data.frame(sample = sample_i,
                                            iterate = iterate_j,
                                            N = rarefy_N,
                                            num_species = num_unique)
                             })
    })
})



ggplot(rarefaction_curves,
       aes(x = N,
           y = num_species,
           group = sample,
           color = sample)) +
    geom_smooth()



## and if I do presence abscence?

pa_data <- tmp_data > 0
d <- dist(pa_data)
fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim

points <- as.data.frame(fit$points)
colnames(points) <- c("x", "y")
points <- cbind(all_data[, c("Replicate", "SampleType")], points)
ggplot(points, aes(x = x, y = y, color = SampleType)) +
    geom_point()


## and lastly (for now) a correlation heat map:
## c_data <- melt(cor(tmp_data, use = "p"))
## qplot(x=Var1, y=Var2, data=c_data,
##       fill=value, geom="tile") +
##     scale_fill_gradient2(limits=c(-1, 1))
## - terrible idea


## look at the distribution of library sizes and sample otu counts:

## lib size
lib_size <- rowSums(tmp_data)
hist(lib_size)

## 
sample_otu_counts <- rowSums(tmp_data > 0)
hist(sample_otu_counts)

## look at extent sequencing vs distinct OTU's:
qplot(x = lib_size, y = sample_otu_counts, geom = "point")

## RLE

## calculate geometric means along samples
gms <- unlist(llply(1:7234, function(i){
    otu <- tmp_data[,i]
    otu <- otu[otu != 0]
    prod(otu)^(1/101)
}))

tmp_df <- data.frame(tmp_data)
tmp_df <- ldply(1:101, function(j){
    sample <- tmp_df[j,]
    sample /  median(as.numeric(sample / gms))
})
