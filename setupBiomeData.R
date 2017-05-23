## Functions for prepping the microbiome data

#' load_mbiome
#'
#' loads microbiom data given to us by Devin
#' @return list containing abundance-data, all-data,
#'         and meta-data.
load_mbiome <- function(){
    require(biom)
    require(plyr)

    homedir <- "~/Documents/PurdomGroup/data/"
    data <- read_biom(paste0(homedir,"otu.biom"))
    sample_data <- read.csv(paste0(homedir,"Till_metadata_nodepth.txt"),
                            header = TRUE,
                            stringsAsFactors = FALSE)

    ## tmp <- biom_data(data)
    ## the S4 method throws an error...
    tmp <- as.list(data)

    tmp_data <- ldply(tmp$data, function(row){
        as.data.frame(matrix(as.numeric(row), 1, 101))
    })

    tmp_data <- t(tmp_data)

    obs_meta_dat <- observation_metadata(data)
    colnames(tmp_data) <- rownames(data)
    rownames(tmp_data) <- colnames(data)

    biodata <- as.data.frame(tmp_data)
    biodata$SampleID <- rownames(tmp_data)
    
    return(list(abundance = biodata,
                sample_data = sample_data,
                obs_data= obs_meta_dat))
}

#' rarefy
#'
#' downsample the data so that each sample has the same
#' library size as the smallest library in the original data
#' @param data - data.frame of abundances
#' @param by - character "row" or "column" indicating whether
#'             the samples are rows or columns
#' @return a rarefied version of the dataset
rarefy <- function(data, by = "row", seed = 0){
    set.seed(seed)
    if(!(by %in% c("row", "column"))){
        stop("by must be 'row' or 'column")
    }
    if(by == "column"){
        tmp_data <- t(data)
    } else{
        tmp_data <- data
    }
    level <- min(rowSums(tmp_data))
    cNames <- colnames(data)
    downSamples <- ldply(1:nrow(tmp_data), function(i){
        row <- tmp_data[i,]
        vals <- rep(cNames, times = row)
        dwn_sample <- sample(vals,
                             level,
                             replace = TRUE) ## hmmm... replace?
        dwn_sample <- table(dwn_sample)
        missing <- !(cNames %in% names(dwn_sample))
        missing_names <- (cNames)[missing]
        missing_values <- rep(0, length(missing_names))
        names(missing_values) <- missing_names
        dwn_sample <- c(dwn_sample, missing_values)
        dwn_sample <- dwn_sample[order(as.integer(names(dwn_sample)))]
        dwn_sample
    })
    downSamples
}


#' chi_dist
#'
#' compute the chi-squared distance between columns of a matrix
#' @param x
chi_dist <- function(x){
    x_dim <- dim(x)
    chi_mat <- matrix(rep(0, x_dim[1]^2), x_dim[1], x_dim[1])
    for(i in 1:x_dim[1]){
        for(j in i:x_dim[1]){
            v1 <- x[i,]
            v2 <- x[j,]
            zeros <- v1 == 0 & v2 == 0
            v1 <- v1[!zeros]
            v2 <- v2[!zeros]
            dst <- .5 * sum((v1 - v2)^2 / (v1 + v2))
            chi_mat[i,j] <- dst
        }
    }
    chi_mat <- chi_mat + t(chi_mat)
}


#' run_ftests
#'
#' runs f tests using a variety of distance metrics and
#' normaliztion methods.
#' @param abundace - combination of sample data and abundance counts
run_ftests <- function(abundance){
    require(ecodist)
    require(vegan)
    require(edgeR)
    sample_data <- abundance[, 1:6]
    abundance <- abundance[, 7:length(colnames(abundance))]
    
    ## create the list of tests to run
    n_methods <- c("TMM", "RLE", "RAR", "REL", "NON")
    dist_methods <- list(dist, bcdist, chi_dist) #entropy?
    names(dist_methods) <- c("EUC", "BCD", "CHI") # "SHA"?
    tests <- expand.grid(method = n_methods, metric = names(dist_methods))

    dsets <- llply(n_methods, function(m){
        ## normalize data using TMM, RLE, and rarification:
        if(m == "RAR"){
            ret <- rarefy(abundance, by = "row", seed = 122623)
        } else if(m == "TMM" | m == "RLE") {
            d <- DGEList(counts = t(abundance))
            d <- calcNormFactors(d, method = m)
            ret <- abundance * d$samples$norm.factors
        } else  if (m == "REL"){
            row_sums <- rowSums(abundance)
            ret <- abundance / row_sums
        } else {
            ret <- abundance
        }
        as.matrix(ret)
    })

    names(dsets) <- n_methods

    groups <- sample_data %>%
        dplyr::select(Till, Covercrop, SampleType) %>%
            dplyr::distinct()

    ## conduct F test:
    scores <- mdply(tests, function(method, metric){
        fnc <- dist_methods[[metric]]
        dmat <- as.matrix(fnc(dsets[[method]]))
        overall_means <- ldply(unique(sample_data$SampleType), function(type){
            i <- which(sample_data$SampleType == type)
            data.frame(OVR_MN = mean(as.dist(as.matrix(dmat)[i,i])),
                       nsamps = length(i))
        })
        overall_means$SampleType <- unique(sample_data$SampleType)
        dmat <- as.matrix(dmat)
        within_group_vars <- mdply(groups,
                                   function(Till, Covercrop, SampleType){
            i <- which(sample_data$Till == Till &
                       sample_data$Covercrop == Covercrop &
                       sample_data$SampleType == SampleType)
            group <- as.dist(dmat[i,i])
            g_mean <- mean(group)
            sq_err <- sum((group - g_mean)^2)
            data.frame(T = Till, C = Covercrop, ST = SampleType,
                       WGV = sq_err / (dim(dmat)[1] - dim(groups)[1]),
                       GMEAN = g_mean, SIZE = length(i),
                       stringsAsFactors = FALSE)
        }, .inform = TRUE)

        vars <- left_join(within_group_vars, overall_means, by = "SampleType")
        
        vars %>%
            ## 3 = #treatments - 1
            dplyr::mutate(VAR_FROM_GMEAN = (SIZE * (GMEAN - OVR_MN)^2) / 3) %>%
            dplyr::group_by(SampleType) %>%
            dplyr::summarise(F = sum(VAR_FROM_GMEAN) / sum(WGV)) %>%
            dplyr::ungroup()
    })
    
    scores
}
