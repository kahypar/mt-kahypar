#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
working_directory <- args[1]
setwd(working_directory)

install.packages("plyr", repos = "http://cran.us.r-project.org")

library(plyr)

aggreg = function(df) data.frame(min_km1 = min(df$km1),
                                 avg_km1 = mean(df$km1),
                                 min_cut = min(df$cut),
                                 avg_cut = mean(df$cut),
                                 min_imbalance = min(df$imbalance),
                                 avg_imbalance = mean(df$imbalance),
                                 min_total_time = min(df$partitionTime),
                                 avg_total_time = mean(df$partitionTime),
                                 avg_time = mean(df$partitionTime))

gm_mean = function(x, na.rm=TRUE, zero.propagate = FALSE){
  if(any(x < 0, na.rm = TRUE)){
    return(NaN)
  }
  if(zero.propagate){
    if(any(x == 0, na.rm = TRUE)){
      return(0)
    }
    return(exp(mean(log(x), na.rm = na.rm)))
  } else {
    return(exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x)))
  }
}

############## SETUP DATA FRAMES ############## 

# setwd("/home/tobias/mt-kahypar/.github/workflows")

result_folder <- args[2]
algo_1 <- args[3]
algo_2 <- args[4]

algo_1_db <- ddply(read.csv(paste(result_folder, "/", algo_1, ".csv", sep=""), header = TRUE), c("graph", "k", "epsilon"), aggreg)
algo_1_db$algorithm <- algo_1
algo_2_db <- ddply(read.csv(paste(result_folder, "/", algo_2, ".csv", sep=""), header = TRUE), c("graph", "k", "epsilon"), aggreg)
algo_2_db$algorithm <- algo_2

gmean_time_algo_1 <- gm_mean(algo_1_db$avg_time)
gmean_time_algo_2 <- gm_mean(algo_2_db$avg_time)
running_time_ratio <- gmean_time_algo_2 / gmean_time_algo_1
sink("running_time_ratio.txt")
cat(running_time_ratio)
sink()

gmean_quality_algo_1 <- gm_mean(algo_1_db$avg_km1)
gmean_quality_algo_2 <- gm_mean(algo_2_db$avg_km1)
quality_ratio <- gmean_quality_algo_2 / gmean_quality_algo_1
sink("quality_ratio.txt")
cat(quality_ratio)
sink()
