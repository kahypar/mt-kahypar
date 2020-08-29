#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
working_directory <- args[1]
setwd(working_directory)

install.packages("plyr", repos = "http://cran.us.r-project.org")
install.packages("ggplot2", repos = "http://cran.us.r-project.org")
install.packages("scales", repos = "http://cran.us.r-project.org")

library(plyr)
library(ggplot2)
library(scales)

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

sqrt5_trans = function() trans_new("sqrt5", function(x) x^(1/5), function(x) x^5)

running_time_plot <- function(algo_1, algo_2) {
  gmean_time_aggreg = function(df) data.frame(gmean_time = gm_mean(df$avg_time))
  result_gmean <- ddply(result, c("algorithm"), gmean_time_aggreg)
  result <- rbind(algo_1, algo_2)
  running_time = ggplot(result, aes(x=algorithm, y=avg_time, fill=algorithm)) +
    geom_point(size = 0.5, pch = 21, position = position_jitterdodge(jitter.width = 1.0)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.5) +
    geom_text(aes(x = algorithm, y = 0, label=round(gmean_time, 2), group = algorithm), result_gmean, 
              size = 2,  position = position_dodge(width=0.75)) +
    coord_trans(y = "sqrt5") + 
    theme_bw(base_size = 10) +
    labs(x="Commit", y="Running Time [s]") +
    theme(aspect.ratio =2/(1+sqrt(5)),
          legend.position = "none",
          panel.grid.major = element_line(linetype="dotted",size = 0.25, 
                                          color = "grey"),
          panel.grid.minor =element_blank(),
          axis.line = element_line(size = 0.2, color = "black"),
          axis.title.y = element_text(vjust=1.5),
          axis.title.x = element_text(vjust=1.5),
          axis.text.x = element_text(angle = 15, hjust = 1, size = 8))
  
  return(running_time)
}

############## SETUP DATA FRAMES ############## 

#setwd("/home/tobias/mt-kahypar/.github/workflows")

result_folder <- args[2]
algo_1 <- args[3]
algo_2 <- args[4]

algo_1_db <- ddply(read.csv(paste(result_folder, "/", algo_1, ".csv", sep=""), header = TRUE), c("graph", "k", "epsilon"), aggreg)
algo_1_db$algorithm <- algo_1
algo_2_db <- ddply(read.csv(paste(result_folder, "/", algo_2, ".csv", sep=""), header = TRUE), c("graph", "k", "epsilon"), aggreg)
algo_2_db$algorithm <- algo_2

png(paste(result_folder,"/","running_time.png",sep=""), width = 800, height = 600)
print(running_time_plot(algo_1_db, algo_2_db))
dev.off()
