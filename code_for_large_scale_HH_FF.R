#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")

#BiocManager::install("sRACIPE")

library(sRACIPE)
library(progress)
library("ggplot2")
library(vegan)
library(dplyr)
library("factoextra")
library('kSamples')
library('SuppDists')
library("philentropy")

#compute_H_F_k.R
#offset = k*15

#first get prefix of the working directory:
#dir_prefix = '/home/huang.liji/multinodenet/'
dir_prefix = '/work/lulab/lijia/multinodenets/'

st_temp <- paste0(dir_prefix, 'list_new')
all_file = read.table(st_temp)

len = 15

HH=numeric(len);
FF=numeric(len);

simulation_results <- list()

#compute accumulated time for each circuit
time <- Sys.time()

offset = 450

for (j in (offset+1):(offset+len)) {
  st <- paste0(dir_prefix, all_file[j,1])
  demoCircuit <- readRDS(st)
  
  rSet <- sRACIPE::sracipeSimulate(circuit = demoCircuit, numModels = 10000,plots = FALSE, integrateStepSize = 0.1, simulationTime = 30)
  
  simulation_results[[j-offset]] <- rSet
  
  str_list <- strsplit(as.character(all_file[j,1]), '/')
  node_st <- substr(str_list[[1]][2], 4, 5)
  num_node <- as.integer(node_st)
  
  d=num_node; N=10000;k=10;
  adj_to_assay = function(rSet){
    gex <-assay(rSet)
    # normalize gex
    log_gex <- log2(gex)
    log_gex <<- log_gex
    tlog_gex <-log2(t(assay(rSet)))
    tlog_gex <- tlog_gex[rowSums(tlog_gex) != -Inf,]
    tlog_gex <<- tlog_gex
  }
  
  adj_to_assay_gex <- adj_to_assay(rSet)
  #get one data example
  D <- as.matrix(dist(adj_to_assay_gex, p=2));
  #get the Euclidean distance
  R <- apply(D, 2, sort);
  #get the assay after sorting the distances
  H1 <- (1/N)*sum(log10(k/(R[k, ])^k));
  H <- log10(N)+log10((pi^(d/2))/gamma(1+(d/2)))-(1/N)*log10(k)*N+(d/N)*sum(log10(R[k, ]));
  #get the robustness
  adj_to_PCA <- function(rSet){
    gex <- assay(rSet)
    # normalize gex
    log_gex <- log2(gex)
    log_gex <<- log_gex
    tlog_gex <- log2(t(assay(rSet)))
    tlog_gex <- tlog_gex[rowSums(tlog_gex) != -Inf,]
    tlog_gex <<- tlog_gex
    tlog_gex_pca <- prcomp(tlog_gex, center = TRUE, scale = TRUE)
    tlog_gex_pca <<- tlog_gex_pca
  }# racipe simulation when given an adjacency matrix
  adj_to_PCA_gex <- adj_to_PCA(rSet)
  eigenv <- (adj_to_PCA_gex$sdev)^2
  #########  getting flexibility ###########
  rset  <-rSet
  params <- sracipeParams(rset)
  #get all the pramams
  #  v <- c('G_A', 'G_B', 'G_C', 'G_D')
  n = d
  F0 <- matrix(nrow=n)
  for (i in 1:n){
    G <- c(i)
    sub <- params[ , G] <= 10
    #get parameters after knocking down
    gex.sub <- adj_to_PCA_gex$x[sub,]
    pca_gex_sub <- gex.sub
    tmp1 = 0;tmp2 = 0;tmp3 = 0;tmp4 = 0;
    tmp1 = ks.test(pca_gex_sub[,1],adj_to_PCA_gex$x[,1])
    tmp2 = ks.test(pca_gex_sub[,2],adj_to_PCA_gex$x[,2])
    tmp3 = ks.test(pca_gex_sub[,3],adj_to_PCA_gex$x[,3])
    tmp4 = ks.test(pca_gex_sub[,4],adj_to_PCA_gex$x[,4])
    F0[i] = tmp1$statistic*eigenv[1] + tmp2$statistic*eigenv[2] + tmp3$statistic*eigenv[3] + tmp4$statistic*eigenv[4]
  }
  F <- sum(F0)*1/d
  F <- unname(F)
  HH[j-offset] <- H;
  FF[j-offset] <- F;
  st0 = paste0(j,'-th iteration')
  print(st0)
  print(difftime(Sys.time(), time))
}

df <- data.frame(HH, FF)

final_data <- list(simulation_results, df)

#st_temp <- paste0(dir_prefix, '/results/XXX.rds')
#saveRDS(final_data, st_temp)
#write.csv(df, 'multi_node.results.csv', row.names=FALSE)
st_temp <- paste0(dir_prefix, '/results/result_31.rds')
saveRDS(final_data, st_temp)
write.csv(df, paste0(dir_prefix, '/results/csv_results_31.csv'), row.names=FALSE)
