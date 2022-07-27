library(sRACIPE)
library(progress)
library("ggplot2")
library(vegan)
library(dplyr)
library("factoextra")
library('kSamples')
library('SuppDists')
library("philentropy")
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}
HH=matrix(nrow=1000);
FF=matrix(nrow=1000);
file_name = "/home/huang.liji/alluniquemode/output/rsetlist1"
temp_rdata <- loadRData(file_name)
for (j in 1:1000) {
  #########  getting robustness ###########
  d=4; N=10000;k=10;
  # 4 node gene network, 10000 models and k is what we can change
  adj_to_assay = function(rSet){
    gex <-assay(rSet)
    # normalize gex
    log_gex <- log2(gex)
    log_gex <<- log_gex
    tlog_gex <-log2(t(assay(rSet)))
    tlog_gex <- tlog_gex[rowSums(tlog_gex) != -Inf,]
    tlog_gex <<- tlog_gex
  }
  
  adj_to_assay_gex <- adj_to_assay(temp_rdata[[j]])
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
  adj_to_PCA_gex <- adj_to_PCA(temp_rdata[[j]])
  eigenv <- (adj_to_PCA_gex$sdev)^2
  #########  getting flexibility ###########
  rset  <-temp_rdata[[j]]
  params <- sracipeParams(rset)
  #get all the pramams
  v <- c('G_A', 'G_B', 'G_C', 'G_D')
  F0 <- matrix(nrow=4)
  for (i in 1:4){
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
  F <- sum(F0)*1/4
  F <- unname(F)
  HH[j] <- H;
  FF[j] <- F;
}

df <- data.frame(HH, FF)
write.csv(df,"/home/huang.liji/scripts/testresult5.csv",row.names=FALSE)

load("/home/huang.liji/alluniquemode/output/tpos.5")

### ###
### ###save combined .Rdata with elementMetaData, HH and FF ###
###setClass("CombData",slots = list(tpos="list",HH="matrix",FF="matrix"))
###final_data <- new("CombData", tpos=tpos.5,HH=HH,FF=FF)

final_data <- list(tpos.5,HH,FF)

save(final_data, file = "/home/huang.liji/scripts/combinedData5.Rdata")
### ###(using function 'loadRData' to load the saved data)