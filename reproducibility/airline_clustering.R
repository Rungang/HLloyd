#clear()
source('HOLloyd.R')
## Input the 'flight_route.RData'
load('flight_route.RData')
r.candidate <- c(3,4,5) # when r = 6, there will be a singleton cluster
rand.seed = 2000
tensor_sample <- as.tensor(air_tensor)

# select rank
set.seed(rand.seed)
select_rank <- function(input.data, r.candidate, linename, portname){
  exp_var <- vector()
  for (r1 in r.candidate){
    for (r2 in r.candidate){
      r_use <- c(r1, r2, r2)
      z.HOSC <- HO.SC(input.data,r_use)
      z.Lloyd.HOSC = HO.Lloyd(input.data, z.HOSC)
      W.hat.Lloyd.HOSC = toweightmatrix(z.Lloyd.HOSC)
      Z.hat.Lloyd.HOSC = tomembershipmatrix(z.Lloyd.HOSC)
      X.hat.Lloyd.HOSC = ttl(ttl(input.data, W.hat.Lloyd.HOSC ,1:3), Z.hat.Lloyd.HOSC, 1:3)
      tensor_mean = ttl(input.data, W.hat.Lloyd.HOSC ,1:3)
      SST <- sum((input.data@data - mean(input.data@data))^2)
      SSR_lloyd <- sum((input.data@data  - X.hat.Lloyd.HOSC@data)^2)
      var_lloyd <- 1 - SSR_lloyd/SST
      d.prod <- length(linename) * (length(portname))^2
      BIC <- log(SSR_lloyd) + log(d.prod)/d.prod * ( r1 * r2^2 + length(linename) * log(r1) + 2*length(portname)*log(r2) )
      exp_var <- rbind(exp_var, c(r_use,var_lloyd, BIC))
    }
  }
  return(exp_var)
}
r.var <- select_rank(as.tensor(air_tensor), r.candidate, pick_linename, pick_portname)
r.var
r_select <- r.var[which.min(r.var[,5] ),1:3] # r_select = c(5,5,5), when it is (5,6,6), there will be a single cluster

# Clustering result based on Spectral initialization + HLloyd
z.HOSC <- HO.SC(tensor_sample,r_select)
z.Lloyd.HOSC = HO.Lloyd(tensor_sample, z.HOSC)

airline_clu_result <- vector()

for (i in c(1:r_select[1])){
  cluster <- select_airl_info[select_airl_info$V4 %in% pick_linename[z.Lloyd.HOSC[[1]]==i],]
  cluster$clu_num <- i
  print(cluster)
  airline_clu_result <- rbind(airline_clu_result, cluster)
}

airport_clu_result <- vector()
for (i in c(1:r_select[2])){
  cluster <- select_airport_info[select_airport_info$V5 %in%pick_portname[z.Lloyd.HOSC[[2]] == i],  ]
  cluster$clu_num <- i
  print(cluster)
  airport_clu_result <- rbind(airport_clu_result, cluster)
}


############### Clustering directly based on HOSVD or CP result
set.seed(rand.seed)
z.HOSVD <- SC(tensor_sample,r_select)

airline_clu_result_hosvd <- vector()
for (i in c(1:r_select[1])){
  cluster <- select_airl_info[select_airl_info$V4 %in% pick_linename[z.HOSVD[[1]]==i],]
  cluster$clu_num <- i
  print(cluster)
  airline_clu_result_hosvd <- rbind(airline_clu_result_hosvd, cluster)
}


############## Cluster mean estimation based on proposed method. 
W.hat.Lloyd.HOSC = toweightmatrix(z.Lloyd.HOSC)
Z.hat.Lloyd.HOSC = tomembershipmatrix(z.Lloyd.HOSC)
X.hat.Lloyd.HOSC = ttl(ttl(tensor_sample, W.hat.Lloyd.HOSC ,1:3), Z.hat.Lloyd.HOSC, 1:3)
tensor_mean = ttl(tensor_sample, W.hat.Lloyd.HOSC ,1:3)


heat.matrix = round(tensor_mean[3,,]@data,3)  # replace your matrix here. '3' is mixture cluster and '1' is US cluster

# Match the second-mode and third-mode clustering order.
heat.matrix <- heat.matrix[,c(3,4,5,1,2)]

colnames(heat.matrix) = c("Cluster 1","Cluster 2","Cluster 3","Cluster 4", "Cluster 5")
rownames(heat.matrix) = c("Cluster 1","Cluster 2","Cluster 3","Cluster 4", "Cluster 5")

colnames(heat.matrix) = c("USA+Canada","SA","EU","China", "USA+UK")
rownames(heat.matrix) = c("USA+Canada","SA","EU","China", "USA+UK")

library(reshape2)
library(tidyverse)
melted_cormat <- melt(heat.matrix)
p1 = ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2( high = "blue", mid = "white",
                        midpoint = 0, limit = c(0,0.22), space = "Lab",
                        name="p-values") + ggtitle("Mixture Airline Cluster Estimated Block Mean")+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),legend.position = "none",
        text = element_text(size=20))+
  coord_fixed() +
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 6) +
  xlab("") + ylab("")
p1







