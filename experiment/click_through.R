source("HOLloyd.R")
library(rTensor)
library(pROC)
source("STATSVD.R")

X = read.csv("Click_through.csv", header = F)
Y = array(0,c(100,50,24,8))
for (i in 1:nrow(X)){
  Y[X[i,1],X[i,2],X[i,3]+1,X[i,4]] = 1
}

# use clustering on all the aggregated data.
set.seed(1)
Y.sum = array(0,c(100,50,24))
for (i in 1:8){
  Y.sum = Y.sum + Y[,,,i]
}
Y.tensor = as.tensor(Y.sum)/8
Y.tensor = as.tensor(Y.sum>0)
Y.clustering.initializer = HO.SC(Y.tensor, c(4,4,4))
Y.clustering = HO.Lloyd(Y.tensor, Y.clustering.initializer)
W.hat = toweightmatrix(Y.clustering)
S.est = ttl(Y.tensor, W.hat ,1:3)


#################
S.est = ttm(Y.tensor, W.hat[[3]] ,3)/8
z1 = c(which(Y.clustering[[1]]==1),which(Y.clustering[[1]]==2),
       which(Y.clustering[[1]]==3),which(Y.clustering[[1]]==4))
z2 = c(which(Y.clustering[[2]]==1),which(Y.clustering[[2]]==2),
       which(Y.clustering[[2]]==3),which(Y.clustering[[2]]==4))
S2 = S.est@data[,,1]
S2 = S2[z1,z2]
levelplot(S2)


#################
# compare HOOI and HLloyd for different dfs.
# Y.tensor = Y.tensor - mean(Y.tensor@data)
X = read.csv("tensor_data_co_clustering.csv", header = F)
Y = array(0,c(100,50,24,8))
for (i in 1:nrow(X)){
  Y[X[i,1],X[i,2],X[i,3]+1,X[i,4]] = 1
}

# use clustering on all the aggregated data.
set.seed(1)
Y.sum = array(0,c(100,50,24))
for (i in 1:8){
  Y.sum = Y.sum + Y[,,,i]
}

Y.tensor = as.tensor(Y.sum-mean(Y.sum))
Y.tensor = as.tensor(Y.sum>0)
Y.array = Y.tensor@data
Y.array[Y.array==0] = 1/8
Y.array[Y.array==1] = 1-1/8
Y.tensor = as.tensor(logit(Y.array))
set.seed(1)
r.set = seq(2,6,1)
RSE = matrix(0, 6, length(r.set))
for (r.ind in 1:length(r.set)) {
  r = r.set[r.ind]
  Y.clustering.initializer = HO.SC(Y.tensor, c(r,r,r))
  Y.clustering = HO.Lloyd(Y.tensor, Y.clustering.initializer)
  W.hat = toweightmatrix(Y.clustering)
  S.est = ttl(Y.tensor, W.hat, 1:3)
  Z.hat = tomembershipmatrix(Y.clustering)
  Y.est = ttl(S.est, Z.hat, 1:3)
  p = dim(Y.tensor)
  sigma.hat.2 = norm(k_unfold(Y.tensor - Y.est,1)@data,'F')^2/prod(p)
  df.HLloyd = r^3 + sum(p)*log2(r)
  error.HLloyd = norm(k_unfold(Y.tensor-Y.est,1)@data, 'F')^2
  RSE[1,r.ind] = (norm(k_unfold(Y.tensor,1)@data, 'F')^2 - 
    norm(k_unfold(Y.tensor-Y.est,1)@data, 'F')^2)/df.HLloyd
  RSE[1,r.ind] = prod(p)/(prod(p)-df.HLloyd)*norm(k_unfold(Y.tensor-Y.est,1)@data, 'F')^2/(norm(k_unfold(Y.tensor,1)@data, 'F')^2)
                    
  
  HOOI.est = HOSVD(Y.tensor, c(r,r,r))
  HOOI.S.est = ttm(ttm(ttm(Y.tensor, t(HOOI.est[[1]]), 1), t(HOOI.est[[2]]), 2), t(HOOI.est[[3]]),3)
  HOOI.Y.est = ttl(HOOI.S.est, HOOI.est, 1:3)
  error.HOOI = norm(k_unfold(Y.tensor-HOOI.Y.est,1)@data, 'F')^2
  sigma.HOOI.2 = norm(k_unfold(Y.tensor - HOOI.Y.est,1)@data,'F')^2/prod(p)
  df.HOOI = (r^3 + sum(dim(Y.tensor)-(r+1)/2)*r)
  RSE[2,r.ind] = (norm(k_unfold(Y.tensor,1)@data, 'F')^2 - 
                    norm(k_unfold(Y.tensor-HOOI.Y.est,1)@data, 'F')^2)/df.HOOI
  RSE[2,r.ind] = prod(p)/(prod(p)-df.HOOI)*norm(k_unfold(Y.tensor-HOOI.Y.est,1)@data, 'F')^2/(norm(k_unfold(Y.tensor,1)@data, 'F')^2)
  
  RSE[3,r.ind] = prod(p)*log((norm(k_unfold(Y.tensor-Y.est,1)@data, 'F')^2 / norm(k_unfold(Y.tensor-HOOI.Y.est,1)@data, 'F')^2)) 
  RSE[4,r.ind] = pchisq(RSE[3,r.ind], df.HOOI-df.HLloyd)
  
  RSE[5,r.ind] = (error.HLloyd-error.HOOI)/error.HOOI*(prod(p)-df.HOOI)/(df.HOOI-df.HLloyd)
  RSE[6,r.ind] = pf(RSE[5,r.ind], df.HOOI-df.HLloyd, prod(p)-df.HOOI)
}


library(ggplot2)
library(reshape2)
library(gridExtra)
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
max.intensity = round(max(as.vector(S.est@data)),3)
heat.matrix = round(S.est@data[,,1],3)
colnames(heat.matrix) = c("I1","I2","I3","I4")#,"I5","I6","I7","I8","I9","I10")
rownames(heat.matrix) = c("U1","U2","U3","U4")#,"U5","U6","U7","U8","U9","U10")
melted_cormat <- melt(heat.matrix)
p1 = ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2( high = "blue", mid = "white", 
                        midpoint = 0, limit = c(0,0.22), space = "Lab", 
                        name="p-values") + ggtitle("6pm-9pm")+
  theme_minimal()+ # minimal theme
  theme(plot.title = element_text(hjust = 0.5),legend.position = "none", 
        text = element_text(size=15))+
  coord_fixed() +
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
  xlab("") + ylab("")

heat.matrix = round(S.est@data[,,2],3)
colnames(heat.matrix) = c("I1","I2","I3","I4")#,"I5","I6","I7","I8","I9","I10")
rownames(heat.matrix) = c("U1","U2","U3","U4")#,"U5","U6","U7","U8","U9","U10")
melted_cormat <- melt(heat.matrix)
p2 = ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2( high = "blue", mid = "white", 
                        midpoint = 0, limit = c(0,0.22), space = "Lab", 
                        name="p-values") + ggtitle("12am-6am")+
  theme_minimal()+ # minimal theme
  theme(plot.title = element_text(hjust = 0.5),legend.position = "none", 
        text = element_text(size=15))+
  coord_fixed() +
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
  xlab("") + ylab("")

heat.matrix = round(S.est@data[,,3],3)
colnames(heat.matrix) = c("I1","I2","I3","I4")#,"I5","I6","I7","I8","I9","I10")
rownames(heat.matrix) = c("U1","U2","U3","U4")#,"U5","U6","U7","U8","U9","U10")
melted_cormat <- melt(heat.matrix)
p3 = ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2( high = "blue", mid = "white", 
                        midpoint = 0, limit = c(0,0.22), space = "Lab", 
                        name="p-values") + ggtitle("9pm-12am")+
  theme_minimal()+ # minimal theme
  theme(plot.title = element_text(hjust = 0.5),legend.position = "none", 
        text = element_text(size=15))+
  coord_fixed() +
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
  xlab("") + ylab("")

heat.matrix = round(S.est@data[,,4],3)
colnames(heat.matrix) = c("I1","I2","I3","I4")#,"I5","I6","I7","I8","I9","I10")
rownames(heat.matrix) = c("U1","U2","U3","U4")#,"U5","U6","U7","U8","U9","U10")
melted_cormat <- melt(heat.matrix)
p4 = ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2( high = "blue", mid = "white", 
                        midpoint = 0, limit = c(0,0.22), space = "Lab", 
                        name="p-values") + ggtitle("6am-6pm")+
  theme_minimal()+ # minimal theme
  theme(plot.title = element_text(hjust = 0.5),legend.position = "none", 
        text = element_text(size=15))+
  coord_fixed() +
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
  xlab("") + ylab("")


grid.arrange(p2, p4, p1, p3, ncol=4, nrow = 1)


# randomly pick four as training, two as validation and two as predictions.
# compare with CP low-rank, Tucker low-rank and generalized low-rank?
index.training = combn(1:8,6)
record.num = ncol(index.training) * 6
training.loss = matrix(0, record.num, 4)
#validation.loss = matrix(0, record.num, 4)
testing.loss = matrix(0, record.num, 4)
r = c(2,2,3)
HLloyd.auc = NULL
HOOI.auc = NULL



# ncol(index.training)
for (i in 1:10) {
  print(i)
  index.testing = setdiff(1:8, index.training[,i])

  Y.training = as.tensor(apply(Y[,,,index.training[,i]], 1:3, mean))
  # Y.testing = as.tensor(apply(Y[,,,index.testing], 1:3, mean))
    
  # vali.loss.min = 100000000
  # r.min = 1
  # for(k in 1:nrow(r.set)){
  #   z.initializer = HO.SC(Y.training, r.set[k,1:3])
  #   z = HO.Lloyd(Y.training, z.initializer)
  #   W.hat = toweightmatrix(z)
  #   Z.hat = tomembershipmatrix(z)
  #   X.hat = ttl(ttl(Y.training, W.hat ,1:3), Z.hat, 1:3)
  #   vali.loss = rTensor::fnorm(X.hat - Y.validation)
  #   if(vali.loss < vali.loss.min){
  #       vali.loss.min = vali.loss
  #       r.min = k
  #   }
  # }
  # r.set[r.min,4] = r.set[k,4] + 1
    z.initializer = HO.SC(Y.training, r)
    z = HO.Lloyd(Y.training, Y.clustering.initializer)
    W.hat = toweightmatrix(z)
    Z.hat = tomembershipmatrix(z)
    X.est = ttl(ttl(Y.training, W.hat ,1:3), Z.hat, 1:3)
    # training.loss[index.count,1] = rTensor::fnorm(X.est - Y.training)
    # validation.loss[index.count,1] = vali.loss.min
    test.1 = as.vector(Y[,,,index.testing[1]])
    test.2 = as.vector(Y[,,,index.testing[2]])
    predict = as.vector(X.est@data)
    HLloyd.auc = c(HLloyd.auc, auc(roc(test.1,predict)), auc(roc(test.2,predict)))
    
    
    # testing.loss[i,1] = rTensor::fnorm(X.est - Y.testing)
    
    
    # compare with HOOI.
    # vali.loss.min = 100000000
    # r.min = 1
    # for(k in 1:nrow(r.set)){
    #   my_U = HOOI(Y.training, r.set[k,], tmax = 2)
    #   P_U_1 = my_U[[1]]%*%t(my_U[[1]])
    #   P_U_2 = my_U[[2]]%*%t(my_U[[2]])
    #   P_U_3 = my_U[[3]]%*%t(my_U[[3]])
    #   X.hat = ttm(ttm(ttm(Y.training, P_U_1, 1), P_U_2, 2), P_U_3, 3)
    #   vali.loss = rTensor::fnorm(X.hat - Y.validation)
    #   if(vali.loss < vali.loss.min){
    #     r.min = k
    #     vali.loss.min = vali.loss
    #   }
      my_U = HOOI(Y.training, r, tmax = 5)
      P_U_1 = my_U[[1]]%*%t(my_U[[1]])
      P_U_2 = my_U[[2]]%*%t(my_U[[2]])
      P_U_3 = my_U[[3]]%*%t(my_U[[3]])
      X.est = ttm(ttm(ttm(Y.training, P_U_1, 1), P_U_2, 2), P_U_3, 3)
      predict = as.vector(X.est@data)
      HOOI.auc = c(HOOI.auc, auc(roc(test.1,predict)), auc(roc(test.2,predict)))
      
      # training.loss[index.count,2] = rTensor::fnorm(X.est - Y.training)
      # validation.loss[index.count,2] = vali.loss.min
      # testing.loss[index.count,2] = rTensor::fnorm(X.est - Y.testing)
      
  
    
    
}











