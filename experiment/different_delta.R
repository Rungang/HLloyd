# In this code, we want to compare the local contraction of Hlloyd algorithm under different initialization contamination rate.

source("HOLloyd.R")
library(doParallel)
library(rTensor)
# Experiment 1-1. Diffirent initialization, different delta
set.seed(100)
r_total = 5
p_total = 50
d = 3
delta.candidate = c(0.3, 0.5, 0.7, 1,2)
perturb.candidate = seq(0,0.6,0.05)
metric = 'ARI'
server = 1
if (server) {
    exp.time = 100 # d = 3 runing time: 90 s
  cluster_num = 20
} else{
  exp.time = 10
  cluster_num = 2
}

r = rep(r_total,d)
p = rep(p_total,d)
final_result <- vector()
cl2 = makeCluster(cluster_num)
registerDoParallel(cl2)

for (delta.index in 1:length(delta.candidate)){
  for (perturb.index in 1:length(perturb.candidate)){
      combined_result <- foreach(i = 1:exp.time) %dopar% {
        library(rTensor)
        library(ssvd)
        library(MASS)
        library(gtools)
        #library(cidr)
        library(mclust)
        library(pracma)
        perturb = perturb.candidate[perturb.index]
        data = TBM.generator(p, r, delta.candidate[delta.index])
        # slightly perturb the labels.
        z.oracle.init = init.perturb(data$labels, perturb)
        z.oracle = HO.Lloyd(data$tensor, z.oracle.init)
        if (strcmp(metric,'ARI')){
            lloyd_error = 1-nARI(z.oracle, data$labels)
            init_error = 1-nARI(z.oracle.init, data$labels)
          } else{
            lloyd_error = nMCR(z.oracle, data$labels)
            init_error = nMCR(z.oracle.init, data$labels)
          }
          error = c(lloyd_error, init_error)
        return(error)
    }
    combined_result <- sapply(combined_result, function(err) err)
    mean_res = apply(combined_result,1,mean)
    sd_res = apply(combined_result,1,sd)
    per_result = c(perturb.candidate[perturb.index], delta.candidate[delta.index], mean_res,sd_res)
    print(per_result)
    final_result = rbind(final_result, per_result)
  }
}
stopCluster(cl2)

colnames(final_result) <- c('epsilon','Delta','lloyd', 'init','lloyd_sd','init_sd')
file_name <- paste("different_delta_winit","p", paste(p_total, sep="", collapse="_"), "r", paste(r_total, sep="", collapse="_"), "d", d ,"delta", paste(delta.candidate, sep="", collapse="_"), "exp_time", exp.time,'metric',metric,sep = "_")
print(file_name)
print(final_result)
write.csv(final_result, file = paste("../results/",file_name, ".csv", sep = ""))
