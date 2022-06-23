# Experiment: tensor recovery
source("HOLloyd.R")
library(doParallel)
set.seed(100)
exp.time = 100 
r.candidate = c(2,3)
p.candidate = seq(40,100,10)
metric = 'mean'
metric2 = 'ARI'
cluster_num = 10
final_result <- vector()
cl2 = makeCluster(cluster_num)
registerDoParallel(cl2)
delta = 2; 
contam_rate = 0.2
d=3


for (p.index in 1:length(p.candidate)){
  for (r.index in 1:length(r.candidate)){
        combined_result <- foreach(i = 1:exp.time ) %dopar% {
          library(rTensor)
          library(ssvd)
          library(MASS)
          library(gtools)
          #library(cidr)
          library(mclust)
          library(pracma)

        r = rep(r.candidate[r.index],d)
        p = rep(p.candidate[p.index],d)
        data = TBM.generator(p, r, delta)
        # Lloyd use oracle initialization
        z.oracle.init = init.perturb(data$labels, contam_rate)
        z.oracle = HO.Lloyd(data$tensor, z.oracle.init)
        W.hat = toweightmatrix(z.oracle)
        Z.hat = tomembershipmatrix(z.oracle)
        X.hat= ttl(ttl(data$tensor, W.hat ,1:d), Z.hat, 1:d)
        # Lloyd use spectral initialization
        z.HOSC = HO.SC(data$tensor, r)
        z.Lloyd = HO.Lloyd(data$tensor, z.HOSC)
        W.hat.HOSC = toweightmatrix(z.Lloyd)
        Z.hat.HOSC = tomembershipmatrix(z.Lloyd)
        X.hat.HOSC = ttl(ttl(data$tensor, W.hat.HOSC ,1:d), Z.hat.HOSC, 1:d)

        if (strcmp(metric2,'ARI')){
          oracle_lloyd = nARI(z.oracle, data$labels)
          spect_lloyd = nARI(z.Lloyd, data$labels)
        } else{
          oracle_lloyd = nMCR(z.oracle, data$labels)
          spect_lloyd = nMCR(z.Lloyd, data$labels)
        }        
        HOOI_result = tucker(data$tensor, r)
        X.HOOI = ttl(HOOI_result$Z,HOOI_result$U,1:d) 
        oLloyd_error = rTensor::fnorm(X.hat - data$signal)
        sLloyd_error = rTensor::fnorm(X.hat.HOSC - data$signal)
        HOOI_error = rTensor::fnorm(X.HOOI - data$signal)
        error = c(oLloyd_error, sLloyd_error, HOOI_error, oracle_lloyd, spect_lloyd)
        return(error)
    }
    combined_result <- sapply(combined_result, function(err) err)
    if (strcmp(metric,'mean')){
      average = apply(combined_result,1,mean)
    } else{
      average = apply(combined_result,1,median)
    }
    sd_res = apply(combined_result,1,sd)
    per_result = c(p.candidate[p.index], r.candidate[r.index], average,sd_res)
    print(per_result)
    final_result = rbind(final_result, per_result)
  }
}

stopCluster(cl2)

colnames(final_result) <- c('p','r','oHLloyd','HLloyd','HOOI','oinit','spectinit','olloyd_sd','slloyd_sd','hooi_sd','oinit_sd','spectinit_sd')
file_name <- paste("tensor_recovery_largesep_winit_clu","p", paste(p.candidate, sep="", collapse="_"), "r", paste(r.candidate, sep="", collapse="_"),"d", d, "simu", exp.time, "metric", metric, "metric2", metric2 ,"delta", delta,'contam_rate', contam_rate,sep = "_")
print(file_name)
print(final_result)
write.csv(final_result, file = paste("../results/",file_name, ".csv", sep = ""))
