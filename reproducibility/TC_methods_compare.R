# In this simulation, we want to do experiment to compare different methods.
source("HOLloyd.R")
library(doParallel)
set.seed(2020)
server = 1
if (server) {
	exp.time = 100 # d=3 runing time: 300 s
	cluster_num = 10
} else{
	exp.time = 10    
	cluster_num = 2
}
r_ind = 5
d_total = c(4)
p_ind = 80
delta.candidate = seq(0.5,0.7,0.03) # d= 3's delta candidate seq(0.5,0.8,0.03)   d= 4's delta candidate seq(0.7,1,0.03)
metric = 'ARI'

final_result <- vector()
cl2 = makeCluster(cluster_num)
registerDoParallel(cl2)


for (d.index in 1:length(d_total)){
	for (delta.index in 1:length(delta.candidate)){
		combined_result <- foreach(i = 1:exp.time) %dopar% {
        library(rTensor)
        library(ssvd)
        library(MASS)
        library(gtools)
        #library(cidr)
        library(mclust)
        library(pracma)
		  p = rep(p_ind,d_total[d.index])
		  r = rep(r_ind,d_total[d.index])
		  data = TBM.generator(p, r, 10*p[1]^(-delta.candidate[delta.index]))
		  #z.SC = SC(data$tensor, r)
		  z.HOSC = HO.SC(data$tensor, r)
		  z.HOSVD = SC(data$tensor, r)
		  z.CP = CP.SC(data$tensor, r)
		  z.Lloyd.HOSC = HO.Lloyd(data$tensor, z.HOSC)
		  if (strcmp(metric,'ARI')){
		  	ohooi_metric = 1-nARI(z.HOSC, data$labels)
	        Lloyd_metric = 1-nARI(z.Lloyd.HOSC, data$labels)
	        hosvd_metric = 1-nARI(z.HOSVD, data$labels)
	        cp_metric = 1-nARI(z.CP, data$labels)
	      } else{
	      	ohooi_metric = nMCR(z.HOSC, data$labels)
	        Lloyd_metric = nMCR(z.Lloyd.HOSC, data$labels)
	        hosvd_metric = nMCR(z.HOSVD, data$labels)
	        cp_metric = nMCR(z.CP, data$labels)
	      }
		  error_metric = c(Lloyd_metric, ohooi_metric, hosvd_metric, cp_metric)
	      return(error_metric)
		}
		combined_result <- sapply(combined_result, function(err) err)
	    mean_res = apply(combined_result,1,mean)
	    sd_res = apply(combined_result,1,sd)
	    per_result = c(d_total[d.index], delta.candidate[delta.index], mean_res,sd_res)
	    print(per_result)
	    final_result = rbind(final_result, per_result)
	}
}

stopCluster(cl2)


colnames(final_result) <- c('d','alpha','HLloyd','ohooi', 'hosvd', 'cp','L_sd','ohooi_sd','hosvd_sd', 'cp_sd')
file_name <- paste("TC_method_compare","d", paste(d_total, sep="", collapse="_"),"p", paste(p_ind, sep="", collapse="_"), "r", paste(r_ind, sep="", collapse="_"), "simu", exp.time,'metric',metric, "delta_min", min(delta.candidate), "delta_max", max(delta.candidate),sep = "_")
print(file_name)
print(final_result)
write.csv(final_result, file = paste("../results/",file_name, ".csv", sep = ""))