# Experiment 2-2.  Phase transition in tensor.
# Compare oracle initialization and spectral intialization. This code is also used in studying the consistency and exact recovery of the HSC and HLoyd algorithm.

source("HOLloyd.R")
library(doParallel)
set.seed(2020)
server = 1
if (server) {
	exp.time = 40 # runing time: 340 s
	cluster_num = 10
} else{
	exp.time = 10
	cluster_num = 2
}

r_total = 5
d = 4
r = rep(r_total,d)
p.candidate = c(50,80,100)
delta_low = 0.7
delta_upper = 1.8
delta.candidate = seq(delta_low,delta_upper,0.03)
metric = 'ARI'

final_result <- vector()
cl2 = makeCluster(cluster_num)
registerDoParallel(cl2)


for (p.index in 1:length(p.candidate)){
	for (delta.index in 1:length(delta.candidate)){
		combined_result <- foreach(i = 1:exp.time) %dopar% {
        library(rTensor)
        library(ssvd)
        library(MASS)
        library(gtools)
        #library(cidr)
        library(mclust)
        library(pracma)
		  p = rep(p.candidate[p.index],d)
		  data = TBM.generator(p, r, 10*p[1]^(-delta.candidate[delta.index]))
		  #z.SC = SC(data$tensor, r)
		  z.HOSC = HO.SC(data$tensor, r)

		  # USE HOOI as initalization
		  #z.HOSC = HOOI.SC(data$tensor, r)
		  #z.Lloyd.SC = HO.Lloyd(data$tensor, z.SC)
		  #  MCR(z.Lloyd.SC, data$labels)
		  z.Lloyd.HOSC = HO.Lloyd(data$tensor, z.HOSC)
 		  z.oracle.init = init.perturb(data$labels, 0.2)
		  z.oracle = HO.Lloyd(data$tensor, z.oracle.init)
		  if (strcmp(metric,'ARI')){
		  	kmeans_error = nARI(z.HOSC, data$labels)
	        Lloyd_error = nARI(z.Lloyd.HOSC, data$labels)
	        oracle_error = nARI(z.oracle, data$labels)
	      } else{
	      	kmeans_error = nMCR(z.HOSC, data$labels)
	        Lloyd_error = nMCR(z.Lloyd.HOSC, data$labels)
	        oracle_error = nMCR(z.oracle, data$labels)
	      }
	      error = c(Lloyd_error, oracle_error,kmeans_error)
	      return(error)
		}
		combined_result <- sapply(combined_result, function(err) err)
	    mean_res = apply(combined_result,1,mean)
	    sd_res = apply(combined_result,1,sd)
	    per_result = c(p.candidate[p.index], delta.candidate[delta.index], mean_res,sd_res)
	    print(per_result)
	    final_result = rbind(final_result, per_result)
	}
}

stopCluster(cl2)


colnames(final_result) <- c('p','alpha','HLloyd',"Oracle",'Kmeans','L_sd', 'O_sd','K_sd')
file_name <- paste("tensor_phase_tran_wkmeans","p", paste(p.candidate, sep="", collapse="_"), "r", paste(r_total, sep="", collapse="_"), "simu", exp.time, 'd', d ,'metric',metric,'delta_low',delta_low,'delta_upper',delta_upper,sep = "_")
print(file_name)
print(final_result)
write.csv(final_result, file = paste("../results/",file_name, ".csv", sep = ""))

# write.table(ARI.oracle, "tensor_pt_oracle.txt", row.names = FALSE, col.names = FALSE)
# write.table(ARI.Lloyd.HOSC, "tensor_pt_Lloyd.txt", row.names = FALSE, col.names = FALSE)
# 