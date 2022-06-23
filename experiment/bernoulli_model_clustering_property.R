# Study the consistency property of the HSC+HLloyd method in the Bernoulli model
source("HOLloyd.R")
library(doParallel)
set.seed(2020)
server = 1
if (server) {
	exp.time = 100 # runing time: 340 s
	cluster_num = 10
} else{
	exp.time = 10
	cluster_num = 2
}
d = 4
r.candidate = c(3,5,7)
p.candidate = c(40,50,60,70,80,90,100)
scale = 10
metric = 'ARI'

final_result <- vector()
cl2 = makeCluster(cluster_num)
registerDoParallel(cl2)


for (p.index in 1:length(p.candidate)){
	for (r.index in 1:length(r.candidate)){
		combined_result <- foreach(i = 1:exp.time) %dopar% {
        library(rTensor)
        library(ssvd)
        library(MASS)
        library(gtools)
        #library(cidr)
        library(mclust)
        library(pracma)
          r = rep(r.candidate[r.index],d)
		  p = rep(p.candidate[p.index],d)
		  data = TBM.generator.bernoulli(p, r, scale)
		  z.HOSC = HO.SC(data$tensor, r)

		  z.Lloyd.HOSC = HO.Lloyd(data$tensor, z.HOSC)
		  if (strcmp(metric,'ARI')){
		  	kmeans_acc = 1-nARI(z.HOSC, data$labels)
	        Lloyd_acc = 1-nARI(z.Lloyd.HOSC, data$labels)
	      } else{
	      	kmeans_acc = nMCR(z.HOSC, data$labels)
	        Lloyd_acc = nMCR(z.Lloyd.HOSC, data$labels)
	      }
	      error = c(Lloyd_acc, kmeans_acc)
	      return(error)
		}
		combined_result <- sapply(combined_result, function(err) err)
	    mean_res = apply(combined_result,1,mean)
	    sd_res = apply(combined_result,1,sd)
	    per_result = c(p.candidate[p.index], r.candidate[r.index], mean_res,sd_res)
	    print(per_result)
	    final_result = rbind(final_result, per_result)
	}
}

stopCluster(cl2)


colnames(final_result) <- c('p','r','HLloyd','Kmeans','L_sd','K_sd')
file_name <- paste("bernoulli_model_clustering","p", paste(p.candidate, sep="", collapse="_"), "r", paste(r.candidate, sep="", collapse="_"), "simu", exp.time, 'd', d ,'metric',metric,'scale', scale,sep = "_")
print(file_name)
print(final_result)
write.csv(final_result, file = paste("../results/",file_name, ".csv", sep = ""))