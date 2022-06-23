# Experiment 2-1.  Phase transition in Matrix.
# Compare oracle initialization and spectral intialization.

source("HOLloyd.R")
library(doParallel)
set.seed(100)
r_total = 5
r = rep(r_total,2)
p.candidate = c(100,200,400)
delta.candidate = seq(0,1,0.05)
metric = 'MCR'
server = 0
if (server) {
  exp.time = 120 # runing time: 215 s
  cluster_num = 40
} else{
  exp.time = 10
  cluster_num = 2
}
final_result <- vector()
cl2 = makeCluster(cluster_num)
registerDoParallel(cl2)


for (p.index in 1:length(p.candidate)){
  for (delta.index in 1:length(delta.candidate)){
    combined_result <- foreach(i = 1:exp.time ) %dopar% {
      library(rTensor)
      library(ssvd)
      library(MASS)
      library(gtools)
      #library(cidr)
      library(mclust)
      library(pracma)
      p = rep(p.candidate[p.index],2)
      data = TBM.generator(p, r, p[1]^(-delta.candidate[delta.index]))
      z.SC = SC(data$tensor, r)
      z.Lloyd = HO.Lloyd(data$tensor, z.SC)
      # slightly perturb the labels.
      z.oracle.init = init.perturb(data$labels, 1/r)
      z.oracle = HO.Lloyd(data$tensor, z.oracle.init)
      if (strcmp(metric,'ARI')){
        Lloyd_error = nARI(z.Lloyd, data$labels)
        oracle_error = nARI(z.oracle, data$labels)
      } else{
        Lloyd_error = nMCR(z.Lloyd, data$labels)
        oracle_error = nMCR(z.oracle, data$labels)
      }
      error = c(Lloyd_error, oracle_error)
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


colnames(final_result) <- c('p','alpha','Lloyd',"Oracle",'L_sd', 'O_sd')
file_name <- paste("matrix_phase_tran","p", paste(p.candidate, sep="", collapse="_"), "r", paste(r_total, sep="", collapse="_"), "simu", exp.time,'metric',metric,sep = "_")
print(file_name)
print(final_result)
write.csv(final_result, file = paste("../results/",file_name, ".csv", sep = ""))



# write.table(ARI.oracle, "matrix_pt_oracle.txt", row.names = FALSE, col.names = FALSE)
# write.table(ARI.Lloyd, "matrix_pt_Lloyd.txt", row.names = FALSE, col.names = FALSE)
