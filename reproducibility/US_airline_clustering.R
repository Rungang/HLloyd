source('HOLloyd.R')
load('US_flight_route.RData')

set.seed(2020)
r.candidate <- c(3)

r.var <- select_rank(as.tensor(US_tensor_data), r.candidate, US_pick_airline_name, US_pick_airport_name)
r.var
r_select <- r.var[which.min(r.var[,5] ),1:3] # r_select = 3, otherwise it will have a singleton
r_select

US_tensor_data <- as.tensor(US_tensor_data)
zus.HOSC <- HO.SC(US_tensor_data,r_select)
zus.Lloyd.HOSC = HO.Lloyd(US_tensor_data, zus.HOSC)


usairline_clu_result <- vector()
for (i in c(1:r_select[1])){
  cluster <- US_airline_infor[US_airline_infor$V4 %in% US_pick_airline_name[zus.Lloyd.HOSC[[1]]==i],]
  cluster$clu_num <- i
  print(cluster)
  usairline_clu_result <- rbind(usairline_clu_result, cluster)
}

usairport_clu_result <- vector()
for (i in c(1:r_select[2])){
  cluster <- US_airport_info[US_airport_info$V5 %in%US_pick_airport_name[zus.Lloyd.HOSC[[2]] == i],  ]
  cluster$clu_num <- i
  print(cluster)
  usairport_clu_result <- rbind(usairport_clu_result, cluster)
}


linefile_name <- paste("us_airline_clustering","line_num", length(US_pick_airline_name),"port_num", length(US_pick_airport_name), "select_r", paste(r_select, sep="", collapse="_"),sep = "_")
print(linefile_name)
portfile_name <- paste("us_airport_clustering","line_num", length(US_pick_airline_name),"port_num", length(US_pick_airport_name), "select_r", paste(r_select, sep="", collapse="_"),sep = "_")
print(portfile_name)
write.csv(usairline_clu_result, file = paste("../results/",linefile_name, ".csv", sep = ""))
write.csv(usairport_clu_result, file = paste("../results/",portfile_name, ".csv", sep = ""))

