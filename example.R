library(HLloyd)

# Simulated data.
set.seed(100)
r = c(5,5,5) # num of clusters
p = c(50,50,50) # tensor dimension
delta = 0.5 # slice seperation
data = TBM.generator(p, r, delta) # generate the model

## matrix method
z.SC = SC(data$tensor, r) 
print(sprintf("Averaged ARI from matrix method is %f", ARI(z.SC, data$labels, mode="averaged")))

## tensor method
z.HOSC = HO.SC(data$tensor, r) # high-order initialization
z.Lloyd.HOSC = HO.Lloyd(data$tensor, z.HOSC) # HLloyd iteration

print(sprintf("Averaged ARI from matrix method is %f", ARI(z.Lloyd.HOSC, data$labels, mode="averaged")))


