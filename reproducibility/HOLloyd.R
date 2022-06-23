library(rTensor)
library(ssvd)
library(MASS)
library(gtools)
#library(cidr)
library(mclust)
library(pracma)


TBM.generator <- function(p, r, delta, sigma = 1){
  # Generate core tensor.
  d = length(p)
  S.array = array(rnorm(prod(r)), dim=r) # Core tensor
  S = as.tensor(S.array)
  # adjust minimal separation.
  #delta.min = 100000
  Sk = k_unfold(S, 1)@data
  delta.min = sqrt(sum((Sk[1,] - Sk[2,])^2))
  for (i in 1:d){
    Sk = k_unfold(S, i)@data
    for (k1 in 1:(r[i]-1)){
      for (k2 in (k1+1):r[i]){
        delta.min = min(delta.min, sqrt(sum((Sk[k1,] - Sk[k2,])^2)))
      }
    }
  }
  S = S * delta / delta.min
  
  # Generate random labels 
  Z = list()
  z = list()
  for (i in 1:d){
    Zi = matrix(0, p[i], r[i])
    zi = rep(1:r[i], each = ceil(p[i]/r[i]))
    zi = permute(zi[1:p[i]])
    for (j in 1:p[i]){
      Zi[j, zi[j]] = 1
    }
    Z = c(Z, list(Zi))
    z = c(z, list(zi))
  }
  X = ttl(S, Z, 1:d)
  E = as.tensor(array(rnorm(prod(p)), dim=p)) # Noise tensor
  Y = X + sigma * E
  return(list("tensor" = Y, "labels" = z, "signal" = X))
}

TBM.generator.bernoulli <- function(p, r, scale, sigma = 1){
  # Generate core tensor.
  d = length(p)
  S.array = array(runif(prod(r)), dim=r) # Core tensor
  S = as.tensor(S.array)
  # do the scaling
  S = S / scale 
  
  # Generate random labels 
  Z = list()
  z = list()
  for (i in 1:d){
    Zi = matrix(0, p[i], r[i])
    zi = rep(1:r[i], each = ceil(p[i]/r[i]))
    zi = permute(zi[1:p[i]])
    for (j in 1:p[i]){
      Zi[j, zi[j]] = 1
    }
    Z = c(Z, list(Zi))
    z = c(z, list(zi))
  }
  X = ttl(S, Z, 1:d)
  Y = as.tensor(array(rbinom(prod(p),1,c(X@data)),dim = p))
  return(list("tensor" = Y, "labels" = z, "signal" = X))
}



TBM.generator.imbalance <- function(p, r, delta, ratio, sigma = 1){
  # Generate core tensor.
  d = length(p)
  S.array = array(rnorm(prod(r)), dim=r) # Core tensor
  S = as.tensor(S.array)
  # adjust minimal separation.
  #delta.min = 100000
  Sk = k_unfold(S, 1)@data
  delta.min = sqrt(sum((Sk[1,] - Sk[2,])^2))
  for (i in 1:d){
    Sk = k_unfold(S, i)@data
    for (k1 in 1:(r[i]-1)){
      for (k2 in (k1+1):r[i]){
        delta.min = min(delta.min, sqrt(sum((Sk[k1,] - Sk[k2,])^2)))
      }
    }
  }
  S = S * delta / delta.min
  
  # Generate random labels 
  Z = list()
  z = list()
  for (i in 1:d){
    zi <- vector()
    Zi = matrix(0, p[i], r[i])
    for (index in 1:r[i]){
      zi <- c(zi, rep(index, each = ceil(p[i]* ratio[index])))
    }
    zi = permute(zi[1:p[i]])
    for (j in 1:p[i]){
      Zi[j, zi[j]] = 1
    }
    Z = c(Z, list(Zi))
    z = c(z, list(zi))
  }
  X = ttl(S, Z, 1:d)
  E = as.tensor(array(rnorm(prod(p)), dim=p)) # Noise tensor
  Y = X + sigma * E
  return(list("tensor" = Y, "labels" = z, "signal" = X))
}

 
toweightmatrix <- function(z){
  d = length(z)
  W = list()
  for (i in 1:d){
    zi = z[[i]]
    pi = length(zi)
    ri = length(unique(zi))
    Wi = matrix(0, ri, pi)
    for (k in 1:ri){
      zi.subset = which(zi==k)
      Wi[k,zi.subset] = 1/length(zi.subset)
    }
    W = c(W, list(Wi))
  }
  return(W)
}

tomembershipmatrix <- function(z){
  d = length(z)
  Z = list()
  for (i in 1:d){
    zi = z[[i]]
    pi = length(zi)
    ri = length(unique(zi))
    Zi = matrix(0, pi, ri)
    for (k in 1:ri){
      zi.subset = which(zi==k)
      Zi[zi.subset,k] = 1
    }
    Z = c(Z, list(Zi))
  }
  return(Z)
}

## This is the spectral clustering in the matrix case. This is also the spectral clustering based on HOSVD in the tensor case.
SC <- function(Y, r){
  try(if(missing("Y")) stop("missing argument: Y is required as the tensor type."))
  try(if(missing("r")) stop("missing argument: r is required as a scalar or vector."))
  try(if(class(Y) != "Tensor") stop("invalid input: Y should be of tensor type."))
  p = dim(Y)
  d = length(p)
  
  if(is.atomic(r) && length(r)==1){
    r = rep(r, d)
  }
  
  z_0 = list();
  for (i in 1:d){
    MY = k_unfold(Y, i)@data
    V_0 = svd(MY)$v[,1:r[i]]
    A.bar = MY %*% V_0
    z_0 = c(z_0, list(kmeanspp(A.bar, r[i])$cluster))
  }

  
  return(z_0)
}

# ## High order clustering based on HOSVD
# HOSVD.SC <- function(Y, r){
#   try(if(missing("Y")) stop("missing argument: Y is required as the tensor type."))
#   try(if(missing("r")) stop("missing argument: r is required as a scalar or vector."))
#   try(if(class(Y) != "Tensor") stop("invalid input: Y should be of tensor type."))
#   p = dim(Y)
#   d = length(p)
  
#   if(is.atomic(r) && length(r)==1){
#     r = rep(r, d)
#   }
  
#   U_t = list();
#   for (i in 1:d){#Initialization
#     MY = k_unfold(Y, i)@data
#     U_t = c(U_t, list(t(svd(MY%*%t(MY))$u[,1:r[i]])))
#   }
  
#   z_0 = list();
#   for (i in 1:d){
#     # A.tilde = k_unfold(S.hat, i)@data
#     # A.bar = t(U_t[[i]]) %*% A.tilde
#     # A.bar = A.bar %*% svd(A.bar)$v[,1:r[i]]
#     Kmatrix = ttl(Y,U_t[-i],(1:d)[-i])
#     Kmatrix = t(U_t[[i]]) %*% U_t[[i]] %*% k_unfold(Kmatrix,i)@data
#     z_0 = c(z_0, list(kmeanspp(Kmatrix, r[i])$cluster))
#   }
  
#   return(z_0)
# }



## One-step HOOI initialization
HO.SC <- function(Y, r){
  try(if(missing("Y")) stop("missing argument: Y is required as the tensor type."))
  try(if(missing("r")) stop("missing argument: r is required as a scalar or vector."))
  try(if(class(Y) != "Tensor") stop("invalid input: Y should be of tensor type."))
  p = dim(Y)
  d = length(p)
  
  if(is.atomic(r) && length(r)==1){
    r = rep(r, d)
  }
  
  U_t = list();
  for (i in 1:d){#Initialization
    MY = k_unfold(Y, i)@data
    U_t = c(U_t, list(t(svd(MY%*%t(MY))$u[,1:r[i]])))
  }
  
 # t = 1; tmax = 5;
  
  for(i in 1:d){
    A = ttl(Y, U_t[-i], (1:d)[-i])
    A_matrix = k_unfold(A, i)@data
    svd.result = svd(A_matrix)
    U_t[[i]] = t(svd.result$u[,1:r[i]])
  }

  #S.hat = ttl(Y, U_t, (1:d))

  
  z_0 = list();
  for (i in 1:d){
    # A.tilde = k_unfold(S.hat, i)@data
    # A.bar = t(U_t[[i]]) %*% A.tilde
    # A.bar = A.bar %*% svd(A.bar)$v[,1:r[i]]
    Kmatrix = ttl(Y,U_t[-i],(1:d)[-i])
    Kmatrix = t(U_t[[i]]) %*% U_t[[i]] %*% k_unfold(Kmatrix,i)@data
    z_0 = c(z_0, list(kmeanspp(Kmatrix, r[i])$cluster))
  }
  
  return(z_0)
}

# use HOOI as the initialization of spectral clustering
HOOI.SC <- function(Y, r){
  try(if(missing("Y")) stop("missing argument: Y is required as the tensor type."))
  try(if(missing("r")) stop("missing argument: r is required as a scalar or vector."))
  try(if(class(Y) != "Tensor") stop("invalid input: Y should be of tensor type."))
  p = dim(Y)
  d = length(p)
  
  if(is.atomic(r) && length(r)==1){
    r = rep(r, d)
  }
  HOOI_result = tucker(Y, r)
  U = HOOI_result$U

  U_t = list();
  for (i in 1:d){#Initialization
    U_t[[i]] = t(U[[i]])
  }
  z_0 = list();
  for (i in 1:d){
    # Kmatrix = ttl(Y,U_t[-i],(1:d)[-i])
    # Kmatrix = t(U_t[[i]]) %*% U_t[[i]] %*% k_unfold(Kmatrix,i)@data
    # z_0 = c(z_0, list(kmeanspp(Kmatrix, r[i])$cluster))
    
    # directly use singular space
    z_0 = c(z_0, list(kmeanspp(t(U_t[[i]]), r[i])$cluster))
  }
  return(z_0)
}


# use CP decomposition as the initialization for spectral clustering
CP.SC <- function(Y, r){
  try(if(missing("Y")) stop("missing argument: Y is required as the tensor type."))
  try(if(missing("r")) stop("missing argument: r is required as a scalar or vector."))
  try(if(class(Y) != "Tensor") stop("invalid input: Y should be of tensor type."))
  p = dim(Y)
  d = length(p)
  
  if(is.atomic(r) && length(r)==1){
    r = rep(r, d)
  }
  CP_result = cp(Y, num_components = r[1])
  U = CP_result$U

  U_t = list();
  for (i in 1:d){#Initialization
    U_t[[i]] = t(U[[i]])
  }
  z_0 = list();
  for (i in 1:d){
    # Kmatrix = ttl(Y,U_t[-i],(1:d)[-i])
    # Kmatrix = t(U_t[[i]]) %*% U_t[[i]] %*% k_unfold(Kmatrix,i)@data
    # z_0 = c(z_0, list(kmeanspp(Kmatrix, r[i])$cluster))
    
    # directly use singular space
    z_0 = c(z_0, list(kmeanspp(t(U_t[[i]]), r[i])$cluster))
  }
  return(z_0)
}





######################
# High-order Lloyd Algorithm
######################
HO.Lloyd <- function(Y, z, t_max = 10){
  d = length(z)
  r = rep(0, d)
  for (i in 1:d){
    r[i] = length(unique(z[[i]]))
  }
  for (iter in 1:t_max){
    # create weight matrix
    W = toweightmatrix(z)
    z = list()
    S.hat = ttl(Y, W, (1:d))
    for (i in 1:d){
      S.hat.matrix = k_unfold(S.hat, i)@data
      A.matrix = k_unfold(ttl(Y, W[-i], (1:d)[-i]) , i)@data
      zi = rep(0, nrow(A.matrix))
      for (j in 1:nrow(A.matrix)){
        dist = 1000000000
        nearest = 1
        for (k in 1:nrow(S.hat.matrix)){
          dist_k = sqrt(sum((A.matrix[j,] - S.hat.matrix[k,])^2))
          if(dist_k < dist){
            nearest = k
            dist = dist_k
          }
        }
        zi[j] = nearest
      }
      z = c(z, list(zi))
    }
  }
  return(z) 
}


######################
# Modified High-order Lloyd Algorithm
######################
Modi.HO.Lloyd <- function(Y, z, t_max = 10){
  d = length(z)
  r = rep(0, d)
  for (i in 1:d){
    r[i] = length(unique(z[[i]]))
  }
  for (iter in 1:t_max){
    # create weight matrix
    Z = tomembershipmatrix(z)
    W = toweightmatrix(z)
    z = list()
    S.hat = ttl(Y, W, (1:d))
    for (i in 1:d){
      S.hat.matrix = k_unfold(ttl(S.hat,Z[-i],(1:d)[-i]), i)@data
      Y.matrix = k_unfold(Y, i)@data
      zi = rep(0, nrow(Y.matrix))
      for (j in 1:nrow(Y.matrix)){
        dist = 1000000000
        nearest = 1
        for (k in 1:nrow(S.hat.matrix)){
          dist_k = sqrt(sum((Y.matrix[j,] - S.hat.matrix[k,])^2))
          if(dist_k < dist){
            nearest = k
            dist = dist_k
          }
        }
        zi[j] = nearest
      }
      z = c(z, list(zi))
    }
  }
  return(z) 
}




## calculate the missclassfication rate of two sets of labels.
ARI <- function(z,z.hat){
  d = length(z)
  ari = 1
  for (i in 1:d){
    ari = min(ari, adjustedRandIndex(z[[i]],z.hat[[i]]))
  }
  return(ari)
}


## New ARI computation: compute the average ARI
nARI <- function(z,z.hat){
  d = length(z)
  ari = rep(0,d)
  for (i in 1:d){
    ari[i] = adjustedRandIndex(z[[i]],z.hat[[i]])
  }
  ave_ari = sum(ari)/d
  return(ave_ari)
}

 

## calculate the missclassfication rate of two sets of labels.
MCR <- function(z,z.hat){
  d = length(z)
  mcr = 0
  for (i in 1:d){
    ri = length(unique(z.hat[[i]]))
    perm.array = permutations(ri,ri)
    MCR.temp = 1
    for (j in 1:nrow(perm.array)){
      z.perm = z.hat[[i]]
      for (k in 1:ri){
        z.perm[which(z.hat[[i]]==k)] = perm.array[j,k]
      }
      MCR.temp = min(MCR.temp, sum(z.perm!=z[[i]])/length(z[[i]]))
    }
    mcr = max(mcr, MCR.temp)
  }
  
  return(mcr)
}

## new MCR computation: compute the average MCR along all modes
nMCR <- function(z,z.hat){
  d = length(z)
  mcr = rep(0,d)
  for (i in 1:d){
    ri = length(unique(z.hat[[i]]))
    perm.array = permutations(ri,ri)
    MCR.temp = 1
    for (j in 1:nrow(perm.array)){
      z.perm = z.hat[[i]]
      for (k in 1:ri){
        z.perm[which(z.hat[[i]]==k)] = perm.array[j,k]
      }
      MCR.temp = min(MCR.temp, sum(z.perm!=z[[i]])/length(z[[i]]))
    }
    mcr[i] = MCR.temp
  }
  ave_mcr = sum(mcr)/d
  return(ave_mcr)
}





## Randomly perturb the labels to generate an oracle version.
init.perturb <- function(z, p){
  d = length(z)
  z.perturb = list()
  for (i in 1:d){
    z.temp = z[[i]]
    flag = runif(length(z[[i]]))
    for (j in 1:length(z.temp)){
      if(flag[j]<=p){
        z.temp[j] = permute(setdiff(1:r[i], z[[i]][j]))[1]
      }
    }
    z.perturb = c(z.perturb, list(z.temp))
  }
  
  return(z.perturb)
}

## Used to select the rank in the real data
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

