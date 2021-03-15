#' @import rTensor
#' @import gtools
#' @import LICORS
#' @import mclust
NULL




#' Generate a tensor block model under permutation and noises.
#' @param p A vector of tensor dimensions.
#' @param r A vector of cluster numbers in each mode.
#' @param delta The slice separation of core tensor.
#' @param sigma Observational noise level.
#' @return Observational tensor; clutering labels; signal tensor.
#' @export
TBM.generator <- function(p, r, delta = 1, sigma = 1){
  try(if(length(p)!=length(r) | any(r>p)) stop("invalid input: Y should be of tensor type."))

  try(if(delta <= 0) stop("Invalid delta."))

  # Generate core tensor.
  d = length(p)
  S.array = array(rnorm(prod(r)), dim=r) # Core tensor
  S = as.tensor(S.array)
  # adjust minimal separation.
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
    zi = rep(1:r[i], each = ceiling(p[i]/r[i]))
    zi = c(zi, rep(1, p[i]-length(zi)))
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



#' Transform a label vector to a weighted matrix
#' Non-exported.
#' @param z vector of clustering labels
#' @return Corresponding weighted matrix
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

#' Transform a label vector to a membership matrix.
#' Non-exported.
#' @param z vector of clustering labels.
#' @return Corresponding membership matrix.
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

#' Matricized tensor spectral clustering using kmeans++.
#' @param Y Observed tensor.
#' @param r Vector of cluster numbers.
#' @return clustering labels.
#' @export
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

#' High-order tensor spectral clustering.
#' @param Y Observed tensor.
#' @param r Vector of cluster numbers.
#' @return clustering labels.
#' @export
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

  for(i in 1:d){#Projection
    A = ttl(Y, U_t[-i], (1:d)[-i])
    A_matrix = k_unfold(A, i)@data
    svd.result = svd(A_matrix)
    U_t[[i]] = t(svd.result$u[,1:r[i]])
  }

  z_0 = list();
  for (i in 1:d){
    Kmatrix = ttl(Y,U_t[-i],(1:d)[-i])
    Kmatrix = t(U_t[[i]]) %*% U_t[[i]] %*% k_unfold(Kmatrix,i)@data
    z_0 = c(z_0, list(kmeanspp(Kmatrix, r[i])$cluster))
  }

  return(z_0)
}


#' High-order Lloyd Algorithm.
#' @param Y Observed tensor.
#' @param z Vector of initialized labels.
#' @param t_max maximum number of iterations.
#' @return Estimated clustering labels.
#' @export
HO.Lloyd <- function(Y, z, t_max = 10){
  try(if(missing("Y")) stop("missing argument: Y is required as the tensor type."))
  try(if(missing("z")) stop("missing argument: z is required as a vector."))
  try(if(class(Y) != "Tensor") stop("invalid input: Y should be of tensor type."))
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
        dist = Inf
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

#' Calculate ARI of two sets of labels.
#' @param z First clustering vector.
#' @param z.hat Second clustering vector.
#' @param mode Averaged/Minimal ARI.
#' @return Calculated ARI.
ARI <- function(z,z.hat,mode=c("averaged","minimal")){
  d = length(z)
  if (mode == "averaged"){
    ari = 1
    for (i in 1:d){
      ari = min(ari, adjustedRandIndex(z[[i]],z.hat[[i]]))
    }
  }
  else{
    ari = rep(0,d)
    for (i in 1:d){
      ari[i] = adjustedRandIndex(z[[i]],z.hat[[i]])
    }
    ari = sum(ari)/d
  }
  return(ari)
}
