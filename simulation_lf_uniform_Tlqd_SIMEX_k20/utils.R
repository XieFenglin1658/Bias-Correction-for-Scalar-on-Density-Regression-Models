# set library path
.libPaths(c("/insomnia001/home/fx2212/miniconda3/envs/fda-43/lib/R/library",
            "/insomnia001/home/fx2212/R/library"))

library(readxl)
library(magrittr)
library(dplyr)


## whether a connectivity matrix qualified for analysis
qualify = function(x){
  if(is.null(x)) return(FALSE)
  else{
    m = x[upper.tri(x)]
    m = m[!is.na(m)]
    if(length(m) < 3003) return(FALSE)
    else return(TRUE)
  }
}


## transformed density estimation
library(fda)
nb = 19
basis1 = create.bspline.basis(c(-1,1), nb)
basis2 = create.bspline.basis(c(0,1), nb)

### original density
dens.eps = 1e-6 ## tolerance of density, do not consider density below it
T0 = function(corr_i, method="spline"){
  if(method == "kde"){
    di = density(corr_i, kernel="gaussian", from=-1, to=1, n=512)
  } else{
    densf = fda::density.fd(corr_i, fd(matrix(0,nb,1), basis1))
    di = list(x=seq(-1,1,length.out=512),
              y=c(exp(eval.fd(seq(-1,1,length.out=512),densf$Wfdobj))/densf$C))
  }
  list(x=di$x, y=di$y)
}

### log density
Tl = function(corr_i, method="spline"){
  if(method == "kde"){
    di = density(corr_i, kernel="gaussian", from=-1, to=1, n=512)
    return(list(x=di$x, y=log(dens.eps+di$y)))
  } else{
    densf = fda::density.fd(corr_i, fd(matrix(0,nb,1), basis1))
    di = list(x=seq(-1,1,length.out=512),
              y=c(eval.fd(seq(-1,1,length.out=512),densf$Wfdobj)-log(densf$C)))
    return(list(x=di$x, y=di$y))
  }
}

### quantile function
Tq = function(corr_i){
  x = seq(0,1,length.out=512)
  y = quantile(corr_i, x)
  list(x=x, y=y)
}

### log quantile density
Tlqd = function(corr_i, method="spline"){
  if(method == "kde"){
    di = density(corr_i, kernel="gaussian", from=min(corr_i)-0.1, to=1, n=512)
  } else{
    densf = fda::density.fd(corr_i, fd(matrix(0,nb,1), basis2))
    di = list(x=seq(0,1,length.out=1024),
              y=c(exp(eval.fd(seq(0,1,length.out=1024),densf$Wfdobj))/densf$C))
    ind1 = min(which(di$y > dens.eps*512)); ind2 = max(which(di$y > dens.eps*512))
    di = list(x = di$x[ind1:ind2], y = di$y[ind1:ind2])
  }
  tryCatch({
    y = suppressWarnings(-dens2lqd(di$y, di$x, N=512))
    x = seq(0,1,length.out=512)
    list(x=x, y=y)
  },
  error = function(e){
    di$y = di$y + dens.eps
    y = suppressWarnings(-dens2lqd(di$y, di$x, N=512))
    x = seq(0,1,length.out=512)
    list(x=x, y=y)
  }
  )
}
