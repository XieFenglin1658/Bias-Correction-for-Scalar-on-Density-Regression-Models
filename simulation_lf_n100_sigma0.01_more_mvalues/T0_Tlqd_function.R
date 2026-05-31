.libPaths(c("/insomnia001/home/fx2212/miniconda3/envs/fda-43/lib/R/library",
            "/insomnia001/home/fx2212/R/library"))

library(readxl)
library(magrittr)
library(dplyr)
library(fda)

### original density
T0 = function(corr_i, nb = 19){
    dens.eps = 1e-6
    basis1 = fda::create.bspline.basis(c(0,1), nb)
    
    densf = fda::density.fd(corr_i, fda::fd(matrix(0,nb,1), basis1))
    di = list(x=seq(0,1,length.out=1024),
              y=c(exp(eval.fd(seq(0,1,length.out=1024),densf$Wfdobj))/densf$C))

    ind1 = min(which(di$y > dens.eps*512)); ind2 = max(which(di$y > dens.eps*512))
    list(x = di$x[ind1:ind2], y = di$y[ind1:ind2])
}


### log quantile density
Tlqd = function(corr_i, nb = 19){
    dens.eps = 1e-6
    basis1 = fda::create.bspline.basis(c(0,1), nb)
    
    densf = fda::density.fd(corr_i, fda::fd(matrix(0,nb,1), basis1))
    di = list(x=seq(0,1,length.out=1024),
              y=c(exp(eval.fd(seq(0,1,length.out=1024),densf$Wfdobj))/densf$C))
    ind1 = min(which(di$y > dens.eps*512)); ind2 = max(which(di$y > dens.eps*512))
    di = list(x = di$x[ind1:ind2], y = di$y[ind1:ind2])
    
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