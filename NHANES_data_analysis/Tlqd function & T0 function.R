### log quantile density
Tlqd = function(corr_i, nb = 19){        
    dens.eps = 1e-6
    data_range = range(corr_i)
    basis1 = fda::create.bspline.basis(data_range, nb)
    
    densf = fda::density.fd(corr_i, fd(matrix(0,nb,1), basis1))
    di = list(x=seq(data_range[1],data_range[2],length.out=1024), 
              y=c(exp(eval.fd(x,densf$Wfdobj))/densf$C)) 
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
  })
}


###T0 function
T0 = function(corr_i, nb = 19){ 
  dens.eps = 1e-6
  data_range = range(corr_i)
  basis1 = fda::create.bspline.basis(data_range, nb)
  
  densf = fda::density.fd(corr_i, fd(matrix(0,nb,1), basis1))
  di = list(x=seq(data_range[1],data_range[2],length.out=1024), 
            y=c(exp(eval.fd(x,densf$Wfdobj))/densf$C))  
  ind1 = min(which(di$y > dens.eps*512)); ind2 = max(which(di$y > dens.eps*512))
  list(x = di$x[ind1:ind2], y = di$y[ind1:ind2])
}
