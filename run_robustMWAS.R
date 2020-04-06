
run_robustMWAS = function(param) {
 
    nfolds = param$nfolds_robustMWAS
    param$covariates                = data.frame( param$covariates[,colnames(param$covariates)[1]],param$covariates[,param$modeloutcome],param$covariates[,param$modelcovariates])
    if ( param$modelPCs > 0 ) param$covariates = add_PCs(param,rowsubset=NULL)
    colnames(param$covariates)[1:2] = c(samples_var,param$modeloutcome)
    param$modelcovariates           =  colnames(param$covariates)[-(1:2)] 
      
    param$modelPCs                = 0
  
    for (j in 1:nfolds) { # j=1
      
      param_robust = param #  param_robust = param
      param_robust$dircoveragenorm =  paste0(param$dircoveragenorm_robust,"/fold_",j,"/")
      
      samples                 = readLines( paste0(param_robust$dircoveragenorm,"Coverage.nmsrow.txt") )
      param_robust$dirmwas    = paste0(param$dirmwas_robust,"/fold_",j,"/")
      
      param_robust$dirmwas   =  paste0(param$dirmwas,"/robust_",param$nfolds_robustMWAS,"_folds_",param$residualize_robustMWAS,"/fold_",j,"/")
      dir.create(param_robust$dirmwas,showWarnings=F,recursive=T)
      
      ind                     = match(samples,param$covariates[,1],nomatch = 0L)
      param_robust$covariates = param$covariates[ind,]
      
      if (param$residualize_robustMWAS=="resid") {
        ind                     = match(samples,rownames(cvrtqr),nomatch = 0L)
        cvrtqr_robust           = cvrtqr[ind,]
        param_robust$modelcovariates    = NULL
        data  = as.matrix ( param_robust$covariates[,param_robust$modeloutcome,drop=F] )
        param_robust$covariates[,param_robust$modeloutcome] = data - tcrossprod(cvrtqr_robust, crossprod(data, cvrtqr_robust))
      } 
      
      file = paste0(param_robust$dirmwas,"/QQ_plot.pdf")   
      if (!file.exists( file)) ramwas5MWAS(param_robust) else message("Robust bulk MWAS was already run, use the 2nd optional \"rerun\" argument for rerunning")
  } 

}  