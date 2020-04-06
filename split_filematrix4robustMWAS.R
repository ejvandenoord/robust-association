
split_filematrix4robustMWAS = function(param) {
  
  calc_foldid = function(param) {  
    
    n_folds = param$nfolds_robustMWAS
    
    calc_foldid_plain = function(n_folds,n_sample) {
      fold_size = rep(floor(n_sample/n_folds),n_folds)
      rest      = n_sample-sum(fold_size)
      if (rest>0) fold_size[1:rest] = fold_size[1:rest] + 1
      foldid    = sample(rep(1:n_folds,fold_size),n_sample) 
      foldid  
    }
    
    calc_foldid_statify = function(n_folds,n_sample,stratvar) {
      
      foldid       = rep(NA,n_sample)
      min_stratvar = min(stratvar)
      
      sel = stratvar == min(stratvar)
      foldid[which(stratvar == min_stratvar )] = calc_foldid_plain(n_folds,sum(sel))
      foldid[which(stratvar != min_stratvar )] = calc_foldid_plain(n_folds,sum(!sel))
      
      foldid  
    }
    
    
    if (param$mrs_gsms=="yes") {
      
      subjectnames        = unlist( lapply(strsplit( param$covariates[,1],"_"), function(l) l[[1]]) ) 
      unique_subjectnames = unique( subjectnames )
      
      foldid = calc_foldid_plain(n_folds,length(unique_subjectnames))
      ind    = match (subjectnames,unique_subjectnames)
      foldid =  foldid[ind]
      
    } else {  
      
      n_sample = length(param$covariates[,1])
      ncat     = length( unique(param$covariates[,param$modeloutcome]) )
      if ( ncat  == 2 ) foldid = calc_foldid_statify(n_folds,n_sample,param$covariates[,param$modeloutcome]) else 
        foldid = calc_foldid_plain(  n_folds,n_sample)
    }
    
    names(foldid) =  param$covariates[,1] 
    foldid 
  }
  # Get covariates + PCs matrix for analysis
  # orthonormalized unless normalize == FALSE
  .getCovariates_mod = function(
    param,
    rowsubset = NULL,
    normalize = TRUE,
    modelhasconstant){
    
    # rowsubset = NULL; normalize = TRUE; modelhasconstant = T
    
    # Named covariates
    cvrt = param$covariates[ param$modelcovariates ];
    
    ### Add PCs as covariates
    if( param$modelPCs > 0 ){
      # e = readRDS(paste0(param$dirpca,"/eigen.rds"));
      eigenvectors = fm.open(
        filenamebase = paste0(param$dirpca, "/eigenvectors"),
        readonly = TRUE);
      PCs = eigenvectors[, seq_len(param$modelPCs)]; # , drop=FALSE
      close(eigenvectors);
      
      if(!is.null( rowsubset ))
        PCs = PCs[rowsubset,];
      mwascvrtqr = cbind(cvrt, PCs);
      colnames(  mwascvrtqr ) = c( param$modelcovariates, paste0("PC",1:param$modelPCs) ) # + EvdO
      # rm(e);
    } else {
      mwascvrtqr = cvrt;
      colnames(  mwascvrtqr ) = param$modelcovariates # + EvdO
    }
    
    # stopifnot( all.equal( tcrossprod(mwascvrtqr), diag(nrow(mwascvrtqr))) );
    if(normalize){
      rez = t(orthonormalizeCovariates(mwascvrtqr, modelhasconstant));
    } else {
      if(modelhasconstant){
        rez = t(cbind(rep(1, nrow(mwascvrtqr)),mwascvrtqr));
      } else {
        rez = t(mwascvrtqr); #;
      }
    }
    
    if( !is.null( rowsubset )) colnames(rez) = param$covariates[rowsubset,1] else colnames(rez) = param$covariates[,1];
    return(rez);
  }
  
  nfolds  = param$nfolds_robustMWAS 
  
  fm      = fm.open(paste(param$dircoverageoriginal,'/','Coverage',sep=""), readonly = TRUE)
  
  # find sample in covaiate file that are in filematrix
  fm_keep = match(param$covariates[,1], rownames(fm), nomatch = 0L)
  samples = rownames(fm)[fm_keep]
  if (length(samples) / nfolds < 10) stop("Less than 10 observations in each fold")
  
  ind              = match(samples,param$covariates[,1])
  param$covariates = param$covariates[ind,]
  if (any(is.na(param$covariates[,param$modelcovariates]))) stop("Missing values in covariates not allowed")
  
  
  foldid                  = calc_foldid(param)
  fold_counts             = table( foldid )
  
  file = paste0(param$dirpca,"/eigenvalues.bmat")
  if (!file.exists(file))  ramwas4PCA(param) 
  cvrtqr    = t(.getCovariates_mod(param, modelhasconstant = TRUE)) 
  
  for (j in 1:nfolds) { # j=1
    
    new_dir      = paste0(param$dircoveragenorm_robust,"/fold_",j,"/")
    dir.create(new_dir, recursive = T, showWarnings = FALSE)
    fo           = fm.create(paste0(new_dir,'/','Coverage'), nrow=fold_counts[j], ncol = ncol(fm), size = fm$size)
    rownames(fo) = samples[ foldid == j ]
    close(fo)
  }
  
  print( "Start splitting" )
  print( param$dircoverageoriginal )
  step1  = 1024*128 # step1 = 10
  runto  = ncol(fm)
  nsteps = ceiling(runto/step1)
  for( part in seq_len(nsteps) ) { # part = 1
    cat( 'copying coverage,', part, 'of', nsteps, '\n')
    fr         = (part-1)*step1 + 1
    to         = min(part*step1, runto)
    data       =  fm[,fr:to][fm_keep,]
    
    for (j in 1:nfolds) { # j = 2
      
      new_dir = paste0(param$dircoveragenorm_robust,"/fold_",j,"/")
      fo      = fm.open(paste0(new_dir,'Coverage',sep=""), readonly = F)
      
      fo[,fr:to] = data[ foldid==j ,] 
      close(fo) 
    }
  }
  close(fm)
  
  for (j in 1:nfolds) {
    
    new_dir  = paste0(param$dircoveragenorm_robust,"/fold_",j,"/")
    
    copy_file(paste0(param$dircoverageoriginal ,'/CpG_locations.bmat'),       paste0(new_dir,'/CpG_locations.bmat'))
    copy_file(paste0(param$dircoverageoriginal ,'/CpG_locations.desc.txt'),   paste0(new_dir,'/CpG_locations.desc.txt'))
    copy_file(paste0(param$dircoverageoriginal ,'/CpG_locations.nmscol.txt'), paste0(new_dir,'/CpG_locations.nmscol.txt'))
    copy_file(paste0(param$dircoverageoriginal ,'/CpG_locations.nmsrow.txt'), paste0(new_dir,'/CpG_locations.nmsrow.txt'))
    
    copy_file(paste0(param$dircoverageoriginal ,'/CpG_chromosome_names.txt'), paste0(new_dir,'/CpG_chromosome_names.txt'))
    
  }
  
}

