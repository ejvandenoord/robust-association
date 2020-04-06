

perform_meta = function(metalocations,dircoveragenorm) {
  
  getMWASinfo = function(path,dircoveragenorm,df_correction) {
    
    rez = list();
    rez$dfs = as.double(readLines(paste0(path,'/DegreesOfFreedom.txt')));
    
    message(  "rez$df 1=",  rez$dfs[1] )
    message(  "rez$df 2=",  rez$dfs[2] )
    
    fm = fm.open(paste0(path,'/Stats_and_pvalues.bmat'));
    
    # colnames(fm)
    if(is.infinite(rez$dfs[2])) {
      tt2z = identity;
      message("D.F. = ",rez$dfs[2]);
    } else
      if(rez$dfs[1] == 1) {
        tt2z = function(tt, dfFull) { lpv = pt(-abs(tt), dfFull, log.p = TRUE); z = -sign(tt) * qnorm(lpv, log.p = TRUE); return(z); };
        message("D.F. = ",rez$dfs[2]);
      } else {
        stop("F-test (categorical outcome) meta not supported yet");
      }
    
    rez$c = fm[,1]
   
    dfFull = rez$dfs[2] 
    rez$z  = tt2z(fm[,2],dfFull);
    rez$p  = fm[,3];
    close(fm);
    
    rez$loc      = fm.load(paste0(dircoveragenorm,'/CpG_locations'));
    rez$chrnames = readLines(paste0(dircoveragenorm,'/CpG_chromosome_names.txt'));
    
     stopifnot( length(rez$z) == NROW(rez$loc) );
 
    return(rez);
  }
  
  param=list()
  param$metalocations=metalocations
  
  nsets = length(param$metalocations);
  
  # Aggregate the z-scores and D.F's
  {
    # accumulating variables
    clist = vector('list', nsets);
    zlist = vector('list', nsets);
    plist = vector('list', nsets);
    names(clist) = paste0('cor.dataset.',seq_along(clist));
    names(zlist) = paste0('z.dataset.',seq_along(zlist));
    names(plist) = paste0('pv.dataset.',seq_along(plist));
    
    if( is.null(param$weights) ) {
      weights2 = integer(nsets);
    } else {
      weights2 = param$weights;
    }
    for( j in seq_len(nsets) ) { # j=1
      dir = param$metalocations[j];
      message('Processing ', dir)
      rez = getMWASinfo(path = dir,dircoveragenorm[j]);
      
      if( is.null(param$weights) )
        weights2[j] = rez$dfs[2]
      
      
      if( j == 1 ) {
        chrset = rez$chrnames;
        locations = rez$loc;
        stopifnot(all( locations[,2] < 1e9L ));
        
        zlist[[j]] = rez$z;
        plist[[j]] = rez$p;
        clist[[j]] = rez$c;
        cumz = rez$z * sqrt(weights2[j]);
        cumc = rez$c * weights2[j];
      } else {
        chrset = c(chrset, rez$chrnames[!(rez$chrnames %in% chrset)]);
        if(any(rez$chrnames != chrset[seq_along(rez$chrnames)]))
          stop("mismatching chromosome names, ask Andrey to fix");
        
        mch = match(locations %*% c(1e9,1) , rez$loc %*% c(1e9,1), nomatch = 0L);
        keep = which(mch > 0L);
        zlist = lapply( zlist, `[`, keep)
        plist = lapply( plist, `[`, keep)
        clist = lapply( clist, `[`, keep)
        locations = locations[keep,];
        
        zlist[[j]] = rez$z[mch]
        plist[[j]] = rez$p[mch]
        clist[[j]] = rez$c[mch]
        cumz = cumz[keep] + rez$z[mch] * sqrt(weights2[j]);
        cumc = cumc[keep] + rez$c[mch] * weights2[j];
      }
    }
  } # zlist, plist, cumz, weights2, chrset
  
   ### Calculate finalz, pv
  {
    message('Calculating meta p-values')
    cumweight2 = sum(weights2);
    finalz = cumz / sqrt(cumweight2);
    finalc = cumc / cumweight2;
    pv = pnorm(abs(finalz), lower.tail = FALSE)*2;
    message('Calculating meta q-values')
    qv = ramwas:::pvalue2qvalue(pv)
    # hist(pv, 100)
  }
  
  
  cov_loc        = as.vector(locations[,1]*1e9 + locations[,2])
  
  result = cbind(cov_loc,locations[,1],locations[,2],finalc,finalz,pv,qv)
  colnames(result) = c("loc","chr","pos","cor","Zstat","pval","qval")
  
  result
  
} # perfrom meta

