


erf    = function(x) 2 * pnorm(2 * x/ sqrt(2)) - 1
erfinv = function(x) qnorm( (x+1)/2 ) / sqrt(2)
pcomb  = function(p) (1-erf(sum(sqrt(2) * erfinv(1-2*p))/sqrt(2*length(p))))/2

tt2z = function(tt,df) { lpv = pt(-abs(tt), df, log.p = TRUE); z = -sign(tt) * qnorm(lpv, log.p = TRUE); return(z); };

calc_summary_results = function(results,preplic) {
  
  summary_results      = rep(NA,5+3*length(preplic))
  names(summary_results)[1:5] = c("mean_cor_train","mean_cor","mean_p","mean_t","nfolds<p.repl")
  summary_results[1:5] = colMeans(results[,c(1:5)],na.rm=T)

  z = 5
  for (i in seq_along(preplic)) { # i=1
    z = z+1
    sel = !is.na( results[,6] )
    summary_results[z] = sum(results[sel,6]< preplic[i]) / sum(sel)
    names(summary_results)[z] = paste0("p_fold1_prepl_",preplic[i])
    
    z = z+1
    sel = !is.na( results[,7] ) & !is.na( results[,6] )
    summary_results[z]   = sum(results[sel,7]< preplic[i]) / sum(sel)
    names(summary_results)[z] = paste0("p_meta_prepl_",preplic[i])
    
    z = z+1
    sel = !is.na( results[,8] ) & !is.na( results[,6] )
    summary_results[z]   = sum(results[sel,8]< preplic[i]) / sum(sel)
    names(summary_results)[z] = paste0("p_all_prepl_",preplic[i])
    
    }

   summary_results
} 



calc_folds = function(nfolds,n_sample) {
  
  set.seed(1) # thos ensures sample fols are used for each effect size
  fold_size = rep(floor(n_sample/nfolds),nfolds)
  rest      = n_sample-sum(fold_size)
  if (rest>0) fold_size[1:rest] = fold_size[1:rest] + 1
  foldid    = sample(rep(1:nfolds,fold_size),n_sample)  
  foldid
}


predict_stats = function(cv.lm.result,out,data) { 
  
  # out=test_out;data=test_data

  if ( directional & pdisc < 1 ) {
     if (sign(cv.lm.result$coefficients[2]) == 1) temp=cor.test(out,data,alternative="greater") else
                                                  temp=cor.test(out,data,alternative="less")
  } else temp=cor.test(out,data)
  
  results = unlist( temp[c("estimate","p.value","statistic")] )

  results
}

run_lm_sim = function(n_sample,pred_data,outcome,nfolds,pdisc) {
  
  # nfolds=nfolds_vector[j]; pdisc=pdisc_vector[m]
  
  # nfolds=10; pdisc=0.01

  foldid  = calc_folds(nfolds,n_sample)

  train_cor = rep(NA,nfolds)
  test_stats = matrix(NA,nfolds,3)
  colnames(test_stats) = c("cor","pval","tstat")
  for (ii in 1:nfolds) { # ii=1

    train_data = pred_data[ which(foldid != ii)] 
    train_out  = outcome[which(foldid != ii)]
    test_data  = pred_data[ which(foldid == ii)]
    test_out   = outcome[which(foldid == ii)]

    train_cor[ii]   = cor( train_data, train_out )
    cv.lm.result    = lm(train_out~train_data)	
    if (summary(cv.lm.result)$coefficients[,"Pr(>|t|)"][2] < pdisc) test_stats[ii,]  = predict_stats(cv.lm.result,test_out,test_data)

   } #   for (i in 1:nfolds)
  
#  fisher_p    = pchisq( -2*sum(log(cor[,"p.value"])), length(cor[,"parameter.df"])*2, lower.tail=FALSE)

  sel = !is.na(test_stats[,1])
  if (any(sel)) {
    if ( directional & pdisc < 1 ) {
         metap   = pcomb(test_stats[sel,"pval"]) } else {
         df    = (table(foldid) - 2)[sel]
         zcomb = sum( tt2z(test_stats[,"tstat"],df) ) / sqrt( length(df) ) 
         metap = pnorm(abs(zcomb),lower.tail=FALSE)*2
         }
  } else metap=p_metap=NA 
  
 
  result      = c(mean(train_cor,na.rm=T),
                  colMeans(test_stats[,c(1,2,3)],na.rm=T),sum(sel),
                  test_stats[1,"pval"],
                  metap,
                  unlist( cor.test(pred_data,outcome)["p.value"] )) 

  result
}

