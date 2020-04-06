
rm(list=ls())
drive       = "G:/My Drive/"
work_dir    = paste0(drive,"SOP\\methylation\\robustMWAS\\github_scripts\\")

GSE        = "GSE42861"
n_tests    = 372348
bonf_thres = 0.05/n_tests 
mad_thres  = 3

factor2dummy = function(data) {  
  
  data_types = sapply(data, class)
  sel        = data_types=="character"
  if ( any(sel) ) {
    cols = which(sel)
    for (j in seq_along(cols)) { # j=1
      d_ = factor( data[,cols[j]] )
      if (j==1) dummies = model.matrix( ~d_ )[,-1] else 
        dummies = cbind(dummies , model.matrix( ~d_ )[,-1] )
    }
    data = cbind(data[,-cols],dummies)
  }  
  
  data  
}    

dev.off()
{ 
  
  print(GSE)
  pdf( paste0(work_dir,GSE,"\\",GSE,"_find_outlier.pdf" ) )  
  par(oma = c(4, 0, 0, 0)) 
  
  files         = list.files( paste0(work_dir,GSE,"//"), full.names = T )
  
  sel_files     = files[ grepl("_data.csv",files) ]
  sel           = grepl("pheno_data.csv",sel_files)
  
  file          = sel_files[ !sel ]
  CpG_data      = read.csv(file,row.names=1,stringsAsFactors=F)

  file          = sel_files[ sel ]
  pheno_data    = read.csv(file,row.names=1,stringsAsFactors=F)

  # put samples in in same order
  ind             = match(row.names( CpG_data ),row.names(pheno_data) )
  pheno_data      = pheno_data[ind,]
  pheno_data      = factor2dummy( pheno_data ) 
  pheno_data      = as.matrix(pheno_data)
  
  modeloutcome    = colnames(pheno_data)[1]
  modelcovariates = colnames(pheno_data)[ !(colnames(pheno_data) %in% modeloutcome) ] 
  celltypes       = colnames(pheno_data)[grepl("Norm",colnames(pheno_data))]
  use_celltypes   = celltypes[-length(celltypes)] # sums to 1 so need to drop one 
  
  file          = files[ grepl("_topresults.csv",files) ]
  sites         = read.csv(file,row.names=1,stringsAsFactors=F)
  
  ind   = match(paste0("chr",sites$chr,"_",sites$position),colnames(CpG_data))
  CpG_data =  CpG_data[,ind] 
  
 
  n_sites     = ncol(CpG_data)
  n_celltypes = length(celltypes)
  results     = data.frame( matrix(NA,n_sites,2))
  colnames(results) = c("all_P","drop_P")
  rownames(results) = rownames(  sites  )
  drop          = ""
  n_drop        = rep(NA,n_sites)

  for (i in 1:n_sites) { # i=1
    
    temp     = strsplit(rownames(sites)[i],"_") 
    target   = temp[[1]][1] 
    title    =  paste0(target,": chr",sites[i,"chr"],"_",sites[i,"position"])
    if (i==1) rowlables = title else
              rowlables = c( rowlables,title)
 
    fit = summary( lm( CpG_data[,i] ~ pheno_data[,modeloutcome] + pheno_data[,use_celltypes] + pheno_data[,modelcovariates] ) ) 
    pval = fit$coefficients[,"Pr(>|t|)"][2]  
   
    if (pval < bonf_thres) {
      
        results[i,1]                 = fit$coefficients[,"Pr(>|t|)"][2]  
        resid_fit                    = resid(fit)
        fit                          = lm(  resid_fit ~ CpG_data[,i] )  
        resid_fit                    = resid(fit)
      
        sel       =  abs(resid_fit - median(resid_fit)) / mad(resid_fit, constant=1) < mad_thres
        drop      = c(drop,rownames(pheno_data)[!sel])
        n_drop[i] = sum(!sel)

        fit  = summary( lm( CpG_data[sel,i] ~ pheno_data[sel,modeloutcome] + pheno_data[sel,use_celltypes] + pheno_data[sel,modelcovariates] ) ) 
       
        results[i,2]                 = fit$coefficients[,"Pr(>|t|)"][2]  

        xlab = paste0("drop=",sum(!sel),": P=",formatC(results[i,1],format = "e",digits = 2)," / ",formatC(results[i,2],format = "e",digits = 2))
        plot( CpG_data[,i],resid_fit,ylab="Residuals", xlab=xlab,main=title )
        abline(0, 0) 
      }
  
  } #   for (i in 1:n_sites)
  
  
  # samples that cause the outlier
  print(table(drop))
  results              = data.frame(sites,results,n_drop)
  row.names( results ) = rowlables
  summary(results)
  
  results = results[!is.na(results$n_drop),]
  write.csv(results, paste0(work_dir,GSE,"\\",GSE,"_find_outlier.csv" ),row.names=T)

  dev.off()
} 
