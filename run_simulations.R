
rm(list=ls())

drive    = "G:\\My Drive\\"
work_dir = paste0(drive,"SOP\\methylation\\robustMWAS\\github_scripts\\")
source(paste0(work_dir,"functions_run_simulations.R")) 
library(MASS)

sim_type="H0_split_half_replication"

set.seed(1)

n_sample  = 250
preplic   = c(1e-3)
preplicH0 = 0.05 


# Type I error
n_sims          = 1e2
nfolds_vector   = c(2,10,50)    
r_vector        = c(0)
poutlier_vector = c(0.005,0.02) 
toutlier_vector = c("dep","indepboth") 
min_error2      = 3
max_error2      = 5
pdisc_vector    = 1.5 # values larger than 1 will not select
directional     = F   # if pdisc < 1 should directional replcation tests be used? 

### power
#n_sims            = 1e4
#nfolds_vector     = c(2,10,50)    
#r_vector          = c(0.3)
#poutlier_vector   = c(0,0.02) 
#toutlier_vector   = c("indepboth") 
#pdisc_vector      = 0.05 # values larger than 1 will not select
#directional       = F   # if pdisc < 1 should directional replcation tests be used? 
#min_error2=3
#max_error2=5
#pdisc_vector   = 1.5 # values larger than 1 will not select

### split half
#n_sims          = 1e5
#nfolds_vector   = c(2)
#r_vector        = c(0)
#poutlier_vector = c(0.005,0.02)
#toutlier_vector = c("dep","indepboth") 
#min_error2      = 3
#max_error2      = 5
#pdisc_vector   = 0.05 # values larger than 1 will not select
#directional   = F   # if pdisc < 1 should directional replcation tests be used? 


# all sims will have same random components (pred_data,error,foldid - same seed used)) to increase comprability 
start_seed = sample(1:(.Machine$integer.max-(n_sims+1)), 1)
print(start_seed)
dir.create(paste0(work_dir,"\\results\\"), showWarnings=F, recursive=T)

n_conditions    = length(r_vector)*length(nfolds_vector)*length(toutlier_vector)*length(poutlier_vector*length(pdisc_vector))
conditions      = matrix(NA,n_conditions ,5)
analysis_labels = rep(NA,n_conditions)
colnames(conditions) = c("r","nfolds","t_outlier","p_outlier","pdisc")
x=0
for (i in seq_along(r_vector) ) 
  for (j in seq_along(nfolds_vector) )
    for (k in seq_along(toutlier_vector) )
      for (l in seq_along(poutlier_vector) ) 
        for (m in seq_along(pdisc_vector) ) { # i=j=k=l=m=1
          
        x=x+1
        conditions[x,1]    = r        = r_vector[i] 
        conditions[x,2]    = nfolds   = nfolds_vector[j]
        conditions[x,3]    = toutlier = toutlier_vector[k]
        conditions[x,4]    = poutlier = poutlier_vector[l]
        conditions[x,5]    = pdisc    = pdisc_vector[m] 

        analysis_labels[x] = paste0("r",r,"_nfolds",nfolds,
                                    "_poutl.",toutlier,"_poutl",poutlier,
                                    "_pdisc",pdisc )
}
rownames(conditions) = analysis_labels

ncol=8
results   = matrix(NA,n_conditions ,ncol)
colnames(results) = c("mean_cor_train","mean_cor","mean_p","mean_t",
                      "nfolds<p.repl","p_fold1","p_meta","p_all")

x=0
for (i in seq_along(r_vector) ) 
  for (j in seq_along(nfolds_vector) )
    for (k in seq_along(toutlier_vector) )
      for (l in seq_along(poutlier_vector) ) 
        for (m in seq_along(pdisc_vector) ) { # i=j=k=l=m=1
          
          x=x+1
          
  #      if ( r_vector[i] == 0 ) n_sims = n_sims * 50
        if ( r_vector[i] == 0 ) preplic = preplicH0
        
        print( analysis_labels[x] )   
        results = matrix(NA,n_sims,ncol)
        for (iter in 1:n_sims )  { # iter=1
          if(iter %% 1000==0) cat(paste0("Simulation: ", iter, "\n"))
          set.seed(iter+start_seed)
          
          pred_data  = rnorm(n_sample)
          error      = rnorm(n_sample)
          r          = r_vector[i]
          outcome    = r*pred_data + sqrt(1-r^2)*error
          
          n_outlier = round( poutlier_vector[l] * n_sample )
          if (n_outlier>0) {
            error2               = runif(n_outlier,min=min_error2,max=max_error2)  
 
            if (toutlier_vector[k] == "dep") {
              pred_data[1:n_outlier] = abs(pred_data[1:n_outlier]) + error2
              outcome[1:n_outlier]   = abs(outcome[1:n_outlier] + error2) }

            if (toutlier_vector[k] == "indepboth") {
               dir = sign( outcome[1:n_outlier] )
               outcome[1:n_outlier] = outcome[1:n_outlier] + dir*error2
               dir = sign( pred_data[1+n_outlier:(2*n_outlier-1)] )
               pred_data[1+n_outlier:(2*n_outlier-1)] = pred_data[1+n_outlier:(2*n_outlier-1)] + dir*error2
             } 
          }         

          results[iter,] = run_lm_sim(n_sample,pred_data,outcome,nfolds_vector[j],pdisc_vector[m]) 

     }
  
      temp =  calc_summary_results(results,preplic) 
      if (x==1) summary_results =  temp else
                summary_results = rbind(summary_results,temp)
      print (temp )
      write.csv(results,paste0(work_dir,"\\results\\",analysis_labels[x],".csv"),row.names=F,quote=F)
  
}
rownames(summary_results) = analysis_labels
temp = data.frame(cbind(conditions,summary_results) )
temp
write.csv(temp,paste0(work_dir,sim_type,"_summary_results.csv"),row.names=T,quote=F)
