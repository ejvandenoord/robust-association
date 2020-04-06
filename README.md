# robust-association

Sripts used for the paper: A method to improve the reproducibility of findings from epigenome- and transcriptome-wide association studies by Edwin JCG van den Oord, Jerry D Guintivano, and Karolina A. Aberg.

Simulations are performed with the script run_simulations.R that uses the functions in functions_run_simulations.R.

The analyses of the real data are performed using the Bioconductor package RaMWAS (PMID: 29447401 and https://www.bioconductor.org/packages/release/bioc/html/ramwas.html). RaMWAS uses data in filematrix format as input. We assume a file matrix has been created and a PCA performed on the methylation data with RaMWAS. To perform the robust MWAS, three additional functions  are required:

1.	“split_filematrix4robustMWAS.R”: This function splits the file matrix in n folds.
2.	“run_robustMWAS.R”. This functioj run the robust MWAS in RaMWAS. 
3.	“perform_meta.R”. This function reads the RaMWAS association output and perfoms the meta analysis.

Finally, the script “find_outlier.R” reanalyzes the top results by removing outliers defined as residual scores (distance between the observed value and the value predicted on the basis of the covariates) with a median absolute deviation (MAD) of >3. The scripts uses all covariates in the the file “GSE42861_pheno_data” to reanalyzes data from the top results in the file “GSE42861_topresults_data.csv”
