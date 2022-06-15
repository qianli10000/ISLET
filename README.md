# ISLET: Individual-Specific ceLl typE referencing Tool 


-------------------
<img align="left" src="vignettes/islet_hex_2.png" width="156" height="180"> We developed a method called `ISLET` (Individual-Specific ceLl typE referencing Tool), to deconvolute mixture samples and obtain the individual-specific and cell-type-specific reference panels. `ISLET` can leverage on multiple observations or temporal measurements of the same subject. `ISLET` adopted a more reasonable assumption that repeated samples from the same subject would share the same reference panel. This unknown panel, treated as missing values, are estimated by an iterative Expectation--Maximization (EM) algorithm in mixed-effect regression framework, when combining all samples from all subjects together. This is the first statistical framework to estimate the subject-level cell-type-specific reference panel, for repeated measures. Our modeling can effectively borrow information across samples within the same subject.


#### `ISLET` is currently being sumitted to Bioconductor.



# 1. Introduction
In clinical samples, the observed bulk sequencing/microarray data are often a mixture of different cell types. Because each unique cell type has its own gene expression profile, the real sequencing/microarray data are the weighted average of signals from multiple pure cell types. In high-throughput data analysis, the mixing proportions will confound with the primary phenotype-of-interest, if not properly accounted for. Over the past several years, researchers have gained substantial interests in using computational methods to deconvolute cell compositions. Under the assumption of a commonly shared feature-by-cell-type reference panel across all samples, deconvolution methods were developed. However, this assumption may not hold. For example, when repeated samples are measured for each subject, assuming a shared reference panel across different time points for each subject is a preferred choice over assuming a shared one across all the samples.

Here, we developed a method called `ISLET` (Individual-Specific ceLl typE referencing Tool), to deconvolute mixture samples and obtain the individual-specific and cell-type-specific reference panels. `ISLET` can leverage on multiple observations or temporal measurements of the same subject. `ISLET` adopted a more reasonable assumption that repeated samples from the same subject would share the same reference panel. This unknown panel, treated as missing values, are estimated by an iterative Expectation--Maximization (EM) algorithm in mixed-effect regression framework, when combining all samples from all subjects together. This is the first statistical framework to estimate the subject-level cell-type-specific reference panel, for repeated measures. Our modeling can effectively borrow information across samples within the same subject.

![Schematic overview of ISLET workflow.](vignettes/fig1.png)

`ISLET` depends on the following packages:

-   `r Biocpkg("BiocParallel")`, for parallel computing implementation,
-   `r CRANpkg("Matrix")`, for large matrices operations in *R*.

# Preparing ISLET input files for deconvolution:
ISLET needs to have the feature by sample matrix for observed values, sample by cell type matrix for cell type proportions, and data frame showing the relationship between samples to subjects. The requirement is the same for both case group and control group, resulting to a total of 6 elements as the input. An example dataset `GE600` is included to show what is required:

ISLET needs to have the feature by sample matrix for observed values, sample by cell type matrix for cell type proportions, and data frame showing the relationship between samples to subjects. The requirement is the same for both case group and control group, resulting to a total of 6 elements as the input. An example dataset `GE600` is included to show what is required:

**Step 1**: Load in example data.

```
library(ISLET)
data(GE600)
ls()
```

`case_obs` is a gene expression value data frame of 1,000 genes by 300 cancer sample, for the cancer group. `ctrl_obs` is a gene expression value data frame of 1,000 genes by 300 control sample, for the control group. 

`case_prop` is a cell proportion matrix of 300 cancer samples by 6 cell types, for the cancer group. `ctrl_prop` is a cell proportion matrix of 300 control samples by 6 cell types, for the control group. 

`case_subject` is a two-column data frame showing the relationship between the 300 cancer samples and their 100 subjects, for the cancer group. The first column is the sample_id. The second column is the subject_id. Each subject has 3 repeated/temporal measurements. `ctrl_subject` is a two-column data frame showing the relationship between the 300 control samples and their 100 subjects, for the control group. The format is the same with `case_subject`.

# Data preparation
This should always be the first step to before using ISLET for deconvolution or testing. This initial step for will prepare your data input ready for downstream deconvolution (function `islet.solve`) and/or differentially expressed gene testing (function `islet.test`).

**Step 2**: Data preparation for downstream ISLET analysis. 

```
study123input = dataprep(case_obs = case_obs, ctrl_obs = ctrl_obs,
         case_prop = case_prop, ctrl_prop = ctrl_prop,
         case_subject = case_subject, ctrl_subject = ctrl_subject)
```

[**Attention**] Here we have strict requirement for the input data. The features (i.e. rows) in `case_obs` and `ctrl_obs` must match each other. The samples (i.e. rows) in `case_prop`/`ctrl_prop` must match samples (i.e. columns) in `case_obs`/`ctrl_obs`, respectively. The subject IDs across case and control groups must be unique, and sorted in advance. 

# Using ISLET for deconvolution

**Step 3**: With the curated data `study123input` from the previous step, now we can use `ISLET` to conduct deconvolution and obtain the individual-specific and cell-type-specific reference panels. This process can be achieved with one line of code:

```
#Use ISLET for deconvolution
res.sol = islet.solve(input = study123input)
```

The `res.sol` is the deconvolution result output is a list of length G, where G is the number of features (i.e. genes). In the example dataset, G=1000. For each of the G elements, it is a list of length 7. Below is the description of the 7 components: 

`case.m` A vector of length K, where K is the number of cell types. In the example dataset, K=6. It contains the mean cell type profiles (i.e. mean cell type gene expression values for K cell types), for case group, on a specific feature. 

`ctrl.m` A vector of length K, where K is the number of cell types. It contains the mean cell type profiles (i.e. mean cell type gene expression values for K cell types), for control group, on a specific feature.

`case.indv` A subject by cell type matrix containing all the feature values (i.e. gene expression values), for case group. It is one of the main products the individual-specific and cell-type-specific solve algorithm. Values are the summation of fixed effect (`case.m`) and individual random effect (not shown explicitly).

`ctrl.indv` A subject by cell type matrix containing all the feature values (i.e. gene expression values), for control group. It is one of the main products the individual-specific and cell-type-specific solve algorithm. Values are the summation of fixed effect (`ctrl.m`) and individual random effect (not shown explicitly).

`var.k` A vector of length K, where K is the number of cell types. It contains the estimated feature value variance for the K cell types, on a specific feature.

`var.0` A scalar. The estimated residual variance, on a specific feature.

`LLK` A scalar. The log-likelihood from the current model. It can be useful for testing purpose such as Likelihood Ratio Test.



# Using ISLET for cell-type-specific differential expression (csDE) genes identification
**Step 4**: Also, with the curated data `study123input` from the previous **Step 2**, now we can test for csDE genes. Notice that Step 3 is optional, if calling csDE genes is the only target of your investigation. This step is done by the following line of code:

```
#Test for csDE genes
res.test = islet.test(input = study123input)
```

The result `res.test` is a matrix of p-values, in the dimension of feature by cell type. Each element is the LRT p-value, by contrasting case group and control group, for one feature in one cell type.

# Session info 

```
## R version 4.2.0 (2022-04-22)
## Platform: x86_64-apple-darwin17.0 (64-bit)
## Running under: macOS Big Sur/Monterey 10.16
## 
## Matrix products: default
## BLAS:   /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRblas.0.dylib
## LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] parallel  stats     graphics  grDevices utils     datasets  methods  
## [8] base     
## 
## other attached packages:
## [1] ISLET_0.99.1        BiocParallel_1.31.8 Matrix_1.4-1       
## [4] BiocStyle_2.25.0   
## 
## loaded via a namespace (and not attached):
##  [1] rstudioapi_0.13     knitr_1.39          magrittr_2.0.3     
##  [4] lattice_0.20-45     R6_2.5.1            rlang_1.0.2        
##  [7] fastmap_1.1.0       stringr_1.4.0       tools_4.2.0        
## [10] grid_4.2.0          xfun_0.31           cli_3.3.0          
## [13] jquerylib_0.1.4     htmltools_0.5.2     yaml_2.3.5         
## [16] digest_0.6.29       bookdown_0.27       BiocManager_1.30.18
## [19] codetools_0.2-18    sass_0.4.1          evaluate_0.15      
## [22] mime_0.12           rmarkdown_2.14      stringi_1.7.6      
## [25] compiler_4.2.0      bslib_0.3.1         jsonlite_1.8.0
```
