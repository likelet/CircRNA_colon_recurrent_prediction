# CircRNA_colon_recurrent_prediction

# Introducntion  
This repo contains the analysis and plot code for constructing cirRNA based signature to predict recurrence of stage II or III colon cancers. For privacy requirement, all patients's information were masked an provided a random number as their identifier in this analysis. 


# Citations

# Authors  

Qi Zhao and Zi-xian Wang from Sun Yat-sen University 


# Sessions for running the command 

```R
> sessionInfo()
R version 3.3.3 (2017-03-06)
Platform: x86_64-apple-darwin13.4.0 (64-bit)
Running under: macOS Sierra 10.12.6

locale:
[1] zh_CN.UTF-8/zh_CN.UTF-8/zh_CN.UTF-8/C/zh_CN.UTF-8/zh_CN.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

loaded via a namespace (and not attached):
 [1] Rcpp_0.12.14        mvtnorm_1.0-6       lattice_0.20-34     tidyr_0.7.2         zoo_1.8-0          
 [6] assertthat_0.2.0    digest_0.6.13       psych_1.7.5         R6_2.2.2            survminer_0.4.3    
[11] plyr_1.8.4          backports_1.1.0     acepack_1.4.1       MatrixModels_0.4-1  ggplot2_3.1.0.9000 
[16] rlang_0.3.0.9000    multcomp_1.4-6      lazyeval_0.2.1      data.table_1.10.4   SparseM_1.77       
[21] rpart_4.1-10        Matrix_1.2-8        checkmate_1.8.2     splines_3.3.3       stringr_1.2.0      
[26] foreign_0.8-67      htmlwidgets_0.8     munsell_0.4.3       broom_0.4.2         pkgconfig_2.0.1    
[31] base64enc_0.1-3     mnormt_1.5-5        htmltools_0.3.6     nnet_7.3-12         tibble_1.3.4       
[36] gridExtra_2.3       km.ci_0.5-2         htmlTable_1.9       codetools_0.2-15    Hmisc_4.0-3        
[41] rms_5.1-1           dplyr_0.7.4         ggpubr_0.1.8.999    MASS_7.3-45         grid_3.3.3         
[46] polspline_1.1.12    nlme_3.1-131        xtable_1.8-2        gtable_0.2.0        magrittr_1.5       
[51] KMsurv_0.1-5        pROC_1.9.1          scales_0.5.0.9000   stringi_1.1.6       reshape2_1.4.3     
[56] bindrcpp_0.2        latticeExtra_0.6-28 survMisc_0.5.4      sandwich_2.3-4      TH.data_1.0-8      
[61] Formula_1.2-1       RColorBrewer_1.1-2  tools_3.3.3         cmprsk_2.2-7        glue_1.2.0         
[66] purrr_0.2.4         yaml_2.1.14         parallel_3.3.3      survival_2.41-3     colorspace_1.3-2   
[71] cluster_2.0.5       knitr_1.20          bindr_0.1           quantreg_5.33      
```