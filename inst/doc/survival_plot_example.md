



This report was automatically generated with the R package **knitr**
(version 1.12.3).


```r
#' Kaplan-Meier plots
#' 
#' This Rscript demonstrates joining BigQuery tables to compare survival rates conditional on somatic mutations. 
#' 
#required packages
require(bigrquery) || install.packages('bigrquery')
```

```
## [1] TRUE
```

```r
library(survival)
source('ggsurv.R')
```

```
## Warning in file(filename, "r", encoding = encoding): cannot open file
## 'ggsurv.R': No such file or directory
```

```
## Error in file(filename, "r", encoding = encoding): cannot open the connection
```

```r
# specify cloud project name
cloudProject = "isb-cgc"

# BigQuery dataset you want to work with 
dataset = "tcga_201510_alpha"

#list tables in the dataset
bigrquery::list_tables(cloudProject,dataset)
```

```
## Auto-refreshing stale OAuth token.
```

```
##  [1] "Annotations"            "Biospecimen_data"      
##  [3] "Clinical_data"          "Copy_Number_segments"  
##  [5] "DNA_Methylation_betas"  "Protein_RPPA_data"     
##  [7] "Somatic_Mutation_calls" "mRNA_BCGSC_HiSeq_RPKM" 
##  [9] "mRNA_UNC_HiSeq_RSEM"    "miRNA_expression"
```

```r
#we will use Somatic_Mutation_calls and Clinical_data for this example
somatic_mutation_table = paste("[",cloudProject,":",dataset,".Somatic_Mutation_calls",']',sep="")
clinical_table = paste("[",cloudProject,":",dataset,".Clinical_data",']',sep="")

#1. First, let's see which protein change is most common across all tumors. Let's limit the result to 20 most frequent protein changes. 
# (GROUP BY ParticipantBarcode and Protein_Change to eliminate duplicate entries)
sqlQuery = paste("SELECT Protein_Change,count(*) as N",
                 "FROM (SELECT ParticipantBarcode,Protein_Change",
                 "FROM", somatic_mutation_table,
                 "GROUP BY ParticipantBarcode,Protein_Change)",
                 "WHERE Protein_Change <> 'null'",
                 "GROUP BY Protein_Change",
                 "ORDER BY N DESC",
                 "LIMIT 20")
result = query_exec(sqlQuery,project = cloudProject)
```

```
## 123.6 megabytes processed
```

```r
#save result as a data frame
resultDF = data.frame(result)
#view results
resultDF
```

```
##    Protein_Change   N
## 1           p.M1I 792
## 2         p.V600E 469
## 3         p.R132H 436
## 4         p.E545K 273
## 5           p.M1V 267
## 6         p.K164K 248
## 7         p.Q269Q 230
## 8         p.S415S 227
## 9          p.L91L 227
## 10        p.S167S 223
## 11        p.S168S 214
## 12       p.H1047R 210
## 13        p.K416E 209
## 14          p.L9L 207
## 15          p.A2A 199
## 16         p.L13L 194
## 17        p.S332S 193
## 18         p.L17L 193
## 19       p.Q2725Q 191
## 20         p.G12D 190
```

```r
#2. Now, for the most frequent protein change, let's count the number of samples per tumor type with this protein change
#Note: The example protein changes used here are unique to their genes, but generally a protein change can be associated with
#more than one genes

top_protein_change = resultDF[1,1]

sqlQuery = paste("SELECT Study, count(*) as n ",
                 "FROM (SELECT ParticipantBarcode, Study ", 
                  "FROM ",somatic_mutation_table,
                  " WHERE Protein_Change='",top_protein_change,"' ",   
                  "GROUP BY ParticipantBarcode,Protein_Change, Study) ",
                  "GROUP BY Study ORDER BY n DESC",sep="")
result = query_exec(sqlQuery,project = cloudProject)
resultDF = data.frame(result)
resultDF
```

```
##    Study   n
## 1   SKCM 129
## 2   LUAD  96
## 3   BLCA  86
## 4   HNSC  64
## 5   LUSC  44
## 6   UCEC  39
## 7   STAD  37
## 8   BRCA  34
## 9   CESC  34
## 10  LIHC  33
## 11  THYM  24
## 12  SARC  23
## 13  KIRP  23
## 14  COAD  21
## 15  ESCA  20
## 16  TGCT  11
## 17   LGG  10
## 18  KIRC  10
## 19   GBM   8
## 20  MESO   7
## 21  DLBC   7
## 22  PRAD   6
## 23  THCA   5
## 24   ACC   5
## 25  PAAD   5
## 26    OV   4
## 27  READ   3
## 28   UCS   2
## 29  CHOL   1
## 30  PCPG   1
```

```r
#let's look at only top 10 studies for survival comparison.
top_ten_studies = resultDF$Study[1:10]


#3. let's see how survival trends look across different tumor types with the most common protein change
sqlQuery = paste("SELECT clin.ParticipantBarcode,vital_status, clin.Study, days_to_last_known_alive ",
                 "FROM ",clinical_table," as clin JOIN (SELECT ParticipantBarcode, Study ", 
                 "FROM ",somatic_mutation_table,
                 " WHERE Protein_Change='",top_protein_change,"' ",   
                 "GROUP BY ParticipantBarcode,Protein_Change, Study) as mut ", 
                 "ON clin.ParticipantBarcode=mut.ParticipantBarcode",sep="")
result = query_exec(sqlQuery,project = cloudProject)
resultDF = data.frame(result)
head(resultDF)
```

```
##   clin_ParticipantBarcode vital_status clin_Study days_to_last_known_alive
## 1            TCGA-DD-A4NR         Dead       LIHC                        9
## 2            TCGA-44-5644        Alive       LUAD                      863
## 3            TCGA-99-8028        Alive       LUAD                     1118
## 4            TCGA-FD-A43X        Alive       BLCA                      110
## 5            TCGA-DK-AA6L         Dead       BLCA                     1163
## 6            TCGA-DK-AA6U        Alive       BLCA                      578
```

```r
#convert vital status to numeric type
resultDF$vital_status[resultDF$vital_status=="Alive"] = 0
resultDF$vital_status[resultDF$vital_status=="Dead"] = 1
class(resultDF$vital_status) = "numeric"

#convert Study type to factor
resultDF$clin_Study = factor(resultDF$clin_Study)

#normalize days_to_last_known_alive within each tumor type to mitigate inherent survival differences across tumor types
for (tumor_type in levels(resultDF$clin_Study))
{
  survival_this_tumor = resultDF$days_to_last_known_alive[resultDF$clin_Study==tumor_type]
  
  min_last_known_alive = min(survival_this_tumor,na.rm = TRUE)
  max_last_known_alive = max(survival_this_tumor,na.rm = TRUE)
  
  resultDF$days_to_last_known_alive[resultDF$clin_Study==tumor_type] = (resultDF$days_to_last_known_alive[resultDF$clin_Study==tumor_type] - min_last_known_alive)/(max_last_known_alive - min_last_known_alive)
  
}

#fit survival curve to data
tcga.surv = survfit(Surv(resultDF$days_to_last_known_alive,resultDF$vital_status)~resultDF$clin_Study,data = resultDF)
#plot survival curve
ggsurv(tcga.surv,main = paste("Protein Change = ",top_protein_change))
```

<img src="figure/survival-plot-example-Rmdauto-report-1.png" title="plot of chunk auto-report" alt="plot of chunk auto-report" style="display: block; margin: auto;" />

```r
#4. Let's only plot top_ten_studies (if there are more than 10 curves in the previous plot)
resultDF = resultDF[resultDF$clin_Study %in% top_ten_studies,]
#fit survival curve to data
tcga.surv = survfit(Surv(resultDF$days_to_last_known_alive,resultDF$vital_status)~resultDF$clin_Study,data = resultDF)
#plot survival curve
ggsurv(tcga.surv,main = paste("Protein Change = ",top_protein_change))
```

<img src="figure/survival-plot-example-Rmdauto-report-2.png" title="plot of chunk auto-report" alt="plot of chunk auto-report" style="display: block; margin: auto;" />

```r
#5. let's add a baseline survival curve to this plot. For this, we'll pull all the samples from top_ten_studies that do not have 
#the top_protein_change
#Make the following simple changes to the query above
#* LEFT JOIN clinical_table to somatic_mutation_table, so that the query result contains all the samples from the clinical data
# regardles of whether they have a mutation or not.
#* Add Protein_Change to list of output columns. It's value will be null for non-matching samples.
#* Add a WHERE clause to get samples only from top_ten_studies

sqlQuery = paste("SELECT clin.ParticipantBarcode,vital_status, clin.Study, days_to_last_known_alive, Protein_Change ",
                 "FROM ",clinical_table," as clin LEFT JOIN (SELECT ParticipantBarcode, Protein_Change, Study ", 
                 "FROM ",somatic_mutation_table,
                 " WHERE Protein_Change='",top_protein_change,"' ",   
                 "GROUP BY ParticipantBarcode,Protein_Change, Study) as mut ", 
                 "ON clin.ParticipantBarcode=mut.ParticipantBarcode ",
                 "WHERE clin.Study in (",paste(shQuote(top_ten_studies),collapse = ","),")",sep="")
result = query_exec(sqlQuery,project = cloudProject)
resultDF = data.frame(result)
head(resultDF)
```

```
##   clin_ParticipantBarcode vital_status clin_Study days_to_last_known_alive
## 1            TCGA-FP-A9TM        Alive       STAD                      189
## 2            TCGA-DD-A4NA        Alive       LIHC                     1008
## 3            TCGA-DD-A4NB        Alive       LIHC                      989
## 4            TCGA-DD-A4ND        Alive       LIHC                     2746
## 5            TCGA-DD-A4NE         Dead       LIHC                      660
## 6            TCGA-DD-A4NF        Alive       LIHC                      942
##   Protein_Change
## 1           <NA>
## 2           <NA>
## 3           <NA>
## 4           <NA>
## 5           <NA>
## 6           <NA>
```

```r
#convert vital status to numeric type
resultDF$vital_status[resultDF$vital_status=="Alive"] = 0
resultDF$vital_status[resultDF$vital_status=="Dead"] = 1
class(resultDF$vital_status) = "numeric"

#make all controls to be in 'All_Controls' study type
resultDF$clin_Study[is.na(resultDF$Protein_Change)] = "All_Controls"

#convert Study type to factor
resultDF$clin_Study = factor(resultDF$clin_Study)

#normalize days_to_last_known_alive within each tumor type to mitigate inherent survival differences across tumor types
for (tumor_type in unique(resultDF$clin_Study))
{
  survival_this_tumor = resultDF$days_to_last_known_alive[resultDF$clin_Study==tumor_type]
  
  min_last_known_alive = min(survival_this_tumor,na.rm = TRUE)
  max_last_known_alive = max(survival_this_tumor,na.rm = TRUE)
  
  resultDF$days_to_last_known_alive[resultDF$clin_Study==tumor_type] = (resultDF$days_to_last_known_alive[resultDF$clin_Study==tumor_type] - min_last_known_alive)/(max_last_known_alive - min_last_known_alive)
  
}

#fit survival curve to data
tcga.surv = survfit(Surv(resultDF$days_to_last_known_alive,resultDF$vital_status)~resultDF$clin_Study,data = resultDF)
#plot survival curve
ggsurv(tcga.surv,main = paste("Protein Change = ",top_protein_change))
```

<img src="figure/survival-plot-example-Rmdauto-report-3.png" title="plot of chunk auto-report" alt="plot of chunk auto-report" style="display: block; margin: auto;" />

```r
#6. Let's look at the survival curve for the tumor type where the top_protein_change is most common.
#For this, we need one simple change to the previous sql query.
#* edit the WHERE claus to return only top_ten_studies[1] records
# Note: Also try top_protein_change='p.V600E' and top_study='THCA' 

sqlQuery = paste("SELECT clin.ParticipantBarcode,vital_status, clin.Study, days_to_last_known_alive, Protein_Change ",
                 "FROM ",clinical_table," as clin LEFT JOIN (SELECT ParticipantBarcode, Protein_Change, Study ", 
                 "FROM ",somatic_mutation_table,
                 " WHERE Protein_Change='",top_protein_change,"'",   
                 " GROUP BY ParticipantBarcode,Protein_Change, Study) as mut", 
                 " ON clin.ParticipantBarcode=mut.ParticipantBarcode",
                 " WHERE clin.Study='",top_ten_studies[1],"'",sep="")

result = query_exec(sqlQuery,project = cloudProject)
resultDF = data.frame(result)
head(resultDF)
```

```
##   clin_ParticipantBarcode vital_status clin_Study days_to_last_known_alive
## 1            TCGA-ER-A199         Dead       SKCM                      279
## 2            TCGA-ER-A196        Alive       SKCM                     1785
## 3            TCGA-ER-A194         Dead       SKCM                     1354
## 4            TCGA-WE-A8K1        Alive       SKCM                     1492
## 5            TCGA-WE-A8K5         Dead       SKCM                     1860
## 6            TCGA-WE-A8K4        Alive       SKCM                      614
##   Protein_Change
## 1           <NA>
## 2           <NA>
## 3           <NA>
## 4           <NA>
## 5           <NA>
## 6           <NA>
```

```r
#convert vital status to numeric type
resultDF$vital_status[resultDF$vital_status=="Alive"] = 0
resultDF$vital_status[resultDF$vital_status=="Dead"] = 1
class(resultDF$vital_status) = "numeric"

#Convert Protein_Change column to a TRUE/FALSE column so we can compare survival trends between cases and controls
resultDF$Protein_Change = as.logical(!is.na(resultDF$Protein_Change))

#fit survival curve to data
tcga.surv = survfit(Surv(resultDF$days_to_last_known_alive,resultDF$vital_status)~resultDF$Protein_Change,data = resultDF)
#plot survival curve
ggsurv(tcga.surv,main = paste("Protein Change = ",top_protein_change,", Study = ",top_ten_studies[1]))
```

<img src="figure/survival-plot-example-Rmdauto-report-4.png" title="plot of chunk auto-report" alt="plot of chunk auto-report" style="display: block; margin: auto;" />

The R session information (including the OS info, R version and all
packages used):


```r
sessionInfo()
```

```
## R version 3.2.4 (2016-03-10)
## Platform: x86_64-apple-darwin13.4.0 (64-bit)
## Running under: OS X 10.11.4 (El Capitan)
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] knitr_1.12.3    bigrquery_0.2.0 ggplot2_2.1.0   survival_2.38-3
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_0.12.4      magrittr_1.5     splines_3.2.4    munsell_0.4.3   
##  [5] colorspace_1.2-6 R6_2.1.2         highr_0.5.1      stringr_1.0.0   
##  [9] httr_1.1.0       plyr_1.8.3       dplyr_0.4.3      tools_3.2.4     
## [13] parallel_3.2.4   grid_3.2.4       gtable_0.2.0     DBI_0.3.1       
## [17] htmltools_0.3.5  openssl_0.9.2    yaml_2.1.13      assertthat_0.1  
## [21] digest_0.6.9     formatR_1.3      curl_0.9.7       evaluate_0.8.3  
## [25] rmarkdown_0.9.5  labeling_0.3     stringi_1.0-1    scales_0.4.0    
## [29] jsonlite_0.9.19
```

```r
Sys.time()
```

```
## [1] "2016-05-03 12:48:18 PDT"
```

