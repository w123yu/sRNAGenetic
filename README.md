# sRNAGenetic

The goal of sRNAGenetic is to analysis the expression changes of sRNA after plant polyploidization. The most important function of the R package sRNAGenetic is the genetic effects analysis of miRNA after plant polyploidization via two methods, and at the same time, it provides various forms of graph related to data characteristics and expression analysis. In terms of two classification methods, one is the calculation of the additive (a) and dominant (d), the other is the evaluation of ELD by comparing the total expression of the miRNA in allotetraploid with the expression level in the parent species.

## Installation

You can install the development version of sRNAGenetic like so:

``` r
BiocManager::install("sRNAGenetic")
```

Loading the R package: sRNAGenetic

``` r
library(sRNAGenetic)
```
## Data Statistics
### Length distribution plot of sRNA

```r
## Recommended to use the "data.table" package for reading data quickly.
## The first column of input data must be their sequences.
P1_sRNA <- srnapredata(srnaseq_dataframe = P1_sRNA_seq, group = "P1")
P2_sRNA <- srnapredata(srnaseq_dataframe = P2_sRNA_seq, group = "P2")
F1_sRNA <- srnapredata(srnaseq_dataframe = F1_sRNA_seq, group = "F1")

## Intergrate all sRNA data
sRNA_data <- rbind(P1_sRNA,P2_sRNA,F1_sRNA)
```
Note: The first column of the input file must be the sRNA sequence. About the output result, the first column is the length of sRNA, the second column of output file is the frequency of occurrence of the length, the third column is the file name representing the grouping information. 

```r
## Length distribution plot
lenplot(file_dataframe = sRNA_data)
```

### Base preference for each position of miRNA
Generally, the “T” base account for the highest percentage of miRNA in the first position. The two R functions (mirnapredata, basepreplot) can be used to describe miRNAs’ base distribution.

```r
## Generate the base frequency data for next drawing
P1_miRNA_data <- mirnapredata(mirnaseq_dataframe = P1_miRNA_count)
P2_miRNA_data <- mirnapredata(mirnaseq_dataframe = P2_miRNA_count)
F1_miRNA_data <- mirnapredata(mirnaseq_dataframe = F1_miRNA_count)
```
Note: The first column of the input file must be the miRNA sequence. About the output result, the first column is the base of miRNA, the second column of output file is the frequency of occurrence of the base, the third column is the position of bases.

```r
## Base preference plot of miRNA
basepreplot(file_dataframe = P1_miRNA_data)
basepreplot(file_dataframe = P2_miRNA_data)
basepreplot(file_dataframe = F1_miRNA_data)
```
## Expression analysis
### Specific expression analysis

miVennPlot: generate the Venn diagram with the specific expression information of miRNAs.

```r
## Venn Diagram
miVennPlot(P1_RPM = P1_miRNA_rpm,
           P2_RPM = P2_miRNA_rpm,
           F1_RPM = F1_miRNA_rpm,rpm_threshold = 1)
```

miVennData: Extract the species-specific miRNAs and the shared miRNAs among parents and offspring.

```r
##Extract the species-specific miRNAs and the shared miRNAs among parents and offspring.
##output_file = "venn_list"
venn_list <- miVennData(P1_RPM = P1_miRNA_rpm,
                        P2_RPM = P2_miRNA_rpm,
                        F1_RPM = F1_miRNA_rpm,
                        rpm_threshold = 1,output_file = "venn_list")
##output_file = "P1_specific"
P1_specific <- miVennData(P1_RPM = P1_miRNA_rpm,
                          P2_RPM = P2_miRNA_rpm,
                          F1_RPM = F1_miRNA_rpm,
                          rpm_threshold = 1,output_file = "P1_specific")
##output_file = "P2_specific"
P2_specific <- miVennData(P1_RPM = P1_miRNA_rpm,
                          P2_RPM = P2_miRNA_rpm,
                          F1_RPM = F1_miRNA_rpm,
                          rpm_threshold = 1,output_file = "P2_specific")
##output_file = "F1_specific"
F1_specific <- miVennData(P1_RPM = P1_miRNA_rpm,
                          P2_RPM = P2_miRNA_rpm,
                          F1_RPM = F1_miRNA_rpm,
                          rpm_threshold = 1,output_file = "F1_specific")
##output_file = "all_common"
all_common <- miVennData(P1_RPM = P1_miRNA_rpm,
                         P2_RPM = P2_miRNA_rpm,
                         F1_RPM = F1_miRNA_rpm,
                         rpm_threshold = 1,output_file = "all_common")
```

### Differential expression analysis
```r
polyDESeq(P1_count = P1_miRNA_count,
          P2_count = P2_miRNA_count,
          F1_count = F1_miRNA_count,
          count_threshold = 5,Pvalue = 0.05)
```


### Filtering low expressed miRNAs

```r
##Get the filitered mirna count table (default: Count >= 5 in at least one sample)
Count5result <- Countfiliter(P1_count = P1_miRNA_count,
                             P2_count = P2_miRNA_count,
                             F1_count = F1_miRNA_count,count_threshold = 5)
```

```r
##Get the filitered mirna rpm table (default: the average rpm >= 1 in three species)
Rpm1result <- Rpmfiliter(P1_RPM = P1_miRNA_rpm,
                         P2_RPM = P2_miRNA_rpm,
                         F1_RPM = F1_miRNA_rpm,rpm_threshold = 1)
```

## Genetic effects analysis of miRNA 

### Method1: |d/a|

```R
##Get the classification results based on the value of |d/a|
DAresult <- GetDAtable(P1_RPM = P1_miRNA_rpm,
                       P2_RPM = P2_miRNA_rpm,
                       F1_RPM = F1_miRNA_rpm,rpm_threshold = 1)
```

### Method2: Twelve bins of expression analysis

```R
##Get the table of 12 expression patterns
Binresult <- Get12Bins(P1_count = P1_miRNA_count,
                       P2_count = P2_miRNA_count,
                       F1_count = F1_miRNA_count,
                       count_threshold = 5,Pvalue = 0.05)
```
