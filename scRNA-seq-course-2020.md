---
title: "scRNA-SEQ-course-2019"
author: "David John & Wesley Abplanalp"
date: "07 05 2020"
---


# 1.) Install required Software
## 1.2.) Install R

1. Windows: <https://cran.rstudio.com/bin/windows/base/>
2. Mac:<https://cran.rstudio.com/bin/macosx/>
3. Linux: Please run the following commands in the terminal
```{shell}
sudo add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/'

sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9

sudo apt update

sudo apt install r-base r-base-core r-recommended
```

## 1.2.) Install Rstudio

To install RStudio please follow the link and install the right RStudio for your operating system.

1. Windows: <https://download1.rstudio.org/desktop/windows/RStudio-1.2.1335.exe>
2. Mac:<https://download1.rstudio.org/desktop/macos/RStudio-1.2.1335.dmg>
3. Linux: <https://www.rstudio.com/products/rstudio/download/#download>

## 1.3.) Install required R-packages

For or training we will need the following R packages, which can be installed by the following R commands
```{r}
# Enter commands in R (or R studio, if installed)

install.packages("tidyr")
install.packages("dplyr")
install.packages("ggpubr")
install.packages("scales")
install.packages("ggplot2")
install.packages("stringr")
install.packages("Seurat")

```


## 1.4.) Test the required R-packages

To test all installed packages, run each command individually. 
If the installed package are loaded witout any error message, the package was installed sucessfully. 
```{r}
library(tidyr)
library(dplyr)
library(ggpubr)
library(scales)
library(stringr)
library(Seurat)
library(ggplot2)
```

# 2.) Download The required Datasets

During this course we will analyse 3 datasets. Two are punlished healthy PBMC's from 10x and one is our data, from a heart failure patient of the university clinic Frankfurt.
1. 2700 PBMCs (v1 Chemistry, healthy) 
* (Details: <https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/pbmc3k>)
2. 5000 PBMCs (V3 Chemistry, healthy) 
* (Details: <https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.2/5k_pbmc_v3>)
3. 1000 CD31 positive cells 
* (V2 Chemistry, heart failure)


## 2.1.) Load the Datasets
##### TODO:
- [ ] Download the datasets, save the folder `scRNA-SEQ_Coursematerial_2020` to your PC:
<https://drive.google.com/open?id=1hY46UrN5hzuHkRC23rOYwARU2R7jaVPi>


## 2.2) Load the RScript
- [ ] Download the Rscript (<https://github.com/djhn75/scRNA-SEQ-ZMM/blob/master/scRNA-seq-course-2020.Rmd>) 
- [ ] Open the file in Rstudio and continue at step 3



# 3.) Run the R-script: 
  <https://github.com/djhn75/scRNA-SEQ-ZMM/blob/master/scRNA-seq-course-2020.Rmd>


# 4.) Student Tasks
## 4.1.) Generate Gene Ontology Graphs with Metascape 
(<http://metascape.org/gp/index.html#/main/step1>)
### 4.1.1. Option A: 
##### Task 1
- [ ] Select Differentially expressed genes between NK and CD8 T
- [ ] Perform a Gene Ontology (GO) Analysis by Metascape with all significant upregulated marker genes for each Cell type (2 GO Analysis)
  
##### Task 2
- [ ] Subset these cell types and recluster this subset


### 4.1.2. Option B:
##### Task 1

- [ ] Select Differentially expressed genes between CD14+ Mono and FCGR3A+ Mono
- [ ] Perform a Gene Ontology (GO) Analysis by Metascape with all significant upregulated marker genes for each Cell type (2 GO Analysis)

##### Task 2
- [ ] Subset these cell types and recluster this subset
  
  
### 4.1.3. Option C:
##### Task 1

- [ ] Select Differentially expressed genes between Naive CD4 T and Memory CD4 T
- [ ] Perform a Gene Ontology (GO) Analysis by Metascape with all significant upregulated marker genes for each Cell type (2 GO Analysis)

##### Task 2
- [ ] Subset these cell types and recluster this subset
  
  
## 4.2.) Questions: Please discuss the following questions:

- [ ] How would you define the difference between the two celltypes?
- [ ] Could you find more than the initial 2 Cluster after reclustering? What is the implication of this? 
