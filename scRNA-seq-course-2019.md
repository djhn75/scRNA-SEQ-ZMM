---
title: "scRNA-SEQ-course-2019"
author: "David John & Wesley Abplanalp"
date: "12 6 2019"
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
1.) 2700 PBMCs (v1 Chemistry, healthy) 
*(Details: <https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/pbmc3k>)
2.) 5000 PBMCs (V3 Chemistry, healthy) 
*(Details: <https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.2/5k_pbmc_v3>)
3.) 1000 CD31 positive cells 
*(V2 Chemistry, heart failure)


## 2.1.) Load the Datasets
##### TODO:
- [ ] Download the datasets, save the folder 'scRNA-SEQ_Coursematerial_2019' to your PC:
<https://drive.google.com/open?id=1qs-Pwk0vmrjZEOUoz5kvtIsmP1KK0Czs>


## 2.2) Load the RScript
-[ ] Download the Rscript (<https://github.com/djhn75/scRNA-SEQ-ZMM/blob/master/scRNA-seq-course-2019.Rmd>) 
-[ ] Open the file in Rstudio and continue at step 3




