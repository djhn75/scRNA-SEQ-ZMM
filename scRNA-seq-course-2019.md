---
title: "scRNA-SEQ-course-2019"
author: "David John & Wesley Abplanalp"
date: "12 6 2019"
---


# 1.) Install required Software
## 1.2.) Install R
1. Windows: <https://cran.rstudio.com/bin/windows/base/>
2. Mac:<https://cran.rstudio.com/bin/macosx/>
<<<<<<< HEAD
3. Linux: Please run the following commands in the terminal
=======
3. Linux: Please run the following commands in the terminal 
>>>>>>> fd0db9047a7e7327a12ca64407617b405c357330
```{shell}
sudo add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/'

sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9

sudo apt update

sudo apt install r-base r-base-core r-recommended
```


## 1.2.) Install Rstudio
To install RStudio please follow the link and install the right RStudio for your operating system.

<<<<<<< HEAD
1. Windows: <https://download1.rstudio.org/desktop/windows/RStudio-1.2.1335.exe>
2. Mac:<https://download1.rstudio.org/desktop/macos/RStudio-1.2.1335.dmg>
3. Linux: <https://www.rstudio.com/products/rstudio/download/#download>

## 1.3.) Install required R-packages
For or training we will need the following R packages, which can be installed by the following R commands
```{r cars}
# Enter commands in R (or R studio, if installed)
# Install the devtools package from Hadley Wickham
install.packages('devtools')
# Replace '2.3.0' with your desired version
devtools::install_version(package = 'Seurat', version = package_version('2.3.4'))

install.packages("tidyr")

install.packages("dplyr")

install.packages("ggpubr")

install.packages("scales")

install.packages("ggplot2")

install.packages("stringr")

```

## 1.4.) Test the required R-packages
```{r cars}
library(tidyr)
library(dplyr)
library(ggpubr)
library(scales)
library(stringr)
library(Seurat)
library(ggplot2)
=======
1. Windows: <https://www.rstudio.com/products/rstudio/download/#download>
2. Mac:<https://www.rstudio.com/products/rstudio/download/#download>
3. Linux: <https://www.rstudio.com/products/rstudio/download/#download>


```{r cars}
summary(cars)
>>>>>>> fd0db9047a7e7327a12ca64407617b405c357330
```

## Including Plots

You can also embed plots, for example:

```{r pressure, echo=FALSE}
plot(pressure)
```

Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot.