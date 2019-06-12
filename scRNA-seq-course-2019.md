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

1. Windows: <https://www.rstudio.com/products/rstudio/download/#download>
2. Mac:<https://www.rstudio.com/products/rstudio/download/#download>
3. Linux: <https://www.rstudio.com/products/rstudio/download/#download>


```{r cars}
summary(cars)
```

## Including Plots

You can also embed plots, for example:

```{r pressure, echo=FALSE}
plot(pressure)
```

Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot.
