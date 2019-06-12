---
title: "scRNA-SEQ-course-2019"
author: "David John & Wesley Abplanalp"
date: "12 6 2019"
---


# 1.) Install required Software
## 1.1.) Install Rstudio
To install Rstudio please follow the link and install the right RStudio for your operating system.
<https://www.rstudio.com/products/rstudio/download/#download>
*Ubuntu
```{shell}
sudo apt-get install rstudio
````
<http://rmarkdown.rstudio.com>.

When you click the **Knit** button a document will be generated that includes both content as well as the output of any embedded R code chunks within the document. You can embed an R code chunk like this:

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

```{r cars}
summary(cars)
```

## Including Plots

You can also embed plots, for example:

```{r pressure, echo=FALSE}
plot(pressure)
```

Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot.
