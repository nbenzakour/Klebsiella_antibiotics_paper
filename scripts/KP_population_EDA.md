# KP_population_EDA
Nouri L. BEN ZAKOUR  
12 July 2018  





```r
# Libraries
#install.packages("tidyverse")

### Setting up environment and loading data

#
library(dplyr)
```

```
## 
## Attaching package: 'dplyr'
```

```
## The following objects are masked from 'package:stats':
## 
##     filter, lag
```

```
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
```

```r
library(ggplot2)
library(knitr)
library(gridExtra)
```

```
## 
## Attaching package: 'gridExtra'
```

```
## The following object is masked from 'package:dplyr':
## 
##     combine
```

```r
### Loading data

bap <- read.csv2("../data/KP_screening_data.csv", header = TRUE, sep = ",", na.strings = "", stringsAsFactors = FALSE)

### Summary
#summary(bap)
```


Population features
-------------------

**Supplemental figure S7**



```r
# sort
bap <- within(bap, ST <- factor(ST, levels=names(sort(table(ST), decreasing=TRUE))))
bubble_chart_3 <- function (data1, arg1, arg2, arg3, minfreq) {
  ### summarise frequencies in data1 for each combination of arg1, arg2 and arg3
  profile_summary <- eval(substitute(as.data.frame(with(data1, table(arg1, arg2, arg3)))))
  ### filter by freq cut-off of unique occurence
  profile_present <- filter(profile_summary, Freq > minfreq)
  ### bubble chart
  p <- eval(substitute(ggplot(profile_present, aes(y = arg2, x = arg3, fill = arg1)) + geom_jitter(aes(size=Freq, alpha=1, colour = arg1), width = 0.12, height = 0.12) + scale_size_area(max_size = 26) + theme(axis.text.x=element_text(angle=30, hjust=1, size=18), axis.title.x=element_text(size=20), axis.text.y=element_text(hjust=1, size=18), axis.title.y=element_text(size=20)) + scale_color_hue(l=60, c=300)))
  return(p)
}

bubble_chart_3(bap, ST, Country, Year, 1)
```

![](../images/KP_unnamed-chunk-2-1.png)<!-- -->

**Exploring further: demographics metrics**



```r
# sort
bap <- within(bap, ST <- factor(ST, levels=names(sort(table(ST), decreasing=TRUE))))
bap30 <- within(bap, ST <- factor(ST, levels=names(sort(table(ST), decreasing=TRUE)[1:30])))

p1 <- ggplot(bap, aes(ST, fill=Host)) + geom_bar(position = position_stack(reverse = TRUE), colour="black", size=0.2) + theme(axis.text.x=element_text(angle=60, hjust=1)) + ggtitle("Distribution of host by ST")
p2 <- ggplot(bap30, aes(ST, fill=Host)) + geom_bar(position = position_stack(reverse = TRUE), colour="black", size=0.2) + theme(axis.text.x=element_text(angle=60, hjust=1)) + ggtitle("Distribution of host by ST - Top 30")

grid.arrange(p1,p2,ncol=1)
```

![](../images/KP_unnamed-chunk-3-1.png)<!-- -->


```r
p1 <- ggplot(bap, aes(ST, fill=Country)) + geom_bar(position = position_stack(reverse = TRUE), colour="black", size=0.2) + theme(axis.text.x=element_text(angle=60, hjust=1)) + ggtitle("Distribution of country by ST")
p2 <- ggplot(bap30, aes(ST, fill=Country)) + geom_bar(position = position_stack(reverse = TRUE), colour="black", size=0.2) + theme(axis.text.x=element_text(angle=60, hjust=1)) + ggtitle("Distribution of country by ST - Top 30")

grid.arrange(p1,p2,ncol=1)
```

![](../images/KP_unnamed-chunk-4-1.png)<!-- -->


```r
bap30 <- within(bap30, Country <- factor(Country, levels=names(sort(table(Country), decreasing=TRUE)[1:30])))
ggplot(bap30, aes(Country, fill=ST)) + geom_bar(position = position_stack(reverse = TRUE), colour="black", size=0.2) + theme(axis.text.x=element_text(angle=60, hjust=1)) + ggtitle("Distribution of ST by country - Top 30") 
```

![](../images/KP_unnamed-chunk-5-1.png)<!-- -->


```r
p1 <- ggplot(bap, aes(ST, fill=Insertion)) + geom_bar(position = position_stack(reverse = TRUE), colour="black", size=0.2) + theme(axis.text.x=element_text(angle=60, hjust=1)) + ggtitle("Distribution of GD/TD insertion by ST")
p2 <- ggplot(bap30, aes(ST, fill=Insertion)) + geom_bar(position = position_stack(reverse = TRUE), colour="black", size=0.2) + theme(axis.text.x=element_text(angle=60, hjust=1)) + ggtitle("Distribution of GD/TD insertion by ST - Top 30")

grid.arrange(p1,p2,ncol=1)
```

![](../images/KP_unnamed-chunk-6-1.png)<!-- -->



```r
### Source cleaned up and summarised before plotting
p1 <- ggplot(bap, aes(ST, fill=Source_summarised)) + geom_bar(position = position_stack(reverse = TRUE), colour="black", size=0.2) + theme(axis.text.x=element_text(angle=60, hjust=1)) + ggtitle("Distribution of source by ST") 
p2 <- ggplot(bap30, aes(ST, fill=Source_summarised)) + geom_bar(position = position_stack(reverse = TRUE), colour="black", size=0.2) + theme(axis.text.x=element_text(angle=60, hjust=1)) + ggtitle("Distribution of source by ST - Top 30") 

grid.arrange(p1,p2,ncol=1)
```

![](../images/KP_unnamed-chunk-7-1.png)<!-- -->


```r
p1 <- bubble_chart_3(bap, ST, Insertion, Country, 1)
p2 <- bubble_chart_3(bap, ST, Insertion, Year, 1)
grid.arrange(p1,p2,ncol=1)
```

![](../images/KP_unnamed-chunk-8-1.png)<!-- -->
