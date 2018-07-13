    # Libraries
    #install.packages("heatmaply")
    #install.packages("tidyverse")

    ### Setting up environment and loading data

    library(dplyr)

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

    library(knitr)
    library(vcd)

    ## Loading required package: grid

    library(grid)

    ### Loading data

    #setwd("~//R_repository_data/scripts")
    bap <- read.csv2("KP_screening_data.csv", header = TRUE, sep = ",", na.strings = "", stringsAsFactors = FALSE)

    ### Summary
    #summary(bap)

Association analyses
--------------------

Extended mosaic plot displays show the standardized residuals of a
loglinear model of the counts from by the color and outline of the
mosaic's tiles. Standardized residuals are often referred to a standard
normal distribution. Negative residuals are drawn in shaded of red and
with broken outlines; positive ones are drawn in blue with solid
outlines.

    bap$ompK36 = as.character(bap$ompK36)
    bap$ompK35 = as.character(bap$ompK35)
    bap$ST = as.character(bap$ST)
    porines <- select(bap,Id,ST,Insertion,ompK36,ompK35,Year,Source,Source_edited,Source_summarised,Host,Country)

**Mosaic plot of the ompK35 presence/absence versus ST**

    porines$ST = as.character(porines$ST)

    porines$Insertion <- as.factor(porines$Insertion)
    levels <- levels(porines$Insertion)
    levels[length(levels) + 1] <- "-"
    porines$Insertion <- factor(porines$Insertion, levels = levels)
    porines$Insertion[is.na(porines$Insertion)] <- "-"
    porines$Insertion = as.character(porines$Insertion)

    ## contingency table
    mytableompK35 <- xtabs(~ST+ompK35, data=porines)

    ## mosaic plot
    mosaicplot(mytableompK35, shade = TRUE, las=2)

![](images/KPassoc_unnamed-chunk-3-1.png)

    # ftable creates flat contingency tables
    summary(mytableompK35)

    ## Call: xtabs(formula = ~ST + ompK35, data = porines)
    ## Number of cases in table: 1557 
    ## Number of factors: 2 
    ## Test for independence of all factors:
    ##  Chisq = 603.7, df = 227, p-value = 5.748e-36
    ##  Chi-squared approximation may be incorrect

    ## mosaic plot with expected freq
    #expected35 <- mosaic(~ST+ompK35, data=porine_data, type = "expected")
    #expected35

**Mosaic plots of the ompK36 presence/absence versus ST**

    mytableompK36 <- xtabs(~ST+ompK36, data=porines)  # 2-Way Frequency Table 
    mosaicplot(mytableompK36, shade = TRUE, las=2)

![](images/KPassoc_unnamed-chunk-4-1.png)

    summary(mytableompK36)  # chi-square test of indepedence

    ## Call: xtabs(formula = ~ST + ompK36, data = porines)
    ## Number of cases in table: 1557 
    ## Number of factors: 2 
    ## Test for independence of all factors:
    ##  Chisq = 207.51, df = 227, p-value = 0.8188
    ##  Chi-squared approximation may be incorrect

    #expected36 <- mosaic(~ST+ompK36, data=porine_data, type = "expected")
    #expected36

**Supplemental figure S8 - Mosaic plot of the ompK35/ompK36
presence/absence versus ST (3-way contigency table)**

    # 3-way contingency table
    mytableompK3536 <- xtabs(~ST+ompK35+ompK36, data=porines)
    mosaicplot(mytableompK3536, shade = TRUE, las=2)

![](images/KPassoc_unnamed-chunk-5-1.png)

    summary(mytableompK3536)

    ## Call: xtabs(formula = ~ST + ompK35 + ompK36, data = porines)
    ## Number of cases in table: 1557 
    ## Number of factors: 3 
    ## Test for independence of all factors:
    ##  Chisq = 936.3, df = 682, p-value = 2.911e-10
    ##  Chi-squared approximation may be incorrect

**Mosaic plot of the ompK35 presence/absence and Insertion in ompK36
versus ST**

    mytableinsert <- xtabs(~ST+Insertion+ompK35, data=porines)  # 2-Way Frequency Table 
    mosaicplot(mytableinsert, shade = TRUE, las=2)

![](images/KPassoc_unnamed-chunk-6-1.png)

    summary(mytableinsert)  # chi-square test of independence

    ## Call: xtabs(formula = ~ST + Insertion + ompK35, data = porines)
    ## Number of cases in table: 1557 
    ## Number of factors: 3 
    ## Test for independence of all factors:
    ##  Chisq = 1654.6, df = 1137, p-value = 6.2e-22
    ##  Chi-squared approximation may be incorrect

**Mosaic plot of the Insertion + presence/absence of ompK36 versus ST
and ompK35**

    # combine insertion and ompK36 presence/absence in one factor
    porine_data2 <- data.frame(porines[,1:5], ompK36_ins=paste(porines[,4],porines[,3],sep="-"))
    mytableinsert2 <- xtabs(~ST+ompK36_ins+ompK35, data=porine_data2)  # 2-Way Frequency Table 
    mosaicplot(mytableinsert2, shade = TRUE, las=2)

![](images/KPassoc_unnamed-chunk-7-1.png)

    summary(mytableinsert2)  # chi-square test of independence

    ## Call: xtabs(formula = ~ST + ompK36_ins + ompK35, data = porine_data2)
    ## Number of cases in table: 1557 
    ## Number of factors: 3 
    ## Test for independence of all factors:
    ##  Chisq = 1989.3, df = 1592, p-value = 3.07e-11
    ##  Chi-squared approximation may be incorrect

**Mosaic plot of the Insertion + presence/absence of ompK36 versus ST
and ompK35 (significant only)**

    # filter significantly different STs based on previous plot
    significant <- c('-', '101', '11', '14', '147', '278', '280', '15', '16', '17', '120', '200', '231', '258', '127', '273', '277', '29', '278', '3', '336', '37', '392', '442', '512', '873')
    porine_data2_significant <- filter(porine_data2, ST %in% significant) %>% na.omit()
    mytableinsert3 <- xtabs(~ST+ompK36_ins+ompK35, data=porine_data2_significant)  # 2-Way Frequency Table 

    mosaicplot(mytableinsert3, shade = TRUE, las=2, cex.axis = 1.5, type="pearson")

![](images/KPassoc_unnamed-chunk-8-1.png)

    summary(mytableinsert3)   # chi-square test of independence

    ## Call: xtabs(formula = ~ST + ompK36_ins + ompK35, data = porine_data2_significant)
    ## Number of cases in table: 1099 
    ## Number of factors: 3 
    ## Test for independence of all factors:
    ##  Chisq = 1081.3, df = 157, p-value = 9.181e-138
    ##  Chi-squared approximation may be incorrect

**Supplemental figure S9 - Mosaic plots of ompK35 and the Insertion +
presence/absence of ompK36 versus year **

    porine_data_date <- select(porines,-Id)
    porine_data_dateins <- data.frame(porine_data_date[,1:5], ompK36_ins=paste(porine_data_date[,3],porine_data_date[,2],sep="-"))

    ## OmpK36_ins by year
    date_ompK36 <- xtabs(~Year+ompK36_ins, data=porine_data_dateins)
    mosaicplot(date_ompK36, shade = TRUE, cex.axis = 1, las=2)

![](images/KPassoc_unnamed-chunk-9-1.png)

    summary(date_ompK36)

    ## Call: xtabs(formula = ~Year + ompK36_ins, data = porine_data_dateins)
    ## Number of cases in table: 1157 
    ## Number of factors: 2 
    ## Test for independence of all factors:
    ##  Chisq = 288.87, df = 48, p-value = 4.041e-36
    ##  Chi-squared approximation may be incorrect

    ## OmpK35 by year
    date_ompK35 <- xtabs(~Year+ompK35, data=porine_data_dateins)
    mosaicplot(date_ompK35, shade = TRUE, cex.axis = 1, las=2)

![](images/KPassoc_unnamed-chunk-9-2.png)

    summary(date_ompK35)

    ## Call: xtabs(formula = ~Year + ompK35, data = porine_data_dateins)
    ## Number of cases in table: 1157 
    ## Number of factors: 2 
    ## Test for independence of all factors:
    ##  Chisq = 157.25, df = 16, p-value = 2.884e-25
    ##  Chi-squared approximation may be incorrect

**Supplemental figure S10 - Mosaic plots of ompK35 and the Insertion +
presence/absence of ompK36 versus country of origin **

    porine_data_country <- select(porines,-Id)
    porine_data_countryins <- data.frame(porine_data_country[,1:5], ompK36_ins=paste(porine_data_country[,3],porine_data_country[,2],sep="-"), porine_data_country[,7:10])

    ## OmpK36_ins by country
    country_ompK36 <- xtabs(~Country+ompK36_ins, data=porine_data_countryins)
    mosaicplot(country_ompK36, shade = TRUE, cex.axis = 1, las=2)

![](images/KPassoc_unnamed-chunk-10-1.png)

    summary(country_ompK36)

    ## Call: xtabs(formula = ~Country + ompK36_ins, data = porine_data_countryins)
    ## Number of cases in table: 1234 
    ## Number of factors: 2 
    ## Test for independence of all factors:
    ##  Chisq = 640.2, df = 108, p-value = 1.632e-76
    ##  Chi-squared approximation may be incorrect

    ## OmpK35 by country
    country_ompK35 <- xtabs(~Country+ompK35, data=porine_data_countryins)
    mosaicplot(country_ompK35, shade = TRUE, cex.axis = 1, las=2)

![](images/KPassoc_unnamed-chunk-10-2.png)

    summary(country_ompK35)

    ## Call: xtabs(formula = ~Country + ompK35, data = porine_data_countryins)
    ## Number of cases in table: 1234 
    ## Number of factors: 2 
    ## Test for independence of all factors:
    ##  Chisq = 338.1, df = 36, p-value = 8.896e-51
    ##  Chi-squared approximation may be incorrect

**Mosaic plots of ompK35 and the Insertion + presence/absence of ompK36
versus source (cleaned up) **

    porine_data_source <- select(porines,-Id)
    porine_data_sourceins <- data.frame(porine_data_country[,1:5], ompK36_ins=paste(porine_data_country[,3],porine_data_country[,2],sep="-"), porine_data_country[,7:10])

    ## OmpK36_ins by Source
    #source_ompK36 <- xtabs(~Source_summarised+ompK36_ins, data=porine_data_sourceins)
    source_ompK36 <- xtabs(~Source_edited+ompK36_ins, data=porine_data_sourceins)
    mosaicplot(source_ompK36, shade = TRUE, cex.axis = 1, las=2)

![](images/KPassoc_unnamed-chunk-11-1.png)

    summary(source_ompK36)

    ## Call: xtabs(formula = ~Source_edited + ompK36_ins, data = porine_data_sourceins)
    ## Number of cases in table: 994 
    ## Number of factors: 2 
    ## Test for independence of all factors:
    ##  Chisq = 421.2, df = 291, p-value = 8.798e-07
    ##  Chi-squared approximation may be incorrect

    ## OmpK35 by Source
    source_ompK35 <- xtabs(~Source_edited+ompK35, data=porine_data_sourceins)
    mosaicplot(source_ompK35, shade = TRUE, cex.axis = 1, las=2)

![](images/KPassoc_unnamed-chunk-11-2.png)

    summary(source_ompK35)

    ## Call: xtabs(formula = ~Source_edited + ompK35, data = porine_data_sourceins)
    ## Number of cases in table: 994 
    ## Number of factors: 2 
    ## Test for independence of all factors:
    ##  Chisq = 203.54, df = 97, p-value = 1.506e-09
    ##  Chi-squared approximation may be incorrect

**Mosaic plots of ompK35 and the Insertion + presence/absence of ompK36
versus source (combined in 6 major categories) **

    porine_data_source <- select(porines,-Id)
    porine_data_sourceins <- data.frame(porine_data_country[,1:5], ompK36_ins=paste(porine_data_country[,3],porine_data_country[,2],sep="-"), porine_data_country[,7:10])

    ## OmpK36_ins by Source
    source_ompK36 <- xtabs(~Source_summarised+ompK36_ins, data=porine_data_sourceins)
    mosaicplot(source_ompK36, shade = TRUE, cex.axis = 1, las=2)

![](images/KPassoc_unnamed-chunk-12-1.png)

    summary(source_ompK36)

    ## Call: xtabs(formula = ~Source_summarised + ompK36_ins, data = porine_data_sourceins)
    ## Number of cases in table: 994 
    ## Number of factors: 2 
    ## Test for independence of all factors:
    ##  Chisq = 79.26, df = 15, p-value = 9.535e-11
    ##  Chi-squared approximation may be incorrect

    ## OmpK35 by Source
    source_ompK35 <- xtabs(~Source_summarised+ompK35, data=porine_data_sourceins)
    mosaicplot(source_ompK35, shade = TRUE, cex.axis = 1, las=2)

![](images/KPassoc_unnamed-chunk-12-2.png)

    summary(source_ompK35)

    ## Call: xtabs(formula = ~Source_summarised + ompK35, data = porine_data_sourceins)
    ## Number of cases in table: 994 
    ## Number of factors: 2 
    ## Test for independence of all factors:
    ##  Chisq = 69.86, df = 5, p-value = 1.095e-13

**Mosaic plot of ompK35, the Insertion + presence/absence of ompK36
versus source (combined in 6 major categories) **

    source_ompK3536 <- xtabs(~Source_summarised+ompK36_ins+ompK35, data=porine_data_sourceins)
    mosaicplot(source_ompK3536, shade = TRUE, las=2)

![](images/KPassoc_unnamed-chunk-13-1.png)

    #summary(source_ompK3536)
