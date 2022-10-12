## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",fig.width=12, fig.height=9, message=FALSE, tidy=TRUE, dpi=75)

## ----setup--------------------------------------------------------------------
library(iphegwas)

## -----------------------------------------------------------------------------
head(ibd,3)
head(bmi,3)
head(Wasisthipratio,3)
head(CrohnsDisease,3)
head(UlcerativeColitis,3)
## Bringing all package data to the environment
ibd <- ibd
bmi <- bmi
Wasisthipratio <- Wasisthipratio
CrohnsDisease <- CrohnsDisease
UlcerativeColitis <- UlcerativeColitis

## -----------------------------------------------------------------------------
phenos <- c("ibd","bmi","CrohnsDisease","UlcerativeColitis","Wasisthipratio")
yy <- fastprocessphegwas(phenos)

## -----------------------------------------------------------------------------
print(phenos)

## -----------------------------------------------------------------------------
landscapefast(yy,sliceval = 7,phenos =phenos)

## -----------------------------------------------------------------------------
landscapefast(yy,sliceval = 7,phenos = iphegwas(phenos))

## -----------------------------------------------------------------------------
iphegwas(phenos,dentogram = TRUE)

## -----------------------------------------------------------------------------
head(hdl,3)
head(ldl,3)
head(trig,3)
head(tchol,3)
## I am changing the name of the dataframe to something meaningful, as the name of the dataframe will be used as phenotype names in the landscape. This also bring all package data to the environment.
HDL <- hdl
LDL <- ldl
TRIGS <- trig
TOTALCHOLESTROL <- tchol

## -----------------------------------------------------------------------------
phenos <- c("HDL", "LDL", "TRIGS", "TOTALCHOLESTROL")
y <- fastprocessphegwas(phenos)

## -----------------------------------------------------------------------------
landscapefast(y,sliceval = 10,phenos =phenos)

## -----------------------------------------------------------------------------
landscapefast(y,sliceval = 7.5,chromosome = 19,phenos =phenos)

## -----------------------------------------------------------------------------
landscapefast(y,sliceval = 7.5,chromosome = 19, geneview = TRUE,phenos =phenos)

## -----------------------------------------------------------------------------
landscapefast(y, sliceval = 30, chromosome = 19,calculateLD= TRUE,mutualLD = TRUE,phenos =phenos)

