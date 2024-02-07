packages <- c("shiny", "ggplot2", "openxlsx", 
              "gridBase", "gridExtra", "DT", 
              "htmlwidgets", "dplyr", 
              "scales")


testin <- function(package){
  if(!(package %in% installed.packages())){
    install.packages(package)
  }
}


#sapply(packages, testin)
lapply(packages, require, character.only = TRUE)
options(scipen=999)

library(ggplot2)
library(shiny)
library(openxlsx)
library(gridBase)
library(gridExtra)
library(DT)
library(htmlwidgets)
library(dplyr)
library(scales)



source("Rwanda_HIV_MSM_Functions_Quarterly.R")
params <<- openxlsx::read.xlsx("./Inputs/MSM_Model_Parameters_PrEP.xlsx", sheet=1, colNames=F)
counts <<- read.table("./Inputs/2017_Counts.txt")

nms = as.character(params[,1])
vals = as.numeric(as.character(params[,2]))
for(i in 1:length(nms)){
  assign(nms[i], vals[i])
}

SW_counts = counts[counts[,2] == "FSW", 1]
MP_counts = counts[counts[,2] == "SC", 1]
GP_counts = counts[counts[,2] == "GP", 1]
MSM_counts = counts[counts[,2] == "MSM", 1]



