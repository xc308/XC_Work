#********************************#
# Space-time Exploratory Analysis
#********************************#

getwd()
source("~/OneDrive - University of Exeter/XC_PhD/Data/Processed/XC_WORK/Data_for_process_for_explore.RData")


install.packages(c("fields",  "CCA", "dplyr", "tidyr",
                   "ggplot2", "gstat", "sp", "spacetime"))
yes

install.packages("gstat")
install.packages("spacetime")

library(fields)
library(CCA)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gstat)  # unsucess
library(sp)
library(spacetime) # unsucess


install.packages(c("devtools", "usethis"))
library(usethis)
library(devtools)

install_github("andrewzm/STRbook")
library(STRbook)




