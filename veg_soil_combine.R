# 
# Author: Jake Diamond
# Purpose: Combine elevation data with veg and chem data
# Date: December 6, 2018
# 

# Set Working Directory
# setwd("E:/Dropbox/Dropbox/Projects/EAB/Data")
setwd("C:/Users/diamo/Dropbox/Projects/EAB/Data")

# Load Libraries
library(tidyverse)
library(broom)
library(lubridate)

# Load data
elev <- read.csv("relative_elevations.csv",
                 stringsAsFactors = FALSE)
elev$X <- NULL
huho <- read.csv("hu.ho.csv", stringsAsFactors = FALSE)

# Data cleaning
huho$point <- paste(huho$plot, huho$position, sep = ".")
huho[huho$site == "B1", "site"] <- "L1"
huho[huho$site == "B3", "site"] <- "L2"
huho[huho$site == "B6", "site"] <- "L3"

# Join data
df <- left_join(huho, elev, by = c("site", "plot" , "point"))