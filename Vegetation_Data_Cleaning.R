# 
# Author: Jake Diamond
# Purpose: To load in and clean up vegetation field data for EAB
# Date: November 16, 2016
# 

# Set Working Directory
setwd("C:/Users/Jake/Dropbox/Projects/EAB/Data")

# Load Libraries
library(zoo)
library(dplyr)
library(tidyr)
library(readxl)

# Function to read in all worksheets
read_excel_allsheets <- function(filename) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  names(x) <- sheets
  x
}

# Load data
df <- read_excel_allsheets("FieldDataforR.xlsx")
df <- do.call(rbind.data.frame, df)
df <- tibble::rownames_to_column(df)
colnames(df)[1] <- "id"
df$site <- substr(df$id, 1, 2)
df$plot <- na.locf(df$plot)
df$position <- na.locf(df$position)
df$point <- paste(df$plot, df$position, sep = ".") 
df$depth <- na.locf(df$depth)

# Load hummock hollow sheet
huho <- read.csv("hu.ho2.csv")

# Load species updates
spec <- read.csv("VegLookup.csv")

# Combine data
df <- left_join(df, select(huho, 
                           -c(id, site, plot, position)), 
                by = "point")

df <- left_join(df, spec)
df$species <- NULL
colnames(df)[15] <- "species"

# Define mosses
moss <- c(
  "Calliergon cordifolium",
  "Climacium dendroides",
  "Funaria hygrometrica",
  "Hypnum cupressiforme",
  "Pleurozium schreberi",
  "Rhizomnium magnifolium",
  "Sphagnum angustifolium",
  "Sphagnum wulfianum",
  "Thuidium delicatulum",
  "Moss1",
  "Lemna minor",
  "Liverwort")

# Create moss column
df$moss <- ifelse(df$species %in% moss, 1, 0)

# Save dataframe to preserve formatting (csv not working)
saveRDS(df, file = "veg_data.rds")


# Write to .csv
write.csv(df, "veg_data_for_analysis.csv")
