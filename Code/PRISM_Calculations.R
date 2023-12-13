# Create Yearly and Summer and Winter Precip
# Colton Padilla
# 6/26/2023

# Last Updated: 6/26/2023

#-----------------#
# Housekeeping ####
#-----------------#

# Packages 
library(tidyverse)
library(sf)
library(terra)

# Load in precip data
prism <- rast("./Data/GIS_Files/1999_2022_monthly_precip_NM.tif")

#----------------------------------------#
# Create Yearly Precipitation Rasters ####
#----------------------------------------#

# Create year object
year <- seq(from = 1999, to = 2022)

# Now Create a for loop to sum each year
for(i in 1:length(year)){
  
  # Filter to the year
  df <- prism[[grepl(year[i], names(prism))]]
  
  # Now sum it all
  df <- app(df, sum)
  names(df) <- year[i]
  
  # Now write it into a stack
  if(i == 1){
    x <- df
  } else
    x <- c(x, df)
}

# Write the raster
writeRaster(x,
            "./Data/GIS_Files/Calendar_Year_Precip.tif",
            gdal = "COMPRESS=DEFLATE",
            overwrite = T)


#------------------#
# Summer Precip ####
#------------------#