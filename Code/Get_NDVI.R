# Landsat NDVI Data
# Colton Padilla
# 7/12/2023

# Last Updated: 7/12/2023

#-----------------#
# Housekeeping ####
#-----------------#

# Packages 
library(tidyverse)
library(sf)
library(MODIStsp)
library(terra)

#-------------------------#
# Begin Modis Download ####
#-------------------------#

# Get products
MODIStsp_get_prodlayers("(M*D13A2)")

# Download the data
MODIStsp(
  gui = F,
  out_folder = "./Data/GIS_Files/",
  out_folder_mod = "$tempdir",
  selprod = "Vegetation_Indexes_16Days_1Km (M*D13A2)",
  bandsel = "NDVI",
  user = "coltonp2018",
  password = "Muledeerandelk21*",
  start_date = "1999.01.01",
  end_date = "2022.12.31",
  spatmeth = "file",
  spafile = "./Data/GIS_Files/NM_populations_v2_1.shp",
  out_format = "GTiff"
)

#------------------------#
# Process raster data ####
#------------------------#

# Create a list of files
files <- list.files("./Data/GIS_Files/NM_populations_v2_1/VI_16Days_1Km_v6/NDVI/",
                    pattern = ".tif")

# Paste directory in front of files
files <- paste0("./Data/GIS_Files/NM_populations_v2_1/VI_16Days_1Km_v6/NDVI/",
                files)

# Now create raster
ndvi <- rast(files)

# Rename the layers
names(ndvi)
names(ndvi) <- gsub("MYD13A2_NDVI_", "", names(ndvi)) %>%
  gsub("MOD13A2NDVI", "",.) %>%
  gsub("\\_", "", .)
names(ndvi)

# Load the raster and project it then rewrite
ndvi <- rast("./Data/GIS_Files/2000_2022_NDVI_NM.tif") 

# Get the NM shapefile
nm <- st_read("./Data/GIS_Files/cb_2018_us_state_500k.shp") %>%
  dplyr::filter(NAME == "New Mexico") %>%
  dplyr::select(NAME, geometry) %>%
  st_transform(., crs = crs(ndvi))

# Now crop the raster to populations
ndvi <- mask(ndvi, nm)

# Now write out the raster
writeRaster(ndvi,
            "./Data/GIS_Files/2000_2022_NDVI_NM.tif",
            gdal = "COMPRESS=DEFLATE",
            overwrite = T)
