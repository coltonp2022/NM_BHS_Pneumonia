# Prism dowloads
# Colton Padilla
# 6/26/2023

# Last Updated: 7/18/2023

#-----------------#
# Housekeeping ####
#-----------------#

# Packages
library(terra)
library(sf)
library(prism)

# Get the NM shapefile
pops <- st_read("./Data/GIS_Files/NM_populations_v2_1.shp") %>%
  dplyr::select(NAME, geometry)

# set the download folder
prism_set_dl_dir("./Data/GIS_Files/PRISM_FIX")

#-------------------------#
# Download monthly PPT ####
#-------------------------#

## Get monthly data
#get_prism_monthlys("ppt",
#                   years = 1999:2022,
#                   mon = 1:12,
#                   keepZip = F)
#
## Check for corruption
#prism_archive_verify("ppt", 
#                     temp_period = "monthly",
#                     keepZip = F)

# Get a raster stack of that data
stack <- pd_stack(prism_archive_subset(
  "ppt",
  temp_period = "monthly"
))

# Turn to terra
stack <- rast(stack)
names(stack)

# Gsub out of names
names(stack) <- gsub("PRISM_ppt_stable_4kmM3_", "", names(stack)) %>%
  gsub("PRISM_ppt_provisional_4kmM3_", "", .) %>%
  gsub("_bil", "", .)
names(stack)

# Now arrange by names
stack <- stack[[order(names(stack))]]
names(stack)

# Save this stack
writeRaster(stack,
            "./Data/GIS_FILES/1999_2022_monthly_precip_USA.tif",
            gdal = "COMPRESS=DEFLATE",
            overwrite = T)

# Project pops
pops <- st_transform(pops, crs = crs(stack))

# Now crop the raster
stack <- crop(stack, pops)

# Write again
writeRaster(stack,
            "./Data/GIS_FILES/1999_2022_monthly_precip_NM.tif",
            gdal = "COMPRESS=DEFLATE",
            overwrite = T)

#-------------------------------#
# Download Mean Temperatures ####
#-------------------------------#

# Temperatures
prism_set_dl_dir("./Data/GIS_Files/PRISM_FIX/Tmean")

## Get monthly data
#get_prism_monthlys("tmean",
#                   years = 1999:2022,
#                   mon = 1:12,
#                   keepZip = F)
#
## Check for corruption
#prism_archive_verify("tmean", 
#                     temp_period = "monthly",
#                     keepZip = F)
#

# Get a raster stack of that data
stack <- pd_stack(prism_archive_subset(
  "tmean",
  temp_period = "monthly"
))

# Turn to terra
stack <- rast(stack)
names(stack)

# Gsub out of names
names(stack) <- gsub("PRISM_tmean_stable_4kmM3_", "", names(stack)) %>%
  gsub("PRISM_tmean_provisional_4kmM3_", "", .) %>%
  gsub("_bil", "", .)
names(stack)

# Now arrange by names
stack <- stack[[order(names(stack))]]
names(stack)

# Save this stack
writeRaster(stack,
            "./Data/GIS_FILES/1999_2022_monthly_tmean_USA.tif",
            gdal = "COMPRESS=DEFLATE",
            overwrite = T)

# Now crop the raster
stack <- crop(stack, pops)

# Write again
writeRaster(stack,
            "./Data/GIS_FILES/1999_2022_monthly_tmean_NM.tif",
            gdal = "COMPRESS=DEFLATE",
            overwrite = T)

