# GIS Extractions
# Colton Padilla
# 6/26/2023

# Last Updated: 7/17/2023

#-----------------#
# Housekeeping ####
#-----------------#

# Packages
library(tidyverse)
library(sf)
library(terra)
library(pbapply)

rm(list = ls())

# Load in NDVI
ndvi <- rast("./Data/GIS_Files/2000_2022_Monthly_NDVI_NM.tif")

# Load in the PRISM data
prism <- rast("./Data/GIS_Files/1999_2022_monthly_precip_NM.tif")

# Yearly
yearly <- rast("./Data/GIS_Files/Calendar_Year_Precip.tif")

# Load in temp
tmean <- rast("./Data/GIS_Files/1999_2022_monthly_tmean_NM.tif")

## ONLY RUN IF NEEDED
#
## Rename ndvi layers
#names(ndvi)
#names(ndvi) <- as.Date(names(ndvi), "%Y%j") %>%
#  format.Date(., "%Y%m")
#
## Order ndvi
#ndvi <- ndvi[[order(names(ndvi))]]
#
## Average ndvi by year months
#names <- unique(names(ndvi))
#ndvi <- do.call(c, lapply(1:length(names), function(i){
#  # Subset to those rasters
#  x <- ndvi[[names(ndvi) == names[i]]]
#  # Average them
#  x <- app(x, fun = mean, cores = 7)
#  # Add names
#  names(x) <- names[i]
#  print(i)
#  #return
#  return(x)
#}))
#
## Write
#writeRaster(ndvi,
#            "./Data/GIS_Files/2000_2022_Monthly_NDVI_NM.tif",
#            gdal = "COMPRESS=DEFLATE",
#            overwrite = T)

#-----------------------#
# Raster projections ####
#-----------------------#

# ONLY RERUN IF NEEDED

## Project all to ndvi crs
#prism <- project(prism, crs(ndvi))
#yearly <- project(yearly, crs(ndvi))
#tmean <- project(tmean, crs(ndvi))
#
## Write
#writeRaster(prism,
#            "./Data/GIS_Files/1999_2022_monthly_precip_NM.tif",
#            gdal = "COMPRESS=DEFLATE",
#            overwrite = T)
#
## Write
#writeRaster(yearly,
#            "./Data/GIS_Files/Calendar_Year_Precip.tif",
#            gdal = "COMPRESS=DEFLATE",
#            overwrite = T)
#
## Write
#writeRaster(tmean,
#            "./Data/GIS_Files/1999_2022_monthly_tmean_NM.tif",
#            gdal = "COMPRESS=DEFLATE",
#            overwrite = T)

#-----------------------------------------------------#
# Load in lamb ewe data and join population shapes ####
#-----------------------------------------------------#

# Lamb Ewe Data
lamb <- read_csv("./Data/Clean/Lamb_Population_Size.csv") %>% 
  arrange(Herd, Year) %>%
  group_by(Herd, Year) %>%
  slice_tail(n = 1)


# Pull in the population shapefiles
pop <- st_read("./Data/GIS_Files/NM_populations_v2_1.shp") %>%
  st_transform(crs = crs(ndvi)) %>%
  mutate(Area = as.numeric(st_area(geometry) / 100000))

# Now join geometry by abreviation
lamb <- right_join(lamb, 
                  pop %>% dplyr::select(Abrev, Area, geometry), 
                  by = c("Herd" = "Abrev")) %>%
  st_as_sf(sf_column_name = "geometry")

# Remove empty polygons for now
lamb <- lamb[!st_is_empty(lamb),, drop = FALSE]

# Add subspecies
unique(lamb$Herd)
lamb$Subspec <- ifelse(lamb$Herd %in% c("PEC", "WH", "SFR", "TC",
                                        "DC", "RGG", "MZ", "CUL",
                                        "LAT", "JZ"), "Rocky",
                       "Desert")


#----------------------------#
# Lamb month manipulation ####
#----------------------------#

# Months and numbers
months <- data.frame(
  month = c("JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", 
            "OCT", "NOV", "DEC"),
  num = c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")
)

# If else new column
for(i in 1:nrow(months)){
  if(i == 1){
    lamb$MonthNum <- ifelse(lamb$Month == months$month[i],
                            months$num[i],
                            lamb$Month)
  } else
    lamb$MonthNum <- ifelse(lamb$Month == months$month[i],
                            months$num[i],
                            lamb$MonthNum)
}

# Now paste year and month together
lamb$YearMonth <- paste0(lamb$Year, lamb$MonthNum)

# Add in rocky time frames
lamb$Timeline <- ifelse(lamb$Herd %in% c("PEC", "WH", "CUL", "LAT"),
                        6, 
                        ifelse(lamb$Herd %in% c("SFR", "TC", "DC", "RGG",
                                                "MZ", "JZ"), 5, 1))

#----------------------------------#
# 3 Months prior to parturition ####
#----------------------------------#

# month year
year <- unique(lamb$Year) %>% sort()

# Subspecies
sub <- unique(lamb$Subspec)

# Get timelines
time <- unique(lamb$Timeline)[2:3]

# Lapply
lamb <- do.call(rbind, pblapply(1:length(year), function(i){
    do.call(rbind, lapply(1:length(sub), function(j){
        # Subset data to a year and subspecies
        df <- lamb %>% filter(Year == year[i] & Subspec == sub[j])
        
        # If subspec
        if(unique(df$Subspec) == "Desert"){
          # Now get previous year
          prev <- year[i] - 1
          # Paste the three months
          months <- paste0(prev, 10:12)
          
          # Subset the raster
          x <- prism[[names(prism) %in% months]]
          
          # Sum
          x <- app(x, fun = sum, na.rm = T)
          
          # Now extract
          df <- terra::extract(x, vect(df), mean, bind = T) %>%
            as.data.frame(.) %>%
            dplyr::select(Herd, Year, ncol(df))  %>%
            rename(PrePart_Precip = ncol(.))
        } else {
          # Get two dataframes with two timelines
          df1 <- df %>% filter(Timeline == time[1])
          df2 <- df %>% filter(Timeline == time[2])
          
          # Now create two separate months
          months1 <- paste0(year[i], 0, 3:5)
          months2 <- paste0(year[i], 0, 2:4)
          
          # Subset the raster
          x1 <- prism[[names(prism) %in% months1]]
          x2 <- prism[[names(prism) %in% months2]]
          
          # Sum
          x1 <- app(x1, fun = sum, na.rm = T)
          x2 <- app(x2, fun = sum, na.rm = T)
          
          # Now extract
          df1 <- terra::extract(x1, vect(df1), mean, bind = T) %>%
            as.data.frame(.) %>%
            dplyr::select(Herd, Year, ncol(df1))  %>%
            rename(PrePart_Precip = ncol(.))
          df2 <- terra::extract(x2, vect(df2), mean, bind = T) %>%
            as.data.frame(.) %>%
            dplyr::select(Herd, Year, ncol(df2))  %>%
            rename(PrePart_Precip = ncol(.))
          # Bind them back together
          df <- rbind(df1, df2)
        }
        # Return df
        return(df)
  }))
})) %>% right_join(lamb, ., by = c("Herd", "Year"))

#------------------------------------#
# Previous Recovery Period Precip ####
#------------------------------------#

# Lapply
lamb <- do.call(rbind, pblapply(1:length(year), function(i){
  do.call(rbind, lapply(1:length(sub), function(j){
    # Subset data to a year and subspecies
    df <- lamb %>% filter(Year == year[i] & Subspec == sub[j])
    
    # If subspec
    if(unique(df$Subspec) == "Desert"){
      # Now get previous year
      prev <- year[i] - 1
      # Paste the recovery months
      months <- paste0(prev, 0,  7:9)
      
      # Subset the raster
      x <- prism[[names(prism) %in% months]]
      
      # Sum
      x <- app(x, fun = sum, na.rm = T)
      
      # Now extract
      df <- terra::extract(x, vect(df), mean, bind = T) %>%
        as.data.frame(.) %>%
        dplyr::select(Herd, Year, ncol(df))  %>%
        rename(Recovery_Precip = ncol(.))
    } else {
      # Get two dataframes with two timelines
      df1 <- df %>% filter(Timeline == time[1])
      df2 <- df %>% filter(Timeline == time[2])
      
      # Now create two separate months
      months1 <- c("09", "10", "11")
      months1 <- paste0(year[i] - 1, months1)
      months2 <- c("08", "09", "10")
      months2 <- paste0(year[i] - 1, months2)
      
      # Subset the raster
      x1 <- prism[[names(prism) %in% months1]]
      x2 <- prism[[names(prism) %in% months2]]
      
      # Sum
      x1 <- app(x1, fun = sum, na.rm = T)
      x2 <- app(x2, fun = sum, na.rm = T)
      
      # Now extract
      df1 <- terra::extract(x1, vect(df1), mean, bind = T) %>%
        as.data.frame(.) %>%
        dplyr::select(Herd, Year, ncol(df1))  %>%
        rename(Recovery_Precip = ncol(.))
      df2 <- terra::extract(x2, vect(df2), mean, bind = T) %>%
        as.data.frame(.) %>%
        dplyr::select(Herd, Year, ncol(df2))  %>%
        rename(Recovery_Precip = ncol(.))
      # Bind them back together
      df <- rbind(df1, df2)
    }
    # Return df
    return(df)
  }))
})) %>% right_join(lamb, ., by = c("Herd", "Year"))

#---------------------------#
# Postparturition Precip ####
#---------------------------#

# Lapply
lamb <- do.call(rbind, pblapply(1:length(year), function(i){
  do.call(rbind, lapply(1:length(sub), function(j){
    # Subset data to a year and subspecies
    df <- lamb %>% filter(Year == year[i] & Subspec == sub[j])
    
    # If subspec
    if(unique(df$Subspec) == "Desert"){
      # Paste the recovery months
      months <- paste0(year[i], 0,  3:5)
      
      # Subset the raster
      x <- prism[[names(prism) %in% months]]
      
      # Sum
      x <- app(x, fun = sum, na.rm = T)
      
      # Now extract
      df <- terra::extract(x, vect(df), mean, bind = T) %>%
        as.data.frame(.) %>%
        dplyr::select(Herd, Year, ncol(df))  %>%
        rename(PostPart_Precip = ncol(.))
    } else {
      # Get two dataframes with two timelines
      df1 <- df %>% filter(Timeline == time[1])
      df2 <- df %>% filter(Timeline == time[2])
      
      # Now create two separate months
      months1 <- c("06", "07", "08")
      months1 <- paste0(year[i], months1)
      months2 <- c("05", "06", "07")
      months2 <- paste0(year[i], months2)
      
      # Subset the raster
      x1 <- prism[[names(prism) %in% months1]]
      x2 <- prism[[names(prism) %in% months2]]
      
      # Sum
      x1 <- app(x1, fun = sum, na.rm = T)
      x2 <- app(x2, fun = sum, na.rm = T)
      
      # Now extract
      df1 <- terra::extract(x1, vect(df1), mean, bind = T) %>%
        as.data.frame(.) %>%
        dplyr::select(Herd, Year, ncol(df1))  %>%
        rename(PostPart_Precip = ncol(.))
      df2 <- terra::extract(x2, vect(df2), mean, bind = T) %>%
        as.data.frame(.) %>%
        dplyr::select(Herd, Year, ncol(df2))  %>%
        rename(PostPart_Precip = ncol(.))
      # Bind them back together
      df <- rbind(df1, df2)
    }
    # Return df
    return(df)
  }))
})) %>% right_join(lamb, ., by = c("Herd", "Year"))

#----------------------------------------------#
# 3 Months prior to parturition temperature ####
#----------------------------------------------#

# Lapply
lamb <- do.call(rbind, pblapply(1:length(year), function(i){
  do.call(rbind, lapply(1:length(sub), function(j){
    # Subset data to a year and subspecies
    df <- lamb %>% filter(Year == year[i] & Subspec == sub[j])
    
    # If subspec
    if(unique(df$Subspec) == "Desert"){
      # Now get previous year
      prev <- year[i] - 1
      # Paste the three months
      months <- paste0(prev, 10:12)
      
      # Subset the raster
      x <- tmean[[names(tmean) %in% months]]
      
      # Sum
      x <- app(x, fun = mean, na.rm = T)
      
      # Now extract
      df <- terra::extract(x, vect(df), mean, bind = T) %>%
        as.data.frame(.) %>%
        dplyr::select(Herd, Year, ncol(df))  %>%
        rename(PrePart_Temp = ncol(.))
    } else {
      # Get two dataframes with two timelines
      df1 <- df %>% filter(Timeline == time[1])
      df2 <- df %>% filter(Timeline == time[2])
      
      # Now create two separate months
      months1 <- paste0(year[i], 0, 3:5)
      months2 <- paste0(year[i], 0, 2:4)
      
      # Subset the raster
      x1 <- tmean[[names(tmean) %in% months1]]
      x2 <- tmean[[names(tmean) %in% months2]]
      
      # Sum
      x1 <- app(x1, fun = mean, na.rm = T)
      x2 <- app(x2, fun = mean, na.rm = T)
      
      # Now extract
      df1 <- terra::extract(x1, vect(df1), mean, bind = T) %>%
        as.data.frame(.) %>%
        dplyr::select(Herd, Year, ncol(df1))  %>%
        rename(PrePart_Temp = ncol(.))
      df2 <- terra::extract(x2, vect(df2), mean, bind = T) %>%
        as.data.frame(.) %>%
        dplyr::select(Herd, Year, ncol(df2))  %>%
        rename(PrePart_Temp = ncol(.))
      # Bind them back together
      df <- rbind(df1, df2)
    }
    # Return df
    return(df)
  }))
})) %>% right_join(lamb, ., by = c("Herd", "Year"))

#------------------------------------#
# Previous Recovery Period Precip ####
#------------------------------------#

# Lapply
lamb <- do.call(rbind, pblapply(1:length(year), function(i){
  do.call(rbind, lapply(1:length(sub), function(j){
    # Subset data to a year and subspecies
    df <- lamb %>% filter(Year == year[i] & Subspec == sub[j])
    
    # If subspec
    if(unique(df$Subspec) == "Desert"){
      # Now get previous year
      prev <- year[i] - 1
      # Paste the recovery months
      months <- paste0(prev, 0,  7:9)
      
      # Subset the raster
      x <- tmean[[names(tmean) %in% months]]
      
      # Sum
      x <- app(x, fun = mean, na.rm = T)
      
      # Now extract
      df <- terra::extract(x, vect(df), mean, bind = T) %>%
        as.data.frame(.) %>%
        dplyr::select(Herd, Year, ncol(df))  %>%
        rename(Recovery_Temp = ncol(.))
    } else {
      # Get two dataframes with two timelines
      df1 <- df %>% filter(Timeline == time[1])
      df2 <- df %>% filter(Timeline == time[2])
      
      # Now create two separate months
      months1 <- c("09", "10", "11")
      months1 <- paste0(year[i] - 1, months1)
      months2 <- c("08", "09", "10")
      months2 <- paste0(year[i] - 1, months2)
      
      # Subset the raster
      x1 <- tmean[[names(tmean) %in% months1]]
      x2 <- tmean[[names(tmean) %in% months2]]
      
      # Sum
      x1 <- app(x1, fun = mean, na.rm = T)
      x2 <- app(x2, fun = mean, na.rm = T)
      
      # Now extract
      df1 <- terra::extract(x1, vect(df1), mean, bind = T) %>%
        as.data.frame(.) %>%
        dplyr::select(Herd, Year, ncol(df1))  %>%
        rename(Recovery_Temp = ncol(.))
      df2 <- terra::extract(x2, vect(df2), mean, bind = T) %>%
        as.data.frame(.) %>%
        dplyr::select(Herd, Year, ncol(df2))  %>%
        rename(Recovery_Temp = ncol(.))
      # Bind them back together
      df <- rbind(df1, df2)
    }
    # Return df
    return(df)
  }))
})) %>% right_join(lamb, ., by = c("Herd", "Year"))

#--------------------------------#
# Postparturition Temperature ####
#--------------------------------#

# Lapply
lamb <- do.call(rbind, pblapply(1:length(year), function(i){
  do.call(rbind, lapply(1:length(sub), function(j){
    # Subset data to a year and subspecies
    df <- lamb %>% filter(Year == year[i] & Subspec == sub[j])
    
    # If subspec
    if(unique(df$Subspec) == "Desert"){
      # Paste the recovery months
      months <- paste0(year[i], 0,  3:5)
      
      # Subset the raster
      x <- tmean[[names(tmean) %in% months]]
      
      # Sum
      x <- app(x, fun = mean, na.rm = T)
      
      # Now extract
      df <- terra::extract(x, vect(df), mean, bind = T) %>%
        as.data.frame(.) %>%
        dplyr::select(Herd, Year, ncol(df))  %>%
        rename(PostPart_Temp = ncol(.))
    } else {
      # Get two dataframes with two timelines
      df1 <- df %>% filter(Timeline == time[1])
      df2 <- df %>% filter(Timeline == time[2])
      
      # Now create two separate months
      months1 <- c("06", "07", "08")
      months1 <- paste0(year[i], months1)
      months2 <- c("05", "06", "07")
      months2 <- paste0(year[i], months2)
      
      # Subset the raster
      x1 <- tmean[[names(tmean) %in% months1]]
      x2 <- tmean[[names(tmean) %in% months2]]
      
      # Sum
      x1 <- app(x1, fun = mean, na.rm = T)
      x2 <- app(x2, fun = mean, na.rm = T)
      
      # Now extract
      df1 <- terra::extract(x1, vect(df1), mean, bind = T) %>%
        as.data.frame(.) %>%
        dplyr::select(Herd, Year, ncol(df1))  %>%
        rename(PostPart_Temp = ncol(.))
      df2 <- terra::extract(x2, vect(df2), mean, bind = T) %>%
        as.data.frame(.) %>%
        dplyr::select(Herd, Year, ncol(df2))  %>%
        rename(PostPart_Temp = ncol(.))
      # Bind them back together
      df <- rbind(df1, df2)
    }
    # Return df
    return(df)
  }))
})) %>% right_join(lamb, ., by = c("Herd", "Year"))

#-----------------------------------------------------#
# Previous year lamb season Postparturition Precip ####
#-----------------------------------------------------#

# Lapply
lamb <- do.call(rbind, pblapply(1:length(year), function(i){
  do.call(rbind, lapply(1:length(sub), function(j){
    # Subset data to a year and subspecies
    df <- lamb %>% filter(Year == year[i] & Subspec == sub[j])
    
    # If subspec
    if(unique(df$Subspec) == "Desert"){
      # Paste the recovery months
      months <- paste0(year[i] - 1, 0,  3:5)
      
      # Subset the raster
      x <- prism[[names(prism) %in% months]]
      
      # Sum
      x <- app(x, fun = sum, na.rm = T)
      
      # Now extract
      df <- terra::extract(x, vect(df), mean, bind = T) %>%
        as.data.frame(.) %>%
        dplyr::select(Herd, Year, ncol(df))  %>%
        rename(PrevLamb_Precip = ncol(.))
    } else {
      # Get two dataframes with two timelines
      df1 <- df %>% filter(Timeline == time[1])
      df2 <- df %>% filter(Timeline == time[2])
      
      # Now create two separate months
      months1 <- c("06", "07", "08")
      months1 <- paste0(year[i] - 1, months1)
      months2 <- c("05", "06", "07")
      months2 <- paste0(year[i] - 1, months2)
      
      # Subset the raster
      x1 <- prism[[names(prism) %in% months1]]
      x2 <- prism[[names(prism) %in% months2]]
      
      # Sum
      x1 <- app(x1, fun = sum, na.rm = T)
      x2 <- app(x2, fun = sum, na.rm = T)
      
      # Now extract
      df1 <- terra::extract(x1, vect(df1), mean, bind = T) %>%
        as.data.frame(.) %>%
        dplyr::select(Herd, Year, ncol(df1))  %>%
        rename(PrevLamb_Precip = ncol(.))
      df2 <- terra::extract(x2, vect(df2), mean, bind = T) %>%
        as.data.frame(.) %>%
        dplyr::select(Herd, Year, ncol(df2))  %>%
        rename(PrevLamb_Precip = ncol(.))
      # Bind them back together
      df <- rbind(df1, df2)
    }
    # Return df
    return(df)
  }))
})) %>% right_join(lamb, ., by = c("Herd", "Year"))

#--------------------------------#
# Postparturition Temperature ####
#--------------------------------#

# Lapply
lamb <- do.call(rbind, pblapply(1:length(year), function(i){
  do.call(rbind, lapply(1:length(sub), function(j){
    # Subset data to a year and subspecies
    df <- lamb %>% filter(Year == year[i] & Subspec == sub[j])
    
    # If subspec
    if(unique(df$Subspec) == "Desert"){
      # Paste the recovery months
      months <- paste0(year[i] - 1, 0,  3:5)
      
      # Subset the raster
      x <- tmean[[names(tmean) %in% months]]
      
      # Sum
      x <- app(x, fun = mean, na.rm = T)
      
      # Now extract
      df <- terra::extract(x, vect(df), mean, bind = T) %>%
        as.data.frame(.) %>%
        dplyr::select(Herd, Year, ncol(df))  %>%
        rename(PrevLamb_Temp = ncol(.))
    } else {
      # Get two dataframes with two timelines
      df1 <- df %>% filter(Timeline == time[1])
      df2 <- df %>% filter(Timeline == time[2])
      
      # Now create two separate months
      months1 <- c("06", "07", "08")
      months1 <- paste0(year[i] - 1, months1)
      months2 <- c("05", "06", "07")
      months2 <- paste0(year[i] - 1, months2)
      
      # Subset the raster
      x1 <- tmean[[names(tmean) %in% months1]]
      x2 <- tmean[[names(tmean) %in% months2]]
      
      # Sum
      x1 <- app(x1, fun = mean, na.rm = T)
      x2 <- app(x2, fun = mean, na.rm = T)
      
      # Now extract
      df1 <- terra::extract(x1, vect(df1), mean, bind = T) %>%
        as.data.frame(.) %>%
        dplyr::select(Herd, Year, ncol(df1))  %>%
        rename(PrevLamb_Temp = ncol(.))
      df2 <- terra::extract(x2, vect(df2), mean, bind = T) %>%
        as.data.frame(.) %>%
        dplyr::select(Herd, Year, ncol(df2))  %>%
        rename(PrevLamb_Temp = ncol(.))
      # Bind them back together
      df <- rbind(df1, df2)
    }
    # Return df
    return(df)
  }))
})) %>% right_join(lamb, ., by = c("Herd", "Year"))

##---------------------#
## Summer mean NDVI ####
##---------------------#
#
## ndvi year
#ndviyear <- year[2:length(year)]
#
## Lapply
#lamb <- do.call(rbind, pblapply(1:length(ndviyear), function(i){
#  
#  # Subset data to a year
#  df <- lamb %>% filter(Year == ndviyear[i])
#  
#  # Now get previous year
#  prev <- ndviyear[i] - 1
#  
#  # Paste summer months to it
#  months <- paste0(prev, 0, 5:9)
#  
#  # Subset the raster
#  x <- ndvi[[names(ndvi) %in% months]]
#  
#  # Sum
#  x <- app(x, fun = mean, na.rm = T)
#  
#  # Now extract
#  df <- terra::extract(x, vect(df), mean, bind = T) %>%
#    as.data.frame(.) %>%
#    dplyr::select(Herd, Year, ncol(df)) %>%
#    rename(Summer_NDVI = ncol(.)) %>%
#    mutate(Summer_NDVI = Summer_NDVI / 10000)
#})) %>%
#  right_join(lamb, ., by = c("Herd", "Year"))
#
#
##-------------------------#
## Previous Winter NDVI ####
##-------------------------#
#
## Lapply
#lamb <- do.call(rbind, pblapply(1:length(ndviyear), function(i){
#  
#  # Subset data to a year
#  df <- lamb %>% filter(Year == ndviyear[i])
#  
#  # Now get previous year
#  prev <- ndviyear[i] - 1
#  
#  # Paste winter months to it
#  months <- paste0(prev, 11:12)
#  
#  # Subset the raster
#  x1 <- ndvi[[names(ndvi) %in% months]]
#  
#  # Paste this year winter months
#  months <- paste0(ndviyear[i], 0, 1:2)
#  
#  # Subset
#  x2 <- ndvi[[names(ndvi) %in% months]]
#  
#  # Bind
#  x <- c(x1, x2)
#  
#  # Sum
#  x <- app(x, fun = mean, na.rm = T)
#  
#  # Now extract
#  df <- terra::extract(x, vect(df), mean, bind = T) %>%
#    as.data.frame(.) %>%
#    dplyr::select(Herd, Year, ncol(df)) %>%
#    rename(Winter_NDVI = ncol(.)) %>%
#    mutate(Winter_NDVI = Winter_NDVI / 10000)
#})) %>%
#  right_join(lamb, ., by = c("Herd", "Year"))
#
##--------------------------#
## Six months prior NDVI ####
##--------------------------#
#
## ndviym
#ndviym <- as.numeric(unique(lamb$YearMonth)) %>% sort()
#
## Lapply
#lamb <- do.call(rbind, pblapply(1:length(ndviym), function(i){
#  
#  # Subset data to a yearmonth
#  df <- lamb %>% filter(YearMonth == ndviym[i])
#  
#  # Get layer number of the yearmonth
#  num <- which(names(ndvi) == ndviym[i]) 
#  low <- as.numeric(num - 7)
#  high <- as.numeric(num - 1)
#  
#  if(low < 1){
#    low <- 1
#  }
#  
#  # Subset the raster
#  x <- subset(ndvi, c(low:high))
#  
#  # Sum
#  x <- app(x, fun = mean, na.rm = T)
#  
#  # Now extract
#  df <- terra::extract(x, vect(df), mean, bind = T) %>%
#    as.data.frame(.) %>%
#    dplyr::select(Herd, YearMonth, ncol(df)) %>%
#    rename(Month6_NDVI = ncol(.)) %>%
#    mutate(Month6_NDVI = Month6_NDVI / 10000)
#  
#  # Return 
#  return(df)
#})) %>%
#  right_join(lamb, ., by = c("Herd", "YearMonth"))

# Drop geometry
lamb <- st_drop_geometry(lamb)

# Recode months
lamb$Month <- as.numeric(lamb$MonthNum)

# Save the data for calculations
write_csv(lamb,
          "./Data/Clean/BHS_Modeling_Data.csv")
