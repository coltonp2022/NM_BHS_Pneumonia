# Create fancy map of precip data

# Packages
library(tidyverse)
library(sf)
library(terra)
library(pbapply)
library(tidyterra)
library(ggspatial)

# Load in the PRISM data
prism <- rast("./Data/GIS_Files/1999_2022_monthly_precip_NM.tif")

# Subset to 2022 August
r <- subset(prism, "202208") %>%
  project(., "epsg:4326", threads = T)

# State Shapefile
nm <- st_read("./Data/GIS_Files/cb_2018_us_state_500k.shp") %>%
  filter(NAME == "New Mexico") %>%
  st_transform(., 4326) %>%
  #st_union() %>%
  st_as_sf() 

# Crop
r <- crop(mask(r, vect(nm)), nm) 
r <- na.omit(r, field = ",", geom = T)
is.na(r)
# Populations
pop <- st_read("./Data/GIS_Files/NM_populations_v2_1.shp") %>%
  st_transform(4326)


# Plot
ggplot() + 
  geom_spatraster(data = r,
                  na.rm = T) +
  geom_sf(data = pop,
          fill = NA,
          aes(color = Subspecies),
          linewidth = 1) +
  geom_sf(data = nm,
          fill = NA,
          linewidth = 1) +
  scale_fill_continuous(na.value = "transparent") +
  annotation_north_arrow(style = north_arrow_fancy_orienteering, # North arrow
                         location = "br", # Bottom right
                         pad_y = unit(0.6, "cm"), # How far from the bottom
                         pad_x = unit(1.2, "cm")) + # How far from the right
  annotation_scale(location = "br", # Bottom Right
                   pad_x = unit(3.2, "cm"), # How far from bottom
                   pad_y = unit(1, "cm")) + # How far from right
  theme_bw() + # Get rid of all lines outside map
  theme(legend.position = c(0.4, 0.05),
        panel.grid = element_blank()) +
  guides(fill = "none",
         color = guide_legend(title.position = "left", # also help with horizontal
                              nrow = 1))

# Now save the map
ggsave("./TablesFigures/Precip_Map.jpg",
       dpi = 600,
       height = 9,
       width = 8)
