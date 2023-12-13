# Map Creation
# Colton Padilla
# 6/2/2023

# Last Updated: 6/2/2023

#-----------------#
# Housekeeping ####
#-----------------#

# Packages
library(tidyverse) # Piping
library(sf) # Shapefiles
library(basemaps) # Basemaps
library(ggnewscale) # Set new fill scales after gg_raster
library(ggrepel) # population labels
library(cowplot)
library(raster)
library(ggspatial)
library(ggpattern)
library(ggpubr)

# Set Plotting Theme
theme_set(theme_bw())
theme_replace(axis.title.y = element_text(margin = margin(r = 10),
                                          angle = 90),
              axis.title.x = element_text(margin = margin(t = 10)),
              strip.text = element_text(margin = margin(3,0,3,0)),
              title = element_text(size = 12),
              axis.text = element_text(size = 12),
              axis.title = element_text(size = 14),
              legend.text = element_text(size = 18),
              legend.title = element_text(size = 20))

# State Shapefile
nm <- st_read("./Data/GIS_Files/cb_2018_us_state_500k.shp") %>%
  filter(NAME == "New Mexico") %>%
  st_transform(., 3857) %>%
  st_union() 

# Check projection
st_crs(nm)

# Load the BHS populations
pop <- st_read("./Data/GIS_Files/Bighorn_Sheep_Ranges.shp") %>%
  st_transform(., 3857) %>%
  dplyr::select(NAME, SPECIES, STATUS, geometry) %>%
  rename(Subspecies = SPECIES)

#-------------------------------#
# Manipulate Population Data ####
#-------------------------------#

# If else
pop$Subspecies <- ifelse(pop$Subspecies == "ROCKY MOUNTAIN",
                         "Rocky Mountain",
                         "Desert")

# Pull out San Andres and remove the historical population
San <- pop[grepl("San Andres", pop$NAME),] %>%
  summarise(geometry = st_union(geometry)) %>%
  mutate(NAME = "San Andres",
         Subspecies = "Desert") 

# Pull out RGG and summarise it 
RGG <- pop[grepl("Rio Grande", pop$NAME),]  %>%
  summarise(geometry = st_union(geometry)) %>%
  mutate(NAME = "Rio Grande Gorge",
         Subspecies = "Rocky Mountain") 

# Now remove those and then bind back
pop <- pop[!(grepl("San Andres", pop$NAME) | grepl("Rio Grande", pop$NAME)),] %>%
  dplyr::select(-STATUS) %>%
  rbind(., San) %>%
  rbind(., RGG)

# Toss Alamo Hueco
pop <- pop[!grepl("Alamo", pop$NAME),] 

# Now remove extra words 
pop$NAME <- gsub(" Mtns", "", pop$NAME) %>%
  gsub(" Range", "", .) %>%
  gsub(" Sand Canyon", "", .) 

# Now add a column with abbreviations
pop$Abrev <- c("DC", "FC", "HAT", "LA", "LAT", "MZ", "PEC", "PE", "SFR",
               "TC", "WH", "CAB", "CUL", "JZ", "SAC", "SA", "RGG")

# Now change the Hatchet name
pop$NAME <- gsub("Hatchet", "Hatchets", pop$NAME)

# Now save the new shapefile
st_write(pop,
         dsn = "./Data/GIS_Files",
         layer = "NM_populations",
         driver = "Esri Shapefile",
         append = F)

# Check distances between SA and others
SA <- pop %>% filter(NAME == "San Andres")
FC <- pop %>% filter(NAME == "Fra Cristobal")
Cab <- pop %>% filter(NAME == "Caballo")

st_distance(SA, FC, by_element = T)
st_distance(SA, Cab, by_element = T)

# Add positives and negatives to population map
pop$NAME <- ifelse(pop$Abrev %in% c("CUL", "RGG", "LAT", "WH", "LA", "FC",
                                    "CAB", "SA", "SAC"),
                   paste0(pop$NAME, " (+)"),
                   paste0(pop$NAME, " (-)"))

#-------------------#
# Create the map ####
#-------------------#

# Get basemap
base <- basemap_raster(nm, "esri", "world_terrain_base")

# Now mask the basemap
base <- mask(base, as_Spatial(nm))
rast <- raster::raster(base)
# Get bounding box of NM
NM_box <- st_bbox(nm)

# Create manual nudges for labels
  # The way I label this has an argument to nudge them a certain distance from 
  # the default. I had to do it manually by changing these distances.
pop$NudgeX <- c(-60000, # DC
                10000, # FC
                50000, # HAT
                20000, # LA
                70000, # LAT
                50000, # MZ
                10000, # PEC
                50000, # PE
                60000, # SFR
                60000, # TC
                70000, # WH
                -30000, # CAB
                40000, # CUL
                -50000, # JZ
                60000, # SAC
                70000, # SA
                -50000) # RGG

# Y distance
pop$NudgeY <- c(-20000, # DC
                30000, # FC
                30000, # HAT
                40000, # LA
                -20000, # LAT
                0, # MZ
                -30000, #PEC 
                30000, # PE
                40000, # SFR
                20000, # TC
                -20000, # WH
                -40000, # CAB
                -10000, # CUL
                0, # JZ
                0, # SAC
                10000, # SA
                0) # RGG

# plot the map
gg_raster(base) + # The original raster
  new_scale_fill() + # Reset the fill scale for the rest of the plot
   # Set the outline weight
  geom_sf(data = pop, # Bring in the pop shapefiles
          aes(fill = Subspecies,
              color = Subspecies), # Set outline color
          linewidth = 0.5) +
  geom_sf(data = nm, # NM outline
          fill = NA, # Remove the fill
          linewidth = 1) +
  geom_label_repel(data = pop, # Now for the labeling, same shapefile
                  aes(label = NAME, # Use the name as the label
                      geometry = geometry), # Set the geometry
                  stat = "sf_coordinates", # Found this online and it works??
                  nudge_x = pop$NudgeX, # Now for nudge X distances
                  nudge_y = pop$NudgeY, # Nudge Y Distance
                  min.segment.length = 0) + # Always have a line segment
  coord_sf(expand = F,
           xlim = c(NM_box[1] - 10000, NM_box[3] + 10000),
           ylim = c(NM_box[2] - 10000, NM_box[4] +10000)) +
  labs(y = NULL, # Remove Y coordinates
       x = NULL, # Remove X coordinates
       fill = "Subspecies",
       color = "Subspecies") + # Reset the legend name
  scale_fill_manual(values = alpha(c("#E69F00", "#56B4E9"), 0.6)) + # Colorblind friendly
  scale_color_manual(values = alpha(c("#E69F00", "#56B4E9"), 0.8)) +
  annotation_north_arrow(style = north_arrow_fancy_orienteering, # North arrow
                         location = "br", # Bottom right
                         pad_y = unit(0.6, "cm"), # How far from the bottom
                         pad_x = unit(1.2, "cm")) + # How far from the right
  annotation_scale(location = "br", # Bottom Right
                   pad_x = unit(3.2, "cm"), # How far from bottom
                   pad_y = unit(1, "cm")) + # How far from right
  theme_bw() + # Get rid of all lines outside map
  theme(legend.position = c(0.4, 0.05), # Set the legend position
        legend.box = "horizontal",
        panel.grid = element_blank()) + # Make it horizontal
  guides(fill = guide_legend(title.position = "left", # also help with horizontal
                             nrow = 1)) # Make sure its a single row

# Now save the map
ggsave("./TablesFigures/Study_Area_Map.jpg",
       dpi = 600,
       height = 9,
       width = 8)

#----------------------------------#
# Create map of only San Andres ####
#----------------------------------#

# Filter pops to just San Andres
SA <- pop %>% filter(NAME == "San Andres")

# Box for plot
plot_box <- st_as_sfc(st_bbox(st_buffer(SA, 20000) %>% 
                                st_buffer(., 20000, singleSide = T))) 

# Box
box <- st_bbox(st_buffer(SA, 20000) %>% 
                 st_buffer(., 20000, singleSide = T))

# Inset base
# Get basemap
inset_base <- basemap_raster(plot_box, "esri", "world_terrain_base")

# Now mask the basemap
inset_base <- mask(inset_base, as_Spatial(plot_box)) %>%
  mask(., as_Spatial(nm))

# plot the map
main <- gg_raster(base) + # The original raster
  new_scale_fill() + # Reset the fill scale for the rest of the plot
  geom_sf(data = nm, # NM outline
          fill = NA, # Remove the fill
          linewidth = 1) + # Set the outline weight
  geom_sf(data = SA, # Bring in the pop shapefiles
          aes(fill = Subspecies), # Fill them by subspecies
          color = "black", # Set outline color
          linewidth = 0.5) + # Outline weight
  geom_sf(data = plot_box,
          fill = NA,
          color = "red",
          linewidth = 1) +
  labs(y = NULL, # Remove Y coordinates
       x = NULL, # Remove X coordinates
       fill = "Subspecies") + # Reset the legend name
  scale_fill_manual(values = alpha(c("#F0E442", "#0072B2"), 0.6)) + # How far from right
  theme_test() +
  theme(legend.position = "none", # Set the legend position
        legend.box = "horizontal",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.ticks.length = unit(0, "pt"),
        plot.margin = margin(0,0,0,0, "cm"),
        panel.background = element_blank()) + # Make it horizontal
  guides(fill = guide_legend(title.position = "left", # also help with horizontal
                             nrow = 1)) # Make sure its a single row

# Draw SA map
inset <- gg_raster(inset_base) + # The original raster
  new_scale_fill() + # Reset the fill scale for the rest of the plot
  coord_sf(expand = F) +
  geom_sf(data = nm, # NM outline
          fill = NA, # Remove the fill
          linewidth = 1) + # Set the outline weight
  geom_sf(data = SA, # Bring in the pop shapefiles
          aes(fill = Subspecies), # Fill them by subspecies
          color = "black", # Set outline color
          linewidth = 0.5) + # Outline weight
  geom_sf(data = plot_box,
          fill = NA,
          color = "black",
          linewidth = 1) +
  coord_sf(expand = F,
           xlim = c(box[1], box[3]),
           ylim = c(box[2], box[4])) + 
  labs(y = NULL, # Remove Y coordinates
       x = NULL, # Remove X coordinates
       fill = "Subspecies") + # Reset the legend name
  scale_fill_manual(values = alpha(c("#F0E442", "#0072B2"), 0.4)) + # How far from right
  theme_bw() +
  theme(legend.position =  "none", # Set the legend position
        legend.box = "horizontal",
        panel.background = element_blank()) + # Make it horizontal
  guides(fill = guide_legend(title.position = "left", # also help with horizontal
                             nrow = 1)) # Make sure its a single row

# Draw the inset
ggdraw() +
  draw_plot(inset) +
  draw_plot(main,
            height = 0.35,
            width = 0.35,
            x = 0.134,
            y = 0.7055) 

# Save
ggsave("./TablesFigures/San_Andres_Map.jpg",
       width = 5,
       height = 9,
       dpi = 600)

#-------------------------------------#
# Get an exposed and unexposed map ####
#-------------------------------------#

# Add exposed category
pop$exp <- c("Negative", "Positive", "Negative", "Positive", "Positive",
             "Negative", "Negative", "Negative", "Negative", "Negative",
             "Positive", "Positive", "Positive", "Negative", "Positive",
             "Positive", "Positive")

# Keep labels of used populations
pop$NAME <- ifelse(!(pop$Abrev %in% c("CAB", "LAT", "FC", "LA", "RGG", "WH")),
                   NA,
                   pop$NAME)

# Change latir nudge
pop$NudgeY <- ifelse(pop$Abrev == "LAT",
                     0,
                     pop$NudgeY)

pop$NudgeX <- ifelse(pop$Abrev == "LAT",
                     90000,
                     pop$NudgeX) 

# plot the map
infect <- gg_raster(base) + # The original raster
  new_scale_fill() + # Reset the fill scale for the rest of the plot
  geom_sf(data = nm, # NM outline
          fill = NA, # Remove the fill
          linewidth = 1) + # Set the outline weight
  geom_sf(data = pop, # Bring in the pop shapefiles
          aes(fill = exp,
              color = exp), # Set outline color
          linewidth = 0.5) + # Outline weight 
  geom_text_repel(data = pop, # Now for the labeling, same shapefile
                  aes(label = NAME, # Use the name as the label
                      geometry = geometry), # Set the geometry
                  stat = "sf_coordinates", # Found this online and it works??
                  nudge_x = pop$NudgeX, # Now for nudge X distances
                  nudge_y = pop$NudgeY, # Nudge Y Distance
                  min.segment.length = 0,
                  size = 5) +
  coord_sf(expand = F,
           xlim = c(NM_box[1] - 10000, NM_box[3] + 10000),
           ylim = c(NM_box[2] - 10000, NM_box[4] +10000)) + # Always have a line segment
  labs(y = NULL, # Remove Y coordinates
       x = NULL, # Remove X coordinates
       fill = "Exposure",
       color = "Exposure") + # Reset the legend name
  scale_fill_manual(values = alpha(c("grey70", "#E5140A"), 0.4)) + # How far from right
  scale_color_manual(values = c("black", "#E5140A"))+
  theme_void() +
  theme(legend.position = c(0.4, 0.05), # Set the legend position
        legend.box = "horizontal") + # Make it horizontal
  guides(fill = guide_legend(title.position = "left", # also help with horizontal
                             nrow = 2)) + # Get rid of all lines outside map
  annotation_north_arrow(style = north_arrow_fancy_orienteering, # North arrow
                         location = "br", # Bottom right
                         pad_y = unit(0.6, "cm"), # How far from the bottom
                         pad_x = unit(1, "cm")) + # How far from the right
  annotation_scale(location = "br", # Bottom Right
                   pad_x = unit(3.5, "cm"), # How far from bottom
                   pad_y = unit(1, "cm")) +
  theme(legend.text = element_text(size = 16),
        legend.title = element_text(size = 18))

# View
infect

# Save
ggsave("./TablesFigures/Exposure_Status.jpg",
       width = 8,
       height = 9,
       dpi = 600)
#--------------------#
# Publication Map ####
#--------------------#

# Scale fill identity
pop$Subspecies <- ifelse(pop$exp == "Positive",
                         paste0(pop$Subspecies, " (Positive)"),
                         paste0(pop$Subspecies, " (Negative)"))

# plot the map
gg_raster(base) + # The original raster
  new_scale_fill() + # Reset the fill scale for the rest of the plot
  geom_sf(data = nm, # NM outline
          fill = NA, # Remove the fill
          linewidth = 1) + # Set the outline weight
  geom_sf(data = pop, # Bring in the pop shapefiles
          aes(fill = Subspecies,
              color = Subspecies),# Set outline color
          linewidth = 0.5) + # Always have a line segment
  coord_sf(expand = F,
           xlim = c(NM_box[1] - 10000, NM_box[3] + 10000),
           ylim = c(NM_box[2] - 10000, NM_box[4] +10000)) +
  labs(y = NULL, # Remove Y coordinates
       x = NULL, # Remove X coordinates
       fill = "Subspecies",
       color = "Subspecies") + # Reset the legend name
  scale_fill_manual(
    name = "Subspecies (Infection Status)",
    values = alpha(c("#E69F00","#E69F00" , "#56B4E9", "#56B4E9"), 0.4),
    labels = c("Desert (Negative)", "Desert (Positive)",
               "Rocky Mountain (Negative)", "Rocky Mountain (Positive)")) + # Colorblind friendly
  scale_color_manual(
    name = "Subspecies (Infection Status)",
    values = c("black", "#E5140A", "black", "#E5140A"),
    labels = c("Desert (Negative)", "Desert (Positive)",
               "Rocky Mountain (Negative)", "Rocky Mountain (Positive)")) +
  annotation_north_arrow(style = north_arrow_fancy_orienteering, # North arrow
                         location = "br", # Bottom right
                         pad_y = unit(0.5, "cm"), # How far from the bottom
                         pad_x = unit(0.4, "cm")) + # How far from right
  theme_bw() +
  theme(legend.position =  c(0.55, 0.05), # Set the legend position
        legend.box = element_blank(),
        legend.background = element_blank(),
        panel.grid = element_blank()) + # Make it horizontal
  guides(fill = guide_legend(title.position = "left", # also help with horizontal
                             nrow = 2))# Make sure its a single row

# Save
ggsave("./TablesFigures/Subspecies_Exposure.jpg",
       width = 8,
       height = 9,
       dpi = 600)
