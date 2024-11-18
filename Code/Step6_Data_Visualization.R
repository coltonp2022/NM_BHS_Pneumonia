# Data Visualization
# Colton Padilla
# 5/23/2023

# Last Updated: 5/23/2023

#-----------------#
# Housekeeping ####
#-----------------#

# Remove everything
rm(list = ls())

# Packages
library(tidyverse)
library(ggpubr)
library(ggpattern)

# Lamb Ewe Data
load("./Data/Clean/Lamb_Ewe_Plotting_Data.RData")

# Load Latir Population Size
pop_size <- read_csv("./Data/Clean/Population_Sizes.csv")

# Now join the population sizes to the lamb data 
pop_size <- full_join(pop_size,
                      lamb %>% select(Herd, Year, Lambs, Status, Month),
                      by = c("Herd", "Year", "Month")) %>%
  arrange(Herd, Year)

# Remove duplicated herd years and fix them
dup <- pop_size %>% 
  group_by(Herd, Year) %>% 
  filter(n() > 1) %>%
  summarise(Lambs = max(Lambs, na.rm = T),
            Est = max(Est, na.rm = T),
            Status = max(Status, na.rm = T))

# Remove and bind back
pop_size <- pop_size %>% 
  group_by(Herd, Year) %>%
  filter(n() == 1) %>%
  rbind(., dup)

# Write this out
write_csv(pop_size,
          "./Data/Clean/Lamb_Population_Size.csv")

# Create a tables and figures folder
if(!dir.exists("./TablesFigures/")){
  dir.create("./TablesFigures")
}

# Set Plotting Theme
theme_set(theme_bw())
theme_replace(axis.title.y = element_text(margin = margin(r = 5),
                                          angle = 90),
              axis.title.x = element_text(margin = margin(t = 10)),
              strip.text = element_text(margin = margin(3,0,3,0)),
              title = element_text(size = 12),
              axis.text = element_text(size = 12),
              axis.title = element_text(size = 14),
              legend.text = element_text(size = 12),
              plot.title = element_text(hjust = 0.5,
                                        vjust = 1))

# Recode Pre and Post
list <- list(ratios, full_ci, full_ci3, spec_ci, spec_ci3)
list2 <- lapply(list, function(x){
  x$Status <- ifelse(x$Status == "Pre", "Pre-exposure", "Post-exposure")
  return(x)
}) 
names(list2) <- c("ratios", "full_ci", "full_ci3", "spec_ci", "spec_ci3")
list2env(list2, .GlobalEnv)

# Same recode
ratios_lat$StatusLat <- ifelse(ratios_lat$StatusLat == "Pre",
                               "Pre-exposure",
                               "Post-exposure")

# Recode
lamb$Status <- ifelse(lamb$Status == "Pre",
                      "Pre-exposure",
                      "Post-exposure")

spec_ci %>% mutate_if(is.numeric, round, 0)

# Calculate pre and post-exposure difference within populations
lamb %>%
  group_by(Herd) %>%
  filter(n_distinct(Status) > 1) %>%
  group_by(Herd, Status, Subspecies) %>%
  summarise(mean = mean(Lambs)) %>%
  group_by(Herd, Subspecies) %>%
  filter(!(Herd %in% c("CUL", "SA"))) %>%
  summarise(perc_change = ((mean[Status == "Pre-exposure"] - 
                              mean[Status == "Post-exposure"]) / 
                             mean[Status == "Pre-exposure"]) * -100) %>%
  mutate(color = if_else(perc_change < 0, "Neg", "Pos"),
         Herd = if_else(Herd == "RGG", paste0(Herd, "*"), Herd)) %>%
  ggplot(.,
         aes(x = Herd,
             y = perc_change,
             fill = color)) +
  geom_bar(stat = "identity",
           color = "black") +
  scale_x_discrete(breaks = c("CAB", "FC", "LA", "LAT", "RGG*", "WH"),
    labels = c("Caballo", "Fra Cristobal", "Ladron",
                              "Latir", "Rio Grande Gorge*", "Wheeler Peak")) +
  facet_wrap(~ Subspecies, scales = "free_x") +
  scale_fill_manual(values = c("#E5140A", "gray70" )) +
  labs(y = "Percent Change in Mean Ratios Post-Exposure",
       x = NULL) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45,
                                   hjust = 1),
        axis.title.y = element_text(margin = margin(r = 10))) 

# Save
ggsave("./TablesFigures/Population_Perc_Decline.jpg",
       height = 7,
       width = 6,
       dpi = 600)

# Save
ggsave("./TablesFigures/Population_Perc_Decline_Paper.jpg",
       height = 7,
       width = 6.5,
       dpi = 600)

# Check ratios
lamb %>% 
  filter(Subspecies == "Rocky",
         Status == "Post-exposure") %>%
  ggplot(., 
         aes(Lambs,
             fill = Herd)) +
  geom_bar(position = "dodge")

# Linear regression of LAT < 2015
summary(lm(
  Lambs ~ Year,
  data = lamb %>% filter(Herd == "LAT", Year < 2016)
))

lamb$extra <- ifelse((lamb$Herd == "LAT" & lamb$Year < 2020), "1", "0")

parasite::boot2samp_t_test(lamb %>% filter(Herd == "LAT", Year >2015),
                           "Lambs", 
                           "extra")

# Plot
ggplot(lamb %>% filter(Herd == "LAT"),
       aes(x = Year,
           y = Lambs)) +
  geom_point() + 
  geom_smooth(method = "lm")

#----------------------------#
# Lamb-Ewe Ratio Plotting ####
#----------------------------#

# Set factor
spec_ci$Status <- factor(spec_ci$Status,
                         levels = c("Pre-exposure", "Post-exposure"))

lamb$Status <- factor(lamb$Status,
                      levels = c("Pre-exposure", "Post-exposure"))

# Change Subspecies
spec_ci$Subspecies <- gsub("Rocky", "Rocky Mountain", spec_ci$Subspecies)
lamb$Subspecies <- gsub("Rocky", "Rocky Mountain", lamb$Subspecies)


# Plot fulls
ggplot(full_ci,
       aes(x = Status, 
           y = Mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = Lower,
                    ymax = Upper),
                width = 0.1)

# Plot by species
CI <- ggplot(spec_ci,
       aes(x = Status, 
           y = Mean,
           fill = Status)) +
  geom_line(aes(group = 1)) +
  geom_errorbar(aes(ymin = Lower,
                    ymax = Upper),
                width = 0.2,
                position = "dodge") +
  geom_point(pch = 21,
             size = 5,
           color = "black") +
  labs(y = "Lambs per 100 ewes",
       color = "Exposure\nstatus",
       x = NULL)  +
  facet_wrap(~ Subspecies,
             scales = "free_x") +
  scale_fill_manual(values = c("gray70", "#E5140A")) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")

# RAW Geom Jitters
RAW <- ggplot(lamb,
       aes(x = Status,
           y = Lambs)) +
  geom_boxplot(aes(fill = Status),
               outlier.shape = NA) +
  geom_jitter(color = "black",
             position = position_jitter(width = 0.15),
             alpha = 0.3)  +
  facet_wrap(~Subspecies,
             scales = "free_x") +
  labs(y = "Lambs per 100 ewes",
       fill = "",
       x = NULL) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_fill_manual(values = c("gray70", "#E5140A"))

# RAW Density ridges
ggplot(lamb,
       aes(x = Lambs)) +
  geom_density(aes(x = Lambs,
                   fill = Status),
               color = "black",
               alpha = 0.6)  +
  facet_wrap(~Subspecies,
             scales = "free") +
  labs(x = "Lambs per 100 ewes",
       fill = "",
       y = "") +
  scale_fill_manual(values = c("gray70","#E5140A")) +
  theme(legend.position = "bottom")

# GGarrange
ggarrange(RAW, CI, ncol = 1, labels = c("A.", "B."),
          common.legend = T, legend = "bottom")

# Save
ggsave("./TablesFigures/Lamb_Ewe_CI_Raw.jpg",
       dpi = 600,
       height = 7,
       width = 5)

# Arrange for paper
ggarrange(RAW + rremove("ylab") + labs(title = "Observed Values"),
          CI+ rremove("ylab") + labs(title = "Confidence Intervals"), 
          ncol = 2, 
          labels = c("A.", "B."),
          common.legend = T, legend = "bottom") %>%
  annotate_figure(., left = text_grob("Lambs per 100 ewes", rot = 90))

# Save
ggsave("./TablesFigures/Lamb_Ewe_CI_Raw_Paper.jpg",
       dpi = 600,
       height = 5,
       width = 8.5)

# Plot pre-post ratios
ratios$Status <- factor(ratios$Status,
                        levels = c("Pre-exposure", "Post-exposure"))

# Code subspecies
rocky <- c("PEC", "WH", "SFR", "TC", "DC", "MZ", "JZ", "RGG", "LAT")

# If else
ratios$Subspecies <- ifelse(ratios$Herd %in% rocky,
                          "Rocky Mountain",
                          "Desert")

# Plot 
ggplot(ratios %>% filter(n() > 1),
       aes(x = Herd,
           y = Ratio,
           fill = Status)) +
  geom_bar(stat = "Identity",
           position = position_dodge2(width = 0, preserve = "single"),
           color = "black") +
  scale_fill_manual(values = c("gray70", "#E5140A")) +
  scale_color_manual(values = c("gray70", "#E5140A")) +
  labs(y = "Lambs per 100 ewes",
       fill = NULL,
       x = NULL) +
  facet_wrap(~Subspecies,
             scales = "free_x") +
  theme(legend.position = "bottom")

# Save
ggsave("./TablesFigures/Pop_Ratios.jpg",
       dpi = 600,
       height = 4,
       width = 6.5)

# Save for presentation
ggsave("./TablesFigures/Pop_Ratios_Presentation.jpg",
       dpi = 600,
       height = 7,
       width = 6)

#-----------------------------------------------#
# Plot RAW RATIOS of all pops with Estimates ####
#-----------------------------------------------#

# Add in exposure dates
pop_lamb_est$expdate <- ifelse(pop_lamb_est$Herd == "CAB", 2019,
                               ifelse(pop_lamb_est$Herd == "FC", 2017,
                                      ifelse(pop_lamb_est$Herd == "LA", 2017,
                                             ifelse(pop_lamb_est$Herd == "RGG", 2020,
                                                    ifelse(pop_lamb_est$Herd == "SA", 2005,
                                                           ifelse(pop_lamb_est$Herd == "WH", 2022, 2022))))))

# Change the name
pop_lamb_est$Herd <- ifelse(pop_lamb_est$Herd == "CAB", "Caballo",
                               ifelse(pop_lamb_est$Herd == "FC", "Fra Cristobal",
                                      ifelse(pop_lamb_est$Herd == "LA", "Ladron",
                                             ifelse(pop_lamb_est$Herd == "RGG", "Rio Grande Gorge",
                                                    ifelse(pop_lamb_est$Herd == "SA", "San Andres",
                                                           ifelse(pop_lamb_est$Herd == "WH", "Wheeler Peak", "Latir"))))))

# Now status
pop_lamb_est$Status <- ifelse(pop_lamb_est$Year < pop_lamb_est$expdate,
                              "Pre-exposure", 
                              "Post-exposure")

# Factor the status
pop_lamb_est$Status <- factor(pop_lamb_est$Status,
                              levels = c("Pre-exposure", "Post-exposure"))

# Unique herds
herd <- unique(pop_lamb_est$Herd)


# Plot function to apply

  # Plot
ggplot(pop_lamb_est,
       aes(x = Year,
           y = Lambs)) +
  geom_vline(aes(xintercept = expdate),
             size = 1,
             linetype = "dashed",
             color = "#E5140A") +
  geom_bar(aes(fill = Status),
           stat = "identity",
           color = "black",
           linewidth = 0.6) +
  geom_line(data = pop_lamb_est %>% drop_na(Est2),
            aes(x = Year,
                y = Est2 / 5),
            linewidth = 0.8) +
  geom_point(aes(x = Year,
                 y = Est2 / 5),
             size = 3)  +
  theme(legend.position = "bottom",
        axis.title.y.right = element_text(margin = margin(l = 10)),
        panel.spacing = unit(0.5, "lines")) +
  scale_y_continuous(name = "Lambs per 100 ewes",
                     limits = c(0, 80),
                     sec.axis = sec_axis(trans = ~.* 5, name = "Population Estimate")) +
  scale_fill_manual(
    breaks = c("Pre-exposure", "Post-exposure"),
    values = c("gray70", "#E5140A")) +
  labs(x = "Year") +
  xlim(2000, 2023) +
  facet_wrap(~Herd) 

# Save 
ggsave("./TablesFigures/Pop_Lamb_Estimates.jpg",
       dpi = 600,
       height = 10,
       width = 12)

#-----------------------------------------------#
# Plot the Latir Populations lamb-ewe ratios ####
#-----------------------------------------------#

# Add coloring 
LAT <- pop_size %>%
  filter(Herd == "LAT") %>%
  mutate(Exp = case_when(
    Year < 2016 ~ "Unknown",
    Year > 2015 & Year < 2022 ~ "Potential epizootic",
    Year > 2021 ~ "Positive"
  ))

# Pattern
LAT$pattern <- ifelse((LAT$Year > 2015 & LAT$Year < 2022),
                      "stripe",
                      "none")

# Change name
LAT$Herd <- "Latir"

# Set factor
LAT$Exp <- factor(LAT$Exp, 
                  levels = c("Unknown", "Potential epizootic",
                             "Positive"))

LAT <- LAT %>% drop_na(Lambs) %>% distinct()

# Plot
ggplot(LAT,
       aes(x = Year,
           y = Lambs,
           fill = Exp,
           pattern = Exp)) +
  geom_bar_pattern(stat = "identity",
    color = "black",
    linewidth = 0.6, 
    pattern_fill = "black",
    pattern_angle = 45,
    pattern_density = 0.1) +
  geom_line(aes(x = Year,
                y = Est / 2.2,
                color = "Population Est."),
            inherit.aes = F,
            linewidth = 0.8) +
  geom_point(aes(x = Year,
                 y = Est / 2.2,
                 color = "Population Est."),
             inherit.aes = F,
             size = 3) +
  labs(y = "Lambs per 100 ewes",
       fill = "",
       x = NULL,
       pattern = "") +
  scale_pattern_manual(values = c(Unknown = "none", `Potential epizootic` = "stripe",
                                  Positive = "none"),
                       guide = "none")+
  scale_fill_manual(values = c("grey70", "#E5140A", "#E5140A")) +
  guides(fill = guide_legend(override.aes = list(pattern = c("none", "stripe",
                                                             "none")))) +
  theme(legend.position = "bottom",
        axis.title.y.right = element_text(margin = margin(l = 10))) +
  facet_wrap(~ Herd) +
  scale_y_continuous(name = "Lambs per 100 ewes",
                     sec.axis = sec_axis(trans = ~. * 2.2, name = "Population Estimate")) +
  scale_color_manual(name = "", values = c("Population Est." = "black"))

# Save
ggsave("./TablesFigures/Latir_Ratios.jpg",
       width = 7,
       height = 6,
       dpi = 600)

# fdsafda
ggsave("./TablesFigures/Latir_Ratios_presentation.jpg",
       width = 6,
       height = 7,
       dpi = 600)

#---------------------------------------------------#
# View lamb ewe ratios for post-exposure < Upper ####
#---------------------------------------------------#
      
# Code 
lamb$exp <- ifelse(lamb$Subspecies == "Rocky Mountain" & lamb$Lambs < 34,
                   1,
                   ifelse(lamb$Subspecies == "Desert" & lamb$Lambs < 31,
                          1, 0))

lamb$exp2 <- ifelse(lamb$Subspecies == "Rocky Mountain" & lamb$Lambs > 34,
                    1,
                    ifelse(lamb$Subspecies == "Desert" & lamb$Lambs > 31,
                           1, 0))

# Calulate percent below upper for post-exposure
lamb %>%
  group_by(Herd, StatusLat) %>%
  summarise(years = sum(exp),
            years2 = sum(exp2),
            perc = years / n(),
            perc2 = years2 / n()) %>%
  ggplot(., 
         aes(x = Herd,
             y = perc2,
             fill = StatusLat)) +
  geom_point(pch = 21,
             size = 4)
                   
# Plot
ggplot(lamb,
       aes(x = Year,
           y = Lambs,
           fill = as.character(exp))) +
  geom_bar(stat = "identity",
           color = "black") +
  scale_fill_manual(values = c("red", "grey70")) +
  labs(x = NULL,
       fill = NULL) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "bottom") +
  facet_wrap(~ Herd,
             ncol = 4) 

# Save 
ggsave("./TablesFigures/All_Raw_Ratios.jpg",
       width = 10,
       height = 14,
       dpi = 600)

#------------------#
# Density plots ####
#------------------#

# Total sample size
samp <- lamb %>% 
  dplyr::group_by(Subspecies, Status) %>%
  summarise(Samp_Size = paste0("n = ", n()))

# Plot raw ratios
ggplot(lamb %>% filter(Herd %in% c("CAB", "FC", "LA", "RGG")),
       aes(x = Lambs,
           fill = Status)) +
  geom_density(color = "black",
               alpha = 0.4) + 
  #geom_label(aes(y = 0.3,
  #               label = Samp_Size)) +
  facet_wrap(~ Herd,
             scales = "free") +
  labs(y = "") +
  theme(legend.position = "bottom",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

# Try a boxplot
box <- ggplot(lamb) +
  geom_boxplot(aes(y = Lambs,
                   x = Status,
                   fill = Status),
              color = "black")  +
  geom_text(data = samp,
             aes(fill = Status,
                 x = Status,
                 y = lamb %>% dplyr::group_by(Subspecies, Status) %>% 
                   summarise(med = median(Lambs)) %>%
                   pull() + 2.25,
                 label = Samp_Size)) +
  facet_wrap(~Subspecies) +
  labs(y = "Lambs per 100 ewes",
       fill = "",
       x = NULL) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") +
  scale_fill_manual(values = c("gray70", "#E5140A")) 

# Save the boxplot
ggsave("./TablesFigures/Raw_Lamb_Ewe_Boxplot.jpg",
       height = 7,
       width = 6,
       dpi = 600)

# Latir change
lamb$Status <- ifelse(lamb$StatusLat == "Post",
                      "Post-exposure",
                      lamb$Status)

ggplot(lamb) +
  geom_boxplot(aes(y = Lambs,
                   x = Status,
                   fill = Status),
               color = "black")  +
  facet_wrap(~Subspecies) +
  labs(y = "Lambs per 100 ewes",
       fill = "",
       x = NULL) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") +
  scale_fill_manual(values = c("gray70", "#E5140A")) 


