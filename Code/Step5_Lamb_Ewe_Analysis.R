# Analysis of Lamb-Ewe Ratios
# Colton Padilla
# 5/23/2023

# Last Updated: 5/23/2023

#-----------------#
# Housekeeping ####
#-----------------#

# Packages
library(tidyverse)
library(parasite)

# Load the data
lamb <- read_csv("./Data/Clean/Lamb-Ewe_Ratios_Clean.csv") %>%
  arrange(Herd, Year)

# Get counts of timeframes by subspecies
lamb %>% dplyr::group_by(Subspecies) %>% count(Month) %>% print(n = 100)

# Choose the later observation
lamb <- lamb %>% 
  group_by(Herd, Year) %>%
  slice(min(n(), 2)) %>%
  ungroup() 

# Remove the Hatchets 2005
lamb <- lamb %>%
  filter(!(Herd == "HAT" & Year == 2005))

# Look for number of herds
lamb %>% count(Herd, Status) %>% filter(Status == "Post") %>% mutate(mean = mean(n))

# Get the longest stretches for each post-exposure population
years <- lamb %>% 
  group_by(Herd) %>%
  filter(Status == "Post") %>%
  summarise(range = max(Year) - min(Year))

years2 <- lamb %>% 
  group_by(Herd) %>%
  filter(Status3 == "Post") %>%
  summarise(range = max(Year) - min(Year))

# Plot ratios < 20
rats <- lamb %>% filter(Lambs < 21)
ggplot(rats, aes(x = Year, y = Lambs)) +
  geom_bar(stat = "identity") + facet_wrap(~ Herd)

# Calculate means by year and get percent declines
lamb %>% 
  group_by(Status, Year, Subspecies) %>%
  summarise(mean = mean(Lambs)) %>% 
  group_by(Year, Subspecies) %>%
  filter(n() > 1) %>%
  arrange(Year) %>%
  summarise(perc = ((mean[2] - mean[1]) / mean[2]) * 100) %>%
  ungroup() %>%
  mutate(med = median(perc))

# Calculate yearly post-exposure stuff
lamb %>% 
  group_by(Year, Subspecies) %>%
  filter(!(n_distinct(Status) == 1)) %>%
  ggplot(.,
         aes(x = Herd,
             y = Lambs,
             fill = Status)) + 
  geom_bar(stat = "identity") +
  facet_wrap(~ Year + Subspecies,
             scales = "free")

# Set levels
lamb$Month <- factor(lamb$Month,
                     levels = c("JAN/AUG", "APR", "MAY",
                                "JUN", "AUG", "SEP", "OCT",
                                "OCT/NOV", "DEC", NA))

lamb %>% 
  filter(Subspecies == "Desert") %>%
  drop_na(Month) %>%
  mutate(ind = if_else(Herd == "SA", "Y", "N")) %>%
  ggplot(.,
         aes(x = Month,
             fill = ind)) +
  geom_histogram(stat = "count",
                 position = "dodge")

# Check out the lamb:ewe ratios in San Andres
lamb %>%
  filter(Herd == "SA") %>%
  group_by(Status) %>%
  summarise(mean = mean(Lambs),
            low = min(Lambs),
            upp = max(Lambs))

#------------------#
# Full Analysis ####
#------------------#

# Plain confidence intervals
full_ci <- bcaCI(lamb,
                 "Lambs",
                 "Status")

# Bootstrap test
full_ci$P_Value <- boot2samp_t_test(lamb,
                 "Lambs",
                 "Status",
                 "greater")

# 3 year confidence
full_ci3 <- bcaCI(lamb,
                  "Lambs",
                  "Status3")

# Bootstrap test
full_ci3$P_Value <- boot2samp_t_test(lamb,
                 "Lambs",
                 "Status3",
                 "greater")

#------------------#
# By Subspecies ####
#------------------#

# Test for normality
shapiro.test(lamb %>% filter(Subspecies == "Desert",
                             Status == "Post") %>% pull(Lambs))
shapiro.test(lamb %>% filter(Subspecies == "Rocky",
                             Status == "Post") %>% pull(Lambs))

# Species
spec <- unique(lamb$Subspecies)

# Loop
spec_ci <- do.call(rbind,lapply(1:length(spec), function(i){
  dat <- lamb %>%
    filter(Subspecies == spec[i]) 
  
  bcaCI(dat,
        "Lambs",
        "Status") %>%
    mutate(Subspecies = spec[i])
}))

# Loop for 3
spec_ci3 <- do.call(rbind,lapply(1:length(spec), function(i){
  dat <- lamb %>%
    filter(Subspecies == spec[i])
  
  bcaCI(dat,
        "Lambs",
        "Status3") %>%
    mutate(Subspecies = spec[i])
}))

# Test for rocky
spec_ci$P_Value[1:2] <- boot2samp_t_test(lamb %>% filter(Subspecies == "Rocky"),
                 "Lambs",
                 "Status",
                 "greater")

# Test for desert
spec_ci$P_Value[3:4] <- boot2samp_t_test(lamb %>% filter(Subspecies == "Desert"),
                 "Lambs",
                 "Status", 
                 "greater")

# Test for rocky
spec_ci3$P_Value[1:2] <- boot2samp_t_test(lamb %>% filter(Subspecies == "Rocky"),
                 "Lambs",
                 "Status3")

# Test for desert
spec_ci3$P_Value[3:4] <- boot2samp_t_test(lamb %>% filter(Subspecies == "Desert"),
                 "Lambs",
                 "Status3")

# Rocky by southern vs northern
boot2samp_t_test(lamb %>% filter(Subspecies == "Rocky", Status == "Pre") %>%
                   mutate(south = if_else(Herd %in% c("SFR", "TC"), "Y", "N")),
                 "Lambs",
                 "south")

bcaCI(lamb %>% filter(Subspecies == "Rocky", Status == "Pre") %>%
        mutate(south = if_else(Herd %in% c("SFR", "TC"), "Y", "N")),
      "Lambs",
      "south")

# GLM
mod <- glmmTMB::glmmTMB(Lambs ~ Status + (1|Herd),
           data = lamb,# %>% filter(Lambs != 0),
           family = glmmTMB::nbinom2)

sjPlot::plot_model(mod, "re")

#-------------------------------------------------------------#
# Test for differences pre- and post exposure in each herd ####
#-------------------------------------------------------------#

# Get herds with both pre and post
herds <- lamb %>%
  group_by(Herd) %>%
  filter(!all(Status == "Pre"),
           !all(Status == "Post")) %>%
  ungroup()

# Get the unique herds
pop <- unique(herds$Herd)[c(1,3,6,7)]

# How many years for each
herds %>% 
  group_by(Herd) %>% 
  filter(Status == "Post") %>%
  summarise(num_post = n())

# Loop through to do bootstrap test
herd_test <- do.call(rbind, lapply(1:length(pop), function(i){
  # filter
  dat <- herds %>% filter(Herd == pop[i])
  
  # Test and return
  boot2samp_t_test(dat,
                   "Lambs",
                   "Status") %>%
    as.data.frame() %>%
    mutate(Herd = pop[i])
}))

#------------------------------------------#
# Summarize the data frame for plotting ####
#------------------------------------------#

# Ratios
ratios <- lamb %>%
  group_by(Herd, Status) %>%
  summarise(Ratio = mean(Lambs))

ratios_lat <- lamb %>%
  group_by(Herd, StatusLat) %>%
  summarise(Ratio = mean(Lambs))

ratios3 <- lamb %>% 
  group_by(Herd, Status3) %>%
  summarise(Ratio = mean(Lambs))

#------------------#
# Save the data ####
#------------------#

lamb <- lamb %>%
  select(-c(Positives,
            Samp_Size))

save(lamb,
     ratios,
     ratios_lat,
     full_ci,
     full_ci3,
     spec_ci,
     spec_ci3,
     file = "./Data/Clean/Lamb_Ewe_Plotting_Data.RData")

