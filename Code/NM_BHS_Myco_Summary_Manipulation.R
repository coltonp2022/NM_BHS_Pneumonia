# NM BHS Mycoplasma Summary Manipulation
# Colton Padilla
# 2/22/2023

#-----------------#
# Housekeeping ####
#-----------------#

# Load packages
library(tidyverse) # Data Manipulation
library(parasite) # Prevalence CI's

# Set the theme for plots
theme_set(theme_bw())
theme_replace(panel.grid.minor = element_blank(),
              panel.grid.major = element_blank(),
              axis.text.x = element_text(angle = 90,
                                         size = 10),
              axis.title = element_text(size = 14,
                                        face = "bold"),
              axis.title.y = element_text(margin = margin(r = 10),
                                          angle = 90),
              axis.title.x = element_text(margin = margin(t = 10)),
              plot.caption.position = "plot",
              plot.caption = element_text(size = 8,
                                          hjust = 1))

# Load the data
bhs <- read_csv("./Data/Clean/NM_BHS_Myco_Data.csv")

#-----------------------#
# Begin manipulation ####
#-----------------------#

# Remove rows
myco <- bhs %>%
  filter(!species == "DS")

# Pull out anything with no total animals
missing <- myco[is.na(myco$`Total Animals`),] # Missing data
myco <- myco[!is.na(myco$`Total Animals`),] # Remove missing data from orig
random <- myco[myco$Test == "M. Ovi, Bluetongue, EHD", ] # Randoms
myco <- myco[!myco$Test == "M. Ovi, Bluetongue, EHD", ] # Remove random from orig

# Now replicate the rows for a similar data set
myco <- myco[rep(seq(nrow(myco)), myco$`Total Animals`),]

# Add test notes column for certain things
myco$test_notes <- ifelse(myco$Type == "Nasal Swab/PCR", "Test type was 'Nasal Swab/PCR'", NA)

# Now recode test values
unique(myco$Type)
myco$Type <- ifelse(myco$Type == "P", "PCR",
                    ifelse(myco$Type == "PCr", "PCR",
                           ifelse(myco$Type == "Nasal Swab/PCR", "PCR", myco$Type)))

# Recode the areas
unique(myco$Range)
myco$Range <- ifelse(myco$Range == "Wheeler", "WHE (RM)",
                     ifelse(myco$Range == "Pecos", "PEC (RM)",
                            ifelse(myco$Range == "Latir", "LAT (RM)",
                                   ifelse(myco$Range == "Red River", "RRiv (RM)",
                                          ifelse(myco$Range == "San Andres", "SAA (DS)", myco$Range)))))
unique(myco$Range)
myco$Range <- ifelse(myco$Range == "Hatchets", "HAT (DS)",
                     ifelse(myco$Range == "Ladrones", "LAD (DS)",
                            ifelse(myco$Range == "Caballos (To Sacs)", "CAB (DS)",
                                   ifelse(myco$Range == "Rio Grande Gorge", "RGG (RM)",
                                          ifelse(myco$Range == "manzano", "MAN (RM)", myco$Range)))))
unique(myco$Range)
myco$Range <- ifelse(myco$Range == "Caballos", "CAB (DS)",
                     ifelse(myco$Range == "San Andres (to  Sacs)", "SAA (DS)",
                            ifelse(myco$Range == "Red Rock (to Pels)", "RRoc (DS)",
                                   ifelse(myco$Range == "Fra Cristobal", "FC (DS)",
                                          ifelse(myco$Range == "Dry Cimarron", "DRY (RM)", myco$Range)))))
unique(myco$Range)
myco$Range <- ifelse(myco$Range == "RR/PE", "RRoc (DS)",
                     ifelse(myco$Range == "San Francisco", "SF (RM)", 
                            ifelse(myco$Range == "RGG", "RGG (RM)", myco$Range)))

# Now get Elisa and PCR Data sets
elisa <- myco %>%
  filter(Type == "ELISA")

pcr <- myco %>%
  filter(Type == "PCR")

# Now for each of these go through and add positivity results
elisa <- elisa %>%
  group_by(Year, Range, `Positive Results`, `Collection occurrence`) %>%
  mutate(Result = if_else(`Positive Animals` > 0, 
                          c(rep(1, unique(`Positive Animals`)),
                            rep(0, unique(`Total Animals`) - unique(`Positive Animals`))),
                          rep(0, unique(`Total Animals`)))) %>%
  relocate(Year, 
           Range, 
           `Positive Results`, 
           Result) %>%
  select(-c(`...12`, `Positive Animals`, `Total Animals`, species))

pcr <- pcr %>%
  group_by(Year, Range, `Positive Results`, `Collection occurrence`) %>%
  mutate(Result = if_else(`Positive Animals` > 0, c(rep(1, unique(`Positive Animals`)),
                                                    rep(0, unique(`Total Animals`) - unique(`Positive Animals`))),
                          rep(0, unique(`Total Animals`)))) %>%
  relocate(Year, Range, `Positive Results`, Result) %>%
  select(-c(`...12`, `Positive Animals`, `Total Animals`, species))

# Now lets make some summary statistics for prevalence by range
range_elisa <- prevCI(elisa, "Result", group = "Range", method = "Blaker")
range_pcr <- prevCI(pcr, "Result", group = "Range", method = "Blaker")

# Now by year
year_elisa <- prevCI(elisa, "Result", group = "Year", method = "Blaker")
year_pcr <- prevCI(pcr, "Result", group = "Year", method = "Blaker")

# Now add test type and bind everything together
elisa_list <- lapply(list(range_elisa,
            year_elisa),
       mutate,
       Type = "ELISA")
pcr_list <- lapply(list(range_pcr,
                   year_pcr),
              mutate,
              Type = "PCR")

range_list <- rbind(elisa_list[[1]], pcr_list[[1]])
year_list <- rbind(elisa_list[[2]], pcr_list[[2]])

# Range Point CI Plot
ggplot(range_list,
       aes(y = Naive_Prev,
           x = Group,
           color = Type)) +
  geom_point(position = position_dodge2(width = 0.9),
           size = 2) +
  geom_errorbar(aes(ymin = Lower,
                    ymax = Upper),
                position = position_dodge2(width = 1),
                size = 1) +
  geom_label(aes(label = N,
                 y = Upper),
             color = "black",
             vjust = -0.2, 
             position = position_dodge2(width = 1)) +
  labs(x = "Range",
       y = "Positivity Estimate",
       color = "Test Type",
       caption = "Error bars represent 95% confidence intervals for \n 
       positivity of the range estimated using Blaker's method") +
  ylim(c(0, max(year_list$Upper) + 0.05))

# Year Point CI plot
ggplot(year_list,
       aes(y = Naive_Prev,
           x = Group,
           color = Type)) +
  geom_point(position = position_dodge2(width = 0.9),
             size = 2) +
  geom_errorbar(aes(ymin = Lower,
                    ymax = Upper),
                position = position_dodge2(width = 1),
                size = 1) +
  geom_label(aes(label = N,
                 y = Upper),
             color = "black",
             vjust = -0.2, 
             position = position_dodge2(width = 1)) +
  labs(x = "Year",
       y = "Positivity Estimate",
       color = "Test Type",
       caption = "Error bars represent 95% confidence intervals for \n 
       positivity of the year estimated using Blaker's method") +
  ylim(c(0, max(year_list$Upper) + 0.05))

# Range Bar plot
ggplot(range_list,
       aes(y = Naive_Prev,
           x = Group,
           fill = Type)) +
  geom_bar(position = "dodge",
           stat = "identity") +
  geom_label(aes(label = N,
                 y = Naive_Prev),
             vjust = -0.2, 
             position = position_dodge(width = 1)) +
  labs(x = "Range",
       y = "Positivity Estimate",
       color = "Test Type") +
  ylim(c(0, max(year_list$Upper) + 0.05))

# Year Bar plot
ggplot(year_list,
       aes(y = Naive_Prev,
           x = Group,
           fill = Type)) +
  geom_bar(position = "dodge",
           stat = "identity")+
  geom_label(aes(label = N,
                 y = Naive_Prev),
             vjust = -0.2, 
             position = position_dodge(width = 1)) +
  labs(x = "Year",
       y = "Positivity Estimate",
       color = "Test Type") +
  ylim(c(0, max(year_list$Naive_Prev) + 0.05))

