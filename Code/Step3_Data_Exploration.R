# Bighorn Sheep Mycoplasma Data Exploration
# Colton Padilla
# 4/10/2023

# Last Updated: 4/24/2023

#-----------------#
# Housekeeping ####
#-----------------#

# Packages
library(tidyverse)
library(parasite)
library(flextable)

# Bring in the database
myco <- read_csv("./Data/Clean/BHS_Myco_Data_Clean.csv")

# Recode one sample
myco$Herd <- ifelse(myco$Herd == "RR", "RRiv", myco$Herd)

# Remove columns with all NAs
all_na <- function(x) any(!is.na(x))
myco <- myco %>% select_if(all_na)

# Remove other unneeded columns
myco <- myco %>% 
  select(-c(CaseID, Report_Date, Lab_Notes, Class,
            Rump_mm))

# Categorize
myco$cat <- paste0(myco$Age2, "(", myco$Sex, ")")

# Create a Year column
myco$Year <- ifelse(is.na(myco$Year), as.Date(myco$Received_Date, "%m/%d/%Y") %>%
  format.Date(., "%Y"),
  myco$Year)

# Now get RGG from the notes section
myco$Herd <- ifelse(grepl("RGG", myco$Notes),
                    "RGG",
                    myco$Herd)

# Recode one age
myco$Age2 <- ifelse(myco$Age2 == "4-5",
                    "A",
                    myco$Age2)

# Write this data frame for analysis
write_csv(myco,
          "./Data/Clean/BHS_Myco_Analysis_Data.csv")

# Count active infections from SA
myco %>%
  filter(Source == "SA" | Herd == "SA") %>%
  group_by(Year, PCR_Result) %>% 
  count(PCR_Result)

# View SA inhibition
myco %>% filter(!is.na(Received_Date),
                !is.na(ELISA_Level),
                (Herd == "SA" | Source == "SA")) %>%
  ggplot(.,
         aes(x = as.numeric(Year),
             y = ELISA_Level)) +
  geom_smooth(method = "lm",
              fill = "darkgreen",
              color = "black") +
  geom_point() +
  theme_bw() 

# Average inhibition
myco %>% filter(Source == "SA" | Herd == "SA") %>%
  drop_na(ELISA_Level) %>%
  group_by(Year) %>%
  summarise(avg_inhib = mean(ELISA_Level))

# Average positivity
myco %>% filter((Source == "SA" |Herd == "SA")) %>%
  drop_na(ELISA_Result) %>%
  group_by(Year) %>%
  ggplot(.,
         aes(x = as.numeric(Year),
             y = ELISA_Result)) +
  geom_point(position = position_jitter(width = 0.1)) +
  geom_smooth(method = "lm",
              color = "black",
              fill = "darkgreen")
 

# Remove indeterminate pcr
myco_pcr <- myco %>% filter(!PCR_Result == "0.5") 

# Remove indeterminate elisa
myco_elisa <- myco %>% filter(!ELISA_Result == "0.5") 

# Remove individuals from RR source
myco <- myco %>% filter(is.na(Source) |Source != "RRoc")

# Recode fras
unique(myco$Source)
myco$Source <- gsub("FRA", "FC", myco$Source)

# Seperate data for only infected pops
# Build timeline table
time <- data.frame(
  Herd = c("SA", "FC", "CAB", "LA", "SAC", "RGG", "RRiv", "WH", "CUL"),
  Year = c(2007, 2017, 2019, 2017, 2018, 2020, 2022, 2022, 2022)
)

# Filter to infected pops
infect <- myco %>% filter(Herd %in% time$Herd |
                            Source %in% time$Herd)

# Remove individuals from <2020 in RGG
infect <- infect %>% 
  filter(!(Herd == "RGG" & Year %in% c(2018, 2019))) %>%
  filter(!Herd == "JZ") %>%
  filter(!(Herd == "RRiv" & Year == 2020))


# Set factors
myco_elisa$cat <- factor(myco_elisa$cat,
                         levels = c("L(F)", "L(M)",
                                    "Y(F)", "Y(M)",
                                    "A(F)", "A(M)"))

myco_pcr$cat <- factor(myco_pcr$cat,
                         levels = c("L(F)", "L(M)",
                                    "Y(F)", "Y(M)",
                                    "A(F)", "A(M)"))





#-----------------------------------#
# Check for age and sex patterns ####
#-----------------------------------#

# PCR Prevalence
pcr_sex <- prevCI(myco_pcr, 
                  "PCR_Result", 
                  method = "Blaker", 
                  group = "cat") %>%
  mutate(Test = "PCR")

# ELISA Prevalence
elisa_sex <- prevCI(myco_elisa, 
                    "ELISA_Result", 
                    method = "Blaker", 
                    group = "cat") %>%
  mutate(Test = "Elisa")

# Bind
sex <- rbind(pcr_sex, elisa_sex)

# Plot
ggplot(sex,
       aes(x = Group,
           y = Naive_Prev,
           color = Test)) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = Lower,
                    ymax = Upper),
                width = 0.2,
                position = position_dodge(width = 0.5)) +
  geom_text(aes(y = Upper + 0.02,
                label = N),
            position = position_dodge(width = 0.5)) +
  theme_bw()

#----------------------#
# Sex by Population ####
#----------------------#

# Pop
pop <- unique(myco$Herd)

# Loop through populations to build sex based estimate
sexpop_pcr <- lapply(1:length(pop), function(i){
  # Create a new frame
  df <- myco_pcr %>% filter(Herd == pop[i])
  if(all(is.na(df$Herd))){
    NULL
  } else {
    pcr <- prevCI(df,
                  "PCR_Result",
                  method = "Blaker", 
                  group = "cat") %>%
      mutate(Test = "PCR",
             Herd = pop[i])
    return(pcr)
  }
})

# Remove null list objects and Rbind
sexpop_pcr <- sexpop_pcr[!sapply(sexpop_pcr, is.null)] %>%
  do.call(rbind, .)

# Loop
sexpop_elisa <- lapply(1:length(pop), function(i){
  # Create a new frame
  df <- myco_elisa %>% filter(Herd == pop[i])
  if(all(is.na(df$Herd))){
    NULL
  } else {
    pcr <- prevCI(df,
                  "ELISA_Result",
                  method = "Blaker", 
                  group = "cat") %>%
      mutate(Test = "ELISA",
             Herd = pop[i])
    return(pcr)
  }
})

# Remove null list objects and Rbind
sexpop <- sexpop_elisa[!sapply(sexpop_elisa, is.null)] %>%
  do.call(rbind, .) %>%
  rbind(., sexpop_pcr)

# Plot
ggplot(sexpop,
       aes(x = Group,
           y = Naive_Prev,
           color = Test)) +
  geom_point(position = position_dodge(width = 0.5),
             size = 2) +
  geom_errorbar(aes(ymin = Lower,
                    ymax = Upper),
                width = 0.2,
                position = position_dodge(width = 0.5),
                linewidth = 1) +
  geom_text(aes(y = Upper + 0.2,
                label = N),
            position = position_dodge(width = 0.5)) +
  theme_bw() +
  facet_wrap(~ Herd)

#-------------------------------#
# Check for Spatial Patterns ####
#-------------------------------#

# Herd PCR
pcr_herd <- prevCI(myco_pcr,
                   "PCR_Result",
                   method = "Blaker", 
                   group = "Herd") %>%
  mutate(Test = "PCR")

# Herd Elisa
elisa_herd <- prevCI(myco_elisa,
                     "ELISA_Result",
                     method = "Blaker",
                     group = "Herd") %>%
  mutate(Test = "Elisa")

# Bind
herd <- rbind(pcr_herd, elisa_herd)

# Plot
ggplot(herd %>% filter(!N < 3),
       aes(x = Group,
           y = Naive_Prev,
           color = Test)) +
  geom_point(position = position_dodge(width = 0.5),
             size = 2) +
  geom_errorbar(aes(ymin = Lower,
                    ymax = Upper),
                width = 0.2,
                position = position_dodge(width = 0.5),
                linewidth = 1) +
  geom_text(aes(y = Upper + 0.05,
                 label = N),
             position = position_dodge(width = 0.5)) +
  theme_bw()

#-----------------------------------------#
# Check for Yearly Population Patterns ####
#-----------------------------------------#

# Get a vector of years
years <- sort(unique(c(unique(myco_pcr$Year), unique(myco_elisa$Year))))

# Apply over years
yearly_pcr <- lapply(1:length(years), function(i){
  # Create a new frame
  df <- myco_pcr %>% filter(Year == years[i])
  if(all(is.na(df$Herd))){
    NULL
  } else {
    pcr <- prevCI(df,
                  "PCR_Result",
                  method = "Blaker", 
                  group = "Herd") %>%
      mutate(Test = "PCR",
             Year = years[i])
    return(pcr)
  }
})

# Remove null list objects and Rbind
yearly_pcr <- yearly_pcr[!sapply(yearly_pcr, is.null)] %>%
  do.call(rbind, .)

# Yearly elisa
yearly_elisa <- lapply(1:length(years), function(i){
  # Create a new frame
  df <- myco_elisa %>% filter(Year == years[i])
  if(all(is.na(df$Herd))){
    NULL
  } else {
    pcr <- prevCI(df,
                  "ELISA_Result",
                  method = "Blaker", 
                  group = "Herd") %>%
      mutate(Test = "ELISA",
             Year = years[i])
    return(pcr)
  }
})

# Remove Null and Bind
yearly <- yearly_elisa[!sapply(yearly_elisa, is.null)] %>%
  do.call(rbind, .) %>%
  rbind(., yearly_pcr)

# Plot by year and herd to see what we have 
ggplot(yearly,
       aes(x = Year,
           y = Naive_Prev,
           color = Test)) +
  geom_point(size = 2,
             position = position_dodge(width = 1)) +
  geom_errorbar(aes(ymin = Lower,
                    ymax = Upper),
                width = 0.2,
                linewidth = 1,
                position = position_dodge(width = 1)) +
  geom_hline(yintercept = 1) +
  geom_label(aes(y = Upper + 0.12,
                 label = N),
             position = position_dodge(width = 1)) +
  facet_wrap(~ Group) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))

#---------------------------#
# Check for Age Patterns ####
#---------------------------#

# PCR Age
pcr_age <- prevCI(infect %>% filter(PCR_Result != 0.5), 
                  "PCR_Result",
                  method = "Blaker",
                  group = "Age2") %>%
  mutate(Test = "PCR")

# Test
FishTest(infect %>% filter(PCR_Result != 0.5) %>%
           drop_na(Age2),
         "PCR_Result",
         group = "Age2")

# Elisa Age
elisa_age <- prevCI(infect %>% filter(ELISA_Result != 0.5),
                    "ELISA_Result",
                    method = "Blaker",
                    group = "Age2") %>%
  mutate(Test = "cELISA")

# Test
FishTest(infect %>% filter(ELISA_Result != 0.5) %>%
           drop_na(Age2),
         column = "ELISA_Result",
         group = "Age2")

# Bind 
age <- rbind(pcr_age, elisa_age)

# Recode the groups
age$Group <- ifelse(age$Group == "Y",
                    "Yearling",
                    ifelse(age$Group == "A", "Adult", "Lamb"))

# Reset factor levels
age$Group <- factor(age$Group,
                    levels = c("Lamb", "Yearling", "Adult"))

# Recode the test
age$Test <- ifelse(age$Test == "cELISA",
                   "Seroprevalence",
                   "Prevalence")

# Plot
ggplot(age %>% drop_na(Group), 
       aes(x = Group,
           y = Naive_Prev)) +
  geom_point(position = position_dodge(width = 0.5),
             size = 2) +
  geom_errorbar(aes(ymin = Lower,
                    ymax = Upper),
                width = 0.2,
                position = position_dodge(width = 0.5),
                linewidth = 1) +
  theme_bw() + 
  facet_wrap(~ Test) +
  labs(y = "Estimated Postitivity",
       x = "Age Class") +
  theme(axis.title.y = element_text(margin = margin(r = 10),
                                    size = 12),
        axis.title.x = element_text(margin = margin(t = 10),
                                    size = 12),
        strip.text = element_text(size = 12))

# Save
ggsave("./TablesFigures/Age_Specific_Disease.jpg",
       height = 7,
       width = 7,
       dpi = 600)


# Reorder
age <- age[order(age$Test, age$Group), ]

# Edit columns and Group the data
age <- age %>% 
  drop_na(Group) %>%
  rename(Prevalence = Naive_Prev,
         `Sample Size` = N) %>%
  as_grouped_data(., "Test")

# Create a table
tab <- flextable(age)

# Write
save_as_docx(tab,
             path = "./TablesFigures/Age_Prevalence.docx")
  
r#------------------------------#
# Check for Weight Patterns ####
#------------------------------#

# PCR Density Plots
ggplot(myco %>% filter(!PCR_Result == "0.5"),
       aes(x = Wt_kg,
           fill = as.character(PCR_Result))) +
  geom_density()

# ELISA Density Plots
ggplot(myco %>% filter(!ELISA_Result == "0.5"),
       aes(x = Wt_kg,
           fill = as.character(ELISA_Result))) +
  geom_density(alpha = 0.5) +
  geom_density(data = myco %>% filter(!ELISA_Result == "0.5"),
               mapping = aes(x = Wt_kg),
               fill = "grey",
               color = "red",
               alpha = 0.5,
               inherit.aes = F)

#---------------#
# Full Cases ####
#---------------#

# If Elisa and PCR
myco_full <- myco[(!is.na(myco$PCR_Result) & !is.na(myco$ELISA_Result)),]
# Remove indeterminate pcr
myco_full <- myco_full %>% filter(!PCR_Result == "0.5")
# Remove indeterminate elisa
myco_full <- myco_full %>% filter(!ELISA_Result == "0.5")

# PrevCI
prev_pcr <- prevCI(myco_full, "PCR_Result", method = "Blaker") %>%
  mutate(Test = "PCR")
prev_elisa <- prevCI(myco_full, "ELISA_Result", method = "Blaker") %>%
  mutate(Test = "Elisa")

# Bind 
prev <- rbind(prev_pcr, prev_elisa)

# Plot
ggplot(prev,
       aes(x = Test,
           y = Naive_Prev)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = Lower,
                    ymax = Upper),
                width = 0.2,
                linewidth = 1) +
  theme_bw()


