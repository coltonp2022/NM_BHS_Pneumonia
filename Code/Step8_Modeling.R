# Modeling lamb:ewe ratios
# Colton Padilla
# 7/31/2023

# Last Updated: 9/25/2023

# Added in range-wide lion removal data for desert models.

#-----------------#
# Housekeeping ####
#-----------------#

# Packages
library(tidyverse)
library(glmmTMB)
library(DHARMa)
library(AICcmodavg)
library(ggpubr)
library(sjPlot)
library(ggeffects)
library(flextable)
#library(parasite)

# Remove other data
rm(list = ls())

# Load the data
lamb <- read_csv("./Data/Clean/BHS_Modeling_Data.csv") %>%
  select(-Timeline)

# Add missing lamb:ewe ratios
miss <- lamb %>% filter(is.na(Lambs)) %>%
  arrange(Herd, Year)
miss$Lambs <- c(61, 33, 57, 71, 71, 27, NA, NA, 18, 56, 27, NA, NA, NA, 13,
                64, 44, 40, 19, 22, 14, 35, 23, 50, 37, 14, 49, NA, 57, NA,
                NA, 75, NA, 82, NA, 14)

# Fix miss
miss <- miss %>% drop_na(Lambs)

# Now add them back
lamb <- lamb %>% filter(!is.na(Lambs)) %>%
  rbind(., miss) %>%
  arrange(Herd, Year) %>%
  mutate(Status = if_else(is.na(Status), "Pre", Status))

lamb$Status <- factor(lamb$Status, levels = c("Pre", 
                                              "Post"))

# Retest the differences
herds <- lamb %>% 
  group_by(Herd) %>%
  filter(!all(Status == "Pre"),
         !all(Status == "Post")) %>%
  ungroup()

# check number of years
herds %>% group_by(Herd, Status) %>% summarise(n = n())

pop <- unique(herds$Herd)[c(1, 3, 4, 6, 7)]

## Loop through to do bootstrap test
#herd_test <- do.call(rbind, lapply(1:length(pop), function(i){
#  # filter
#  dat <- herds %>% filter(Herd == pop[i])
#  
#  # Test and return
#  boot2samp_t_test(dat,
#                   "Lambs",
#                   "Status") %>%
#    as.data.frame() %>%
#    mutate(Herd = pop[i])
#}))
#
## Now test the latir
#LAT <- lamb %>% 
#  filter(Herd == "LAT") %>%
#  mutate(Status = if_else(
#    Year < 2016, "Pre", "Post"
#  ))
#boot2samp_t_test(LAT, "Lambs", "Status")

# Now remove nas
lamb <- lamb %>%
  drop_na(Lambs, Est)

#-------------------------------------#
# Separate the rockies and deserts ####
#-------------------------------------#

# Rocky
rocky <- lamb %>% filter(Subspec == "Rocky") %>%
  select(-c(Month, MonthNum, Area, Subspec, YearMonth)) 
rocky <- rocky[-c(57, 113, 136),] # Remove non matching

# Desert 
desert <- lamb %>% filter(Subspec == "Desert") %>%
  select(-c(Month, MonthNum, Area, Subspec, YearMonth))

#==============================================================================#
# Rocky Mountain Analysis 
#==============================================================================#

#--------------------#
# Scale Variables ####
#--------------------#

# Scale
rocky <- rocky %>%
  mutate(across(c(2, 6:13), scale)) %>%
  mutate(Temp_Class = if_else(PostPart_Temp < -1, "Cold", "Normal/Warm"),
         Temp_Class2 = if_else(PostPart_Temp < 0, "Below Average", "At/Above Average"),
         Lamb_Temp_Class = if_else(PrevLamb_Temp < -1, "Cold", "Normal/Warm"),
         Lamb_Temp_Class2 = if_else(PrevLamb_Temp > 0, "Below Average", "At/Above Average"),
         PostPart_Temp_Class = case_when(
           PostPart_Temp < -1 ~ "Cold",
           PostPart_Temp > 1 ~ "Warm",
           .default = "Normal"
         ),
         PrePart_Temp_Class = case_when(
           PrePart_Temp < -1 ~ "Cold",
           PrePart_Temp > 1 ~ "Warm",
           .default = "Normal"
         ),
         Recovery_Temp_Class = case_when(
           Recovery_Temp < -1 ~ "Cold",
           Recovery_Temp > 1 ~ "Warm",
           .default = "Normal"
         ),
         PrevLamb_Temp_Class = case_when(
           PrevLamb_Temp < -1 ~ "Cold",
           PrevLamb_Temp > 1 ~ "Warm",
           .default = "Normal"
         ))

# Set temperature class factor level
rocky$Lamb_Temp_Class <- factor(rocky$Lamb_Temp_Class,
                                    levels = c("Normal/Warm", "Cold"))
rocky$PrevLamb_Temp_Class <- factor(rocky$PrevLamb_Temp_Class,
                                    levels = c("Normal", "Cold", "Warm"))

#-----------------------#
# View Relationships ####
#-----------------------#

# Get variables
vars <- colnames(rocky)[11:18]

# Loop to view all
for(i in 1:length(vars)){
  print(ggplot(rocky,
           aes_string(vars[i],
               "Lambs")) +
      geom_point() +
      geom_smooth(method = "loess"))
}

# Loop to view interactions
for(i in 1:length(vars)){
  print(ggplot(rocky,
         aes_string(vars[i],
                    "Lambs")) +
    geom_point() +
    geom_smooth(method = "lm"))
}

# Loop to view interactions
for(i in 1:length(vars)){
  print(ggplot(rocky,
               aes_string(vars[i],
                          "Lambs",
                          color = "Status")) +
          geom_point() +
          geom_smooth(method = "lm"))
}

#-----------------------------------------#
# Precipitation and temperature models ####
#-----------------------------------------#

# Model formulas
formulas <- c("Lambs ~ 1",
              "Lambs ~ PrePart_Precip",
              "Lambs ~ PostPart_Precip",
              "Lambs ~ Recovery_Precip",
              "Lambs ~ PrePart_Temp",
              "Lambs ~ PostPart_Temp",
              "Lambs ~ Recovery_Temp",
              "Lambs ~ PrevLamb_Precip",
              "Lambs ~ PrevLamb_Temp",
              "Lambs ~ PrevLamb_Precip + I(PrevLamb_Precip ^ 2)",
              "Lambs ~ PrevLamb_Precip * PrevLamb_Temp_Class",
              "Lambs ~ PrevLamb_Precip * Lamb_Temp_Class",
              "Lambs ~ PrevLamb_Precip * Lamb_Temp_Class2",
              "Lambs ~ PrePart_Precip * PrePart_Temp_Class",
              "Lambs ~ PostPart_Precip * PostPart_Temp_Class",
              "Lambs ~ PostPart_Precip * Temp_Class",
              "Lambs ~ PostPart_Precip * Temp_Class2",
              "Lambs ~ Recovery_Precip * Recovery_Temp_Class",
              "Lambs ~ PostPart_Precip * PostPart_Temp_Class",
              "Lambs ~ PrePart_Precip + I(PrePart_Precip ^ 2)",
              "Lambs ~ PostPart_Precip + I(PostPart_Precip ^ 2)",
              "Lambs ~ Recovery_Precip + I(Recovery_Precip ^ 2)",
              "Lambs ~ PrePart_Precip + I(PrePart_Precip ^ 2) * PrePart_Temp_Class",
              "Lambs ~ PostPart_Precip + I(PostPart_Precip ^ 2) * PostPart_Temp_Class",
              "Lambs ~ Recovery_Precip + I(Recovery_Precip ^ 2) * Recovery_Temp_Class",
              "Lambs ~ PrevLamb_Precip + I(PrevLamb_Precip ^ 2) * PrevLamb_Temp_Class")

# Modnames
names <- gsub("Lambs ~ ", "", formulas)

# Create model list
modlist <- lapply(1:length(formulas), function(i){
  # Create the formula
  x <- as.formula(formulas[i])
  # Create the model
  mod <- glmmTMB(x,
                 data = rocky,
                 family = nbinom2)
  # Return
  return(mod)
})

# Names
names(modlist) <- names

# Check AIC
aictab(modlist)

# Now select informative models
aic <- aictab(modlist)
inform <- aic[which(aic$AICc < aic[aic$Modnames == "1",]$AICc), "Modnames"]
tab <- aic[which(aic$AICc < aic[aic$Modnames == "1",]$AICc),]
tab

# Choose those that improve fit
modlist <- modlist[names(modlist) %in% inform]

# Now summaries
for(i in 1:length(modlist)){
  print(summary(modlist[[i]]))
}

# Plot models to choose interactions
plot_model(modlist[["PrevLamb_Precip + I(PrevLamb_Precip ^ 2)"]],
           type = "pred",
           terms = c("PrevLamb_Precip[all]"))

plot_model(modlist[["PrePart_Precip + I(PrePart_Precip ^ 2)"]],
           type = "pred",
           terms = c("PrePart_Precip[all]"))

# Plot models to choose interactions
plot_model(modlist[["PrevLamb_Precip * Lamb_Temp_Class"]],
           type = "pred",
           terms = c("PrevLamb_Precip", "Lamb_Temp_Class"))

# Plot
plot_model(modlist[["PrevLamb_Precip * Lamb_Temp_Class"]],
           type = "pred",
           terms = c("PrevLamb_Precip", "Lamb_Temp_Class"))

#-------------------------------------#
# Get the best precipitation model ####
#-------------------------------------#

# Model formula
formulas <- c("Lambs ~ PrevLamb_Precip * Lamb_Temp_Class",
              "Lambs ~ PrePart_Precip + I(PrePart_Precip ^ 2)",
              "Lambs ~ PrevLamb_Precip + I(PrevLamb_Precip ^ 2)",
              "Lambs ~ PostPart_Precip",
              "Lambs ~ Recovery_Precip",
              "Lambs ~ PrevLamb_Precip * Lamb_Temp_Class + PrePart_Precip + I(PrePart_Precip ^ 2)",
              "Lambs ~ PrevLamb_Precip * Lamb_Temp_Class + PostPart_Precip",
              "Lambs ~ PrevLamb_Precip * Lamb_Temp_Class + Recovery_Precip",
              "Lambs ~ PrevLamb_Precip * Lamb_Temp_Class + PrePart_Precip",
              "Lambs ~ PrePart_Precip + I(PrePart_Precip ^ 2) + PrevLamb_Precip + I(PrevLamb_Precip ^ 2)",
              "Lambs ~ PrevLamb_Precip + I(PrevLamb_Precip ^ 2) + PrePart_Precip",
              "Lambs ~ PrePart_Precip + I(PrePart_Precip ^ 2) + PostPart_Precip",
              "Lambs ~ PrePart_Precip + I(PrePart_Precip ^ 2) + Recovery_Precip",
              "Lambs ~ PrevLamb_Precip + I(PrevLamb_Precip ^ 2) + PostPart_Precip",
              "Lambs ~ PrevLamb_Precip + I(PrevLamb_Precip ^ 2) + Recovery_Precip",
              "Lambs ~ PrevLamb_Precip * Lamb_Temp_Class + PrePart_Precip + I(PrePart_Precip ^ 2) + PostPart_Precip",
              "Lambs ~ PrePart_Precip + I(PrePart_Precip ^ 2) + PrevLamb_Precip + I(PrevLamb_Precip ^ 2) + PostPart_Precip",
              "Lambs ~ PrevLamb_Precip * Lamb_Temp_Class + PrePart_Precip + I(PrePart_Precip ^ 2) + Recovery_Precip",
              "Lambs ~ PrePart_Precip + I(PrePart_Precip ^ 2) + PrevLamb_Precip + I(PrevLamb_Precip ^ 2) + Recovery_Precip",
              "Lambs ~ PrevLamb_Precip * Lamb_Temp_Class + PrePart_Precip + I(PrePart_Precip ^ 2) + PostPart_Precip + Recovery_Precip",
              "Lambs ~ PrePart_Precip + I(PrePart_Precip ^ 2) + PrevLamb_Precip + I(PrevLamb_Precip ^ 2) + PostPart_Precip + Recovery_Precip",
              "Lambs ~ 1")

# Modnames
names <- gsub("Lambs ~ ", "", formulas)

# Create model list
modlist <- lapply(1:length(formulas), function(i){
  # Create the formula
  x <- as.formula(formulas[i])
  # Create the model
  mod <- glmmTMB(x,
                 data = rocky,
                 family = nbinom2)
  # Return
  return(mod)
})

# Names
names(modlist) <- names

# Check AIC
aictab(modlist)

# Now summaries of top models
summary(modlist[["PrevLamb_Precip * Lamb_Temp_Class + PrePart_Precip + I(PrePart_Precip ^ 2)"]]) # Barely Informative parameter
summary(modlist[["PrevLamb_Precip * Lamb_Temp_Class + PrePart_Precip + I(PrePart_Precip ^ 2) + PostPart_Precip"]]) # Uninformative parameter
summary(modlist[["PrePart_Precip + I(PrePart_Precip ^ 2) + PrevLamb_Precip + I(PrevLamb_Precip ^ 2)"]]) # Barely Informative parameter
summary(modlist[["PrevLamb_Precip * Lamb_Temp_Class + PrePart_Precip"]]) # All very informative. Top model here
summary(modlist[["PrevLamb_Precip * Lamb_Temp_Class + PrePart_Precip + I(PrePart_Precip ^ 2) + Recovery_Precip"]]) # Uninformative parameter
summary(modlist[["PrePart_Precip + I(PrePart_Precip ^ 2) + PrevLamb_Precip + I(PrevLamb_Precip ^ 2) + Recovery_Precip"]]) # Uninformative parameter
summary(modlist[["PrevLamb_Precip * Lamb_Temp_Class + PrePart_Precip + I(PrePart_Precip ^ 2) + PostPart_Precip + Recovery_Precip"]]) # Uninformative parameter
summary(modlist[["PrevLamb_Precip + I(PrevLamb_Precip ^ 2) + PrePart_Precip"]]) # Informative as well
summary(modlist[["PrePart_Precip + I(PrePart_Precip ^ 2) + PrevLamb_Precip + I(PrevLamb_Precip ^ 2) + PostPart_Precip"]]) # Uninformative parameter

mod1 <- plot_model(modlist[["PrevLamb_Precip * Lamb_Temp_Class + PrePart_Precip"]],
                   type = "pred",
                   terms = c("PrevLamb_Precip", "Lamb_Temp_Class")) +
  theme(legend.position = "bottom",
        legend.direction = "horizontal")
mod2 <- plot_model(modlist[["PrevLamb_Precip + I(PrevLamb_Precip ^ 2) + PrePart_Precip"]],
                   type = "pred",
                   terms = "PrevLamb_Precip[all]")
ggarrange(mod1, mod2)

#---------------------------------------------------------------#
# Last Model Set to test for density and status interactions ####
#---------------------------------------------------------------#

# Model formulas
formulas <- c("Lambs ~ Status + PrevLamb_Precip * Lamb_Temp_Class + PrePart_Precip",
              "Lambs ~ Status + Est + PrevLamb_Precip * Lamb_Temp_Class + PrePart_Precip",
              "Lambs ~ Status * Est + PrevLamb_Precip * Lamb_Temp_Class + PrePart_Precip",
              "Lambs ~ Est + PrevLamb_Precip * Lamb_Temp_Class + PrePart_Precip",
              "Lambs ~ PrevLamb_Precip * Lamb_Temp_Class + PrePart_Precip",
              "Lambs ~ Status",
              "Lambs ~ Status + Est",
              "Lambs ~ Status * Est",
              "Lambs ~ Est",
              "Lambs ~ 1")

# Modnames
names <- gsub("Lambs ~ ", "", formulas)

# Create model list
modlist <- lapply(1:length(formulas), function(i){
  # Create the formula
  x <- as.formula(formulas[i])
  # Create the model
  mod <- glmmTMB(x,
                 data = rocky,
                 family = nbinom2)
  # Return
  return(mod)
})

# Names
names(modlist) <- names

# Check AIC
aictab(modlist)

# Recalculate AIC
aic <- data.frame(aictab(cand.set = modlist,
                          modnames = names(modlist)))  %>%
  mutate_if(is.numeric, round, 3) %>%
  rename(Model = Modnames,
         `Log Likelihood` = LL, 
         Weight = AICcWt) %>%
  select(Model, K, `Log Likelihood`, AICc, Delta_AICc, Weight)
aic

aic$Model <- gsub("Est", "Abundance", aic$Model) %>%
  gsub("PrevLamb_Precip", "PLS Precipitation", .) %>%
  gsub("Lamb_Temp_Class", "PLS Temperature", .) %>%
  gsub("PrePart_Precip", "PreP Precipitation", .) %>%
  gsub("Status", "Disease Status", .)

# Remove unneeded columns
aic <- aic %>% 
  select(-c(`Log Likelihood`, AICc))

# Create the table and save it
tab <- flextable(aic) %>%
  set_header_labels(values = list(Delta_AICc = "\u394 AICc"))
  
save_as_docx(tab,
             path = "./TablesFigures/Rocky_AIC_table.docx")

# Summary of top model
summary(modlist[["Status * Est + PrevLamb_Precip * Lamb_Temp_Class + PrePart_Precip"]])
plot_model(modlist[["Status * Est + PrevLamb_Precip * Lamb_Temp_Class + PrePart_Precip"]],
           type = "pred",
           terms = c("PrevLamb_Precip",
                     "Lamb_Temp_Class"))
plot_model(modlist[["Status * Est + PrevLamb_Precip * Lamb_Temp_Class + PrePart_Precip"]],
           type = "pred",
           terms = "PrePart_Precip [all]")
plot_model(modlist[["Status * Est + PrevLamb_Precip * Lamb_Temp_Class + PrePart_Precip"]],
           type = "pred",
           terms = c("Est", "Status"))



# Check fit and collinearity
testDispersion(modlist[["Status * Est + PrevLamb_Precip * Lamb_Temp_Class + PrePart_Precip"]])
performance::check_collinearity(modlist[["Status * Est + PrevLamb_Precip * Lamb_Temp_Class + PrePart_Precip"]])

#-----------------------------------------------#
# Conduct model averaging for the manuscript ####
#-----------------------------------------------#

# Conditionally average the model list
avg_rocky <- MuMIn::model.avg(modlist)

# Create an averaged summary table
summ_rocky <- summary(avg_rocky)$coefmat.subset %>%
  as.data.frame() %>%
  rownames_to_column(var = "Parameter") %>%
  mutate(Parameter = str_remove(Parameter, "cond\\(") %>% str_remove(., "\\)")) %>%
  mutate(across(2:ncol(.), ~round(., 3))) %>%
  mutate(Parameter = case_when(
    Parameter == "(Int)" ~ "Intercept",
    Parameter == "StatusPost" ~ "Pathogen Status",
    Parameter == "Est" ~ "Abundance",
    Parameter == "PrevLamb_Precip" ~ "PLS Precipitation",
    Parameter == "Lamb_Temp_ClassCold" ~ "PLS Temperature",
    Parameter == "PrePart_Precip" ~ "PreP Precipitation",
    Parameter == "Est:StatusPost" ~ "Abundance * Pathogen Status",
    Parameter == "Lamb_Temp_ClassCold:PrevLamb_Precip" ~ "PLS Temperature * PLS Precipitation"
  )) %>%
  rename(Coefficient = Estimate,
         SE = `Std. Error`,
         ASE = `Adjusted SE`,
         Z = `z value`,
         `P-value` = `Pr(>|z|)`)

# Save the model averaged coefficients
save_as_docx(summ_rocky,
             path = "./TablesFigures/Rocky_Model_Summary.docx")

# top Model
mod_rocky <- modlist[["Status * Est + PrevLamb_Precip * Lamb_Temp_Class + PrePart_Precip"]]

# Add random effect
mod2 <- glmmTMB(
  Lambs ~ Status * Est + PrevLamb_Precip * Lamb_Temp_Class + PrePart_Precip + 
    (0 + PrevLamb_Precip * Lamb_Temp_Class|Herd),
  family = nbinom2,
  data = rocky
)

options(scipen = 999)
ranef(mod2)
print(ranef(mod2),simplify=FALSE)

aictab(list(mod_rocky, mod2))
plot_model(mod2, type = "pred", terms = "Est")

#==============================================================================#
# Desert Analysis 
#==============================================================================#

#------------------------------#
# Create the lion covariate ####
#------------------------------#

# Herds
herds <- c("PE", "HAT", "BH", "LH", "LA")

# if else
desert <- desert %>%
  mutate(Lion_REM = if_else(Herd %in% herds & Year > 2000, 1, 
                            if_else(Herd == "CAB" & Year > 2000, 1,
                                    if_else(Herd == "SA" & between(Year, 2001, 2013), 1,
                                            if_else(Herd == "FC" & between(Year, 2005, 2014), 1, 0)))))

# Check
desert %>% 
  select(Herd, Year, Lion_REM) %>%
  arrange(Herd, Year) %>%
  print(n = 124)

table(desert$Lion_REM)

#--------------------#
# Scale Variables ####
#--------------------#

# Scale
desert <- desert %>%
  mutate(across(c(2, 6:13), scale)) %>%
  mutate(Temp_Class = if_else(PostPart_Temp < -1, "Cold", "Normal/Warm"),
         Temp_Class2 = if_else(PostPart_Temp < 0, "Below Average", "At/Above Average"),
         Lamb_Temp_Class = if_else(PrevLamb_Temp < -1, "Cold", "Normal/Warm"),
         Lamb_Temp_Class2 = if_else(PrevLamb_Temp < 0, "Below Average", "At/Above Average"),
         PostPart_Temp_Class = case_when(
           PostPart_Temp < -1 ~ "Cold",
           PostPart_Temp > 1 ~ "Warm",
           .default = "Normal"
         ),
         PrePart_Temp_Class = case_when(
           PrePart_Temp < -1 ~ "Cold",
           PrePart_Temp > 1 ~ "Warm",
           .default = "Normal"
         ),
         Recovery_Temp_Class = case_when(
           Recovery_Temp < -1 ~ "Cold",
           Recovery_Temp > 1 ~ "Warm",
           .default = "Normal"
         ),
         PrevLamb_Temp_Class = case_when(
           PrevLamb_Temp < -1 ~ "Cold",
           PrevLamb_Temp > 1 ~ "Warm",
           .default = "Normal"
         ))


#-----------------------#
# View Relationships ####
#-----------------------#

# Get variables
vars <- colnames(desert)[11:18]

# Loop to view all
for(i in 1:length(vars)){
  print(ggplot(desert,
               aes_string(vars[i],
                          "Lambs")) +
          geom_point() +
          geom_smooth(method = "loess"))
}

# Loop to view interactions
for(i in 1:length(vars)){
  print(ggplot(desert,
               aes_string(vars[i],
                          "Lambs")) +
          geom_point() +
          geom_smooth(method = "lm"))
}

# Loop to view interactions
for(i in 1:length(vars)){
  print(ggplot(desert,
               aes_string(vars[i],
                          "Lambs",
                          color = "Status")) +
          geom_point() +
          geom_smooth(method = "lm"))
}

#----------------------#
# Initial Models ####
#----------------------#

# Toss the Pelocillos 2001 to test for outlier
desert <- desert %>%
  filter(!(Herd == "PE" & Year == 2001))

# Model formulas
# Model formulas
formulas <- c("Lambs ~ 1",
              "Lambs ~ PrePart_Precip",
              "Lambs ~ PostPart_Precip",
              "Lambs ~ Recovery_Precip",
              "Lambs ~ PrePart_Temp",
              "Lambs ~ PostPart_Temp",
              "Lambs ~ Recovery_Temp",
              "Lambs ~ PrevLamb_Precip",
              "Lambs ~ PrevLamb_Temp",
              "Lambs ~ PrevLamb_Precip + I(PrevLamb_Precip ^ 2)",
              "Lambs ~ PrevLamb_Precip * PrevLamb_Temp_Class",
              "Lambs ~ PrevLamb_Precip * Lamb_Temp_Class",
              "Lambs ~ PrevLamb_Precip * Lamb_Temp_Class2",
              "Lambs ~ PrePart_Precip * PrePart_Temp_Class",
              "Lambs ~ PostPart_Precip * PostPart_Temp_Class",
              "Lambs ~ PostPart_Precip * Temp_Class",
              "Lambs ~ PostPart_Precip * Temp_Class2",
              "Lambs ~ Recovery_Precip * Recovery_Temp_Class",
              "Lambs ~ PostPart_Precip * PostPart_Temp_Class",
              "Lambs ~ PrePart_Precip + I(PrePart_Precip ^ 2)",
              "Lambs ~ PostPart_Precip + I(PostPart_Precip ^ 2)",
              "Lambs ~ Recovery_Precip + I(Recovery_Precip ^ 2)",
              "Lambs ~ PrePart_Precip + I(PrePart_Precip ^ 2) * PrePart_Temp_Class",
              "Lambs ~ PostPart_Precip + I(PostPart_Precip ^ 2) * PostPart_Temp_Class",
              "Lambs ~ Recovery_Precip + I(Recovery_Precip ^ 2) * Recovery_Temp_Class",
              "Lambs ~ PrePart_Temp + PrePart_Precip + I(PrePart_Precip ^ 2)")


# Modnames
names <- gsub("Lambs ~ ", "", formulas)

# Create model list
modlist <- lapply(1:length(formulas), function(i){
  # Create the formula
  x <- as.formula(formulas[i])
  # Create the model
  mod <- glmmTMB(x,
                 data = desert,
                 family = nbinom2)
  # Return
  return(mod)
})

# Names
names(modlist) <- names

# Check AIC
aictab(modlist)

summary(modlist[["PrePart_Temp"]])
summary(modlist[["PrePart_Precip + I(PrePart_Precip ^ 2)"]])
summary(modlist[["PrePart_Temp + PrePart_Precip + I(PrePart_Precip ^ 2)"]]) # Weakly informative prepart temp

#--------------------------------------------#
# Finish Environmental covariate modeling ####
#--------------------------------------------#

# Model formulas
formulas <- c("Lambs ~ 1",
              "Lambs ~ PrePart_Temp",
              "Lambs ~ PrePart_Precip + I(PrePart_Precip ^ 2)",
              "Lambs ~ PrePart_Temp + PrePart_Precip + I(PrePart_Precip ^ 2)",
              "Lambs ~ Status + PrePart_Precip + I(PrePart_Precip ^ 2)",
              "Lambs ~ Status + PrePart_Temp + PrePart_Precip + I(PrePart_Precip ^ 2)", # Uninformative param
              "Lambs ~ Status + PrePart_Temp", # Uninformative param
              "Lambs ~ Est + PrePart_Temp",
              "Lambs ~ Est + PrePart_Precip + I(PrePart_Precip ^ 2)",
              "Lambs ~ Est + PrePart_Temp + PrePart_Precip + I(PrePart_Precip ^ 2)",
              "Lambs ~ Lion_REM + PrePart_Temp",
              "Lambs ~ Lion_REM + PrePart_Precip + I(PrePart_Precip ^ 2)",
              "Lambs ~ Lion_REM + PrePart_Temp + PrePart_Precip + I(PrePart_Precip ^ 2)",
              "Lambs ~ Status + Est + PrePart_Temp",
              "Lambs ~ Status + Est + PrePart_Precip + I(PrePart_Precip ^ 2)", # Uninformative Param
              #"Lambs ~ Status + Est + PrePart_Temp + PrePart_Precip + I(PrePart_Precip ^ 2)",
              "Lambs ~ Status + Lion_REM + PrePart_Temp",
              "Lambs ~ Status + Lion_REM + PrePart_Precip + I(PrePart_Precip ^ 2)", # Uninformative Param
              "Lambs ~ Status + Lion_REM + PrePart_Temp + PrePart_Precip + I(PrePart_Precip ^ 2)",
              "Lambs ~ Est + Lion_REM + PrePart_Temp",
              "Lambs ~ Est + Lion_REM + PrePart_Precip + I(PrePart_Precip ^ 2)", 
              "Lambs ~ Est + Lion_REM + PrePart_Temp + PrePart_Precip + I(PrePart_Precip ^ 2)",
              "Lambs ~ Status * Est + PrePart_Temp",
              "Lambs ~ Status * Est + PrePart_Precip + I(PrePart_Precip ^ 2)",
              "Lambs ~ Status * Est + PrePart_Temp + PrePart_Precip + I(PrePart_Precip ^ 2)",
              "Lambs ~ Status + Est + Lion_REM + PrePart_Temp",
              "Lambs ~ Status + Est + Lion_REM + PrePart_Precip + I(PrePart_Precip ^ 2)", 
              "Lambs ~ Status + Est + Lion_REM + PrePart_Temp + PrePart_Precip + I(PrePart_Precip ^ 2)",
              "Lambs ~ Status * Est + Lion_REM + PrePart_Temp",
              "Lambs ~ Status * Est + Lion_REM + PrePart_Precip + I(PrePart_Precip ^ 2)", 
              "Lambs ~ Status * Est + Lion_REM + PrePart_Temp + PrePart_Precip + I(PrePart_Precip ^ 2)",
              "Lambs ~ Lion_REM",
              "Lambs ~ Status",
              "Lambs ~ Est",
              "Lambs ~ Status * Est")

# Modnames
names <- gsub("Lambs ~ ", "", formulas)

# Create model list
modlist <- lapply(1:length(formulas), function(i){
  # Create the formula
  x <- as.formula(formulas[i])
  # Create the model
  mod <- glmmTMB(x,
                 data = desert,
                 family = nbinom2)
  # Return
  return(mod)
})

# Names
names(modlist) <- names

# Check AIC
aictab(modlist)

# Recalculate AIC
aic <- data.frame(aictab(cand.set = modlist,
                          modnames = names(modlist)))  %>%
  mutate_if(is.numeric, round, 3) %>%
  rename(Model = Modnames,
         `Log Likelihood` = LL, 
         Weight = AICcWt) %>%
  select(Model, K, `Log Likelihood`, AICc, Delta_AICc, Weight)
aic

# Summaries
summary(modlist[["Status + PrePart_Temp + PrePart_Precip + I(PrePart_Precip ^ 2)"]]) # Weakly informative temp
summary(modlist[["Status + PrePart_Precip + I(PrePart_Precip ^ 2)"]]) # Highly informative
summary(modlist[["Status + Est + PrePart_Temp + PrePart_Precip + I(PrePart_Precip ^ 2)"]]) # Uninformative param
summary(modlist[["Status + PrePart_Temp"]]) # Informative

# SJplot
plot_model(modlist[["Status + PrePart_Temp + PrePart_Precip + I(PrePart_Precip ^ 2)"]],
           terms = "PrePart_Precip[all]",
           type = "pred")

# gSUB
aic$Model <- gsub("Status", "Disease Status", aic$Model) %>%
  gsub("PrePart_Temp", "PreP Temperature", .) %>%
  gsub("Est", "Abundance", .) %>%
  gsub("PrePart_Precip", "PreP Precipitation", .) %>%
  gsub("\\^ 2", "2", .) %>%
  gsub("I", "", .) %>%
  gsub("\\(", "", .) %>%
  gsub("\\)", "", .)
aic

# Create the table and save it
tab <- aic %>%
  select(-c(`Log Likelihood`, AICc)) %>%
  flextable(.) %>%
  set_header_labels(., values = list(Delta_AICc = "\u394 AICc"))
save_as_docx(tab,
             path = "./TablesFigures/Desert_AIC_table.docx")

# Top model summary
# coefficient tables
tab <- summary(modlist[["Status + PrePart_Temp + PrePart_Precip + I(PrePart_Precip ^ 2)"]])
tab <- data.frame(tab$coefficients$cond) %>%
  rownames_to_column(var = "Coefficient") %>%
  rename(SE = Std..Error, Z = z.value, `P-value` = Pr...z..) %>%
  mutate_if(is.numeric, round, 3)

# Now rename variables
tab$Coefficient <- gsub("StatusPost", "Disease Status (Post-Exposure)", tab$Coefficient) %>%
    gsub("PrePart_Temp", "PreP Temperature", .) %>%
    gsub("PrePart_Precip", "PreP Precipitation", .) %>%
    gsub("\\:", "*", .) %>%
    gsub("I\\(", "", .) %>%
    gsub("\\(", "", .) %>%
    gsub("\\)", "", .)
tab

# Recode
tab$`P-value` <- ifelse(tab$`P-value` < 0.001, "< 0.001", tab$`P-value`)

#---------------------------------#
# Model average for manuscript ####
#---------------------------------#

# Avg
avg_desert <- MuMIn::model.avg(modlist)

# Create an averaged summary table
summ_desert <- summary(avg_desert)$coefmat.subset %>%
  as.data.frame() %>%
  rownames_to_column(var = "Parameter") %>%
  mutate(Parameter = str_remove(Parameter, "cond\\(") %>% str_remove(., "\\)")) %>%
  mutate(across(2:ncol(.), ~round(., 3))) %>%
  mutate(Parameter = case_when(
    Parameter == "(Int)" ~ "Intercept",
    Parameter == "StatusPost" ~ "Pathogen Status",
    Parameter == "PrePart_Precip" ~ "PreP Precipitation",
    Parameter == "Est:StatusPost" ~ "Abundance * Pathogen Status",
    Parameter == "PrePart_Temp" ~ "PreP Temperature",
    Parameter == "I(PrePart_Precip^2)" ~ "PreP Precipitation2",
    Parameter == "Lion_REM" ~ "Lion Removal",
    Parameter == "Est" ~ "Abundance",
    Parameter == "Est:StatusPost" ~ "Abundance * Pathogen Status"
  )) %>%
  rename(Coefficient = Estimate,
         SE = `Std. Error`,
         ASE = `Adjusted SE`,
         Z = `z value`,
         `P-value` = `Pr(>|z|)`)

# Save
save_as_docx(flextable(summ_desert), 
             path = "./TablesFigures/Desert_Model_Summary.docx")

# Mod des
mod_des <- modlist[["Status + PrePart_Temp + PrePart_Precip + I(PrePart_Precip ^ 2)"]]

#----------------#
# Save models ####
#----------------#

save(mod_des, avg_desert, avg_rocky, mod_rocky, rocky, desert, file = "./Data/Top_Models.RData")

