# Add Sheep Data to the Myco Data
# Colton Padilla
# 4/10/2023

# Last Updated: 4/10/2023

#-----------------#
# Housekeeping ####
#-----------------#

# Packages
library(tidyverse)

# Load the capture data
sheep <- read_csv("./Data/Capture Database updated 2023.csv")

# Bring in Myco data
myco <- read_csv("./Data/BHS_Myco_Testing_4_7_2023.csv") %>%
  select(-c(`...12`, `...13`, `...14`))

# Remove duplicate individuals until testing sorted
test_issue <- myco %>% filter((AnimalID == "P1014" |
                                  AnimalID == "P762" |
                                  AnimalID == "R578")) %>%
  mutate(Report_Date = as.Date(Report_Date, format = "%m/%d/%Y"),
         Received_Date = as.Date(Received_Date, format = "%m/%d/%Y"))
myco <- myco %>% filter(!(AnimalID == "P1014" |
                        AnimalID == "P762" |
                          AnimalID == "R578"))

#-------------------------#
# Clean the sheep data ####
#-------------------------#

# Select only the needed columns
sheep <- sheep %>%
  select(`Capture Date`, `Permanent Id`, `Ear Tag`, Source, Herd, `Weight (kg)`,
         Sex, Age, Class, `Rump (mm)`, `Preg?`)

# Rename columns
sheep <- sheep %>%
  rename(Cap_Date = `Capture Date`,
         AnimalID = `Permanent Id`,
         EarTag = `Ear Tag`,
         Wt_kg = `Weight (kg)`,
         Rump_mm = `Rump (mm)`,
         Preg = `Preg?`)

# Remove random strings from Eartag
sheep$EarTag <- gsub("_right", "", sheep$EarTag)
sheep$EarTag <- gsub("_RIGHT", "", sheep$EarTag)

# Sort out double ear tags
sheep$EarTag <- sub("\\/.*", "", sheep$EarTag)
sheep$EarTag <- sub("\\,.*", "", sheep$EarTag)

# Now capitalization
sheep$EarTag <- toupper(sheep$EarTag)

# Now Filter to animal IDs in the myco data
sheep2 <- sheep %>%
  filter(EarTag %in% myco$AnimalID)

# Now make dates
sheep2$Cap_Date <- as.Date(sheep2$Cap_Date, format = "%m/%d/%Y")
myco$Received_Date <- as.Date(myco$Received_Date, format = "%m/%d/%Y")
myco$Report_Date <- as.Date(myco$Report_Date, format = "%m/%d/%Y")

# Remove all Y2, same date, same capture
y2 <- sheep2 %>% filter(EarTag == "Y2")

# Now join the data together and then get the animal captured closest to test date
myco2 <- left_join(myco, 
                   sheep2 %>% filter(!EarTag == "Y2"),
                   by = c("AnimalID" = "EarTag")) %>%
  mutate(date_diff = abs(Cap_Date - Received_Date)) %>%
  group_by(AnimalID) %>%
  filter(date_diff == min(date_diff)) %>%
  select(-c(date_diff, VetPak, PermID)) %>%
  rename(PermID = AnimalID.y)

# Now join all y2 data together
y2_myco <- left_join(myco %>% filter(AnimalID == "Y2"),
                     y2,
                     by = c("AnimalID" = "EarTag", 
                            "PermID" = "AnimalID"),
                     keep = T) %>%
  select(-c(VetPak, PermID, EarTag)) %>%
  rename(AnimalID = AnimalID.x,
         PermID = AnimalID.y)

# Now bind this back together
myco2 <- rbind(myco2, y2_myco)

# Remove these animals from the original myco data
myco <- myco %>% filter(!AnimalID %in% myco2$AnimalID)

# Bind these back to the data
myco <- plyr::rbind.fill(myco, myco2)

# Bind the removed animals back in
myco <- plyr::rbind.fill(myco, test_issue)

# Make the date a character
myco <- myco %>%
  mutate(across(c(Received_Date, Report_Date, Cap_Date), as.character))

# Create a new folder for final data
if(!dir.exists("./Data/Clean")){
  dir.create("./Data/Clean")
}

# Write the data
write_csv(myco,
          "./Data/Clean/BHS_Myco_Data.csv")

# Remove unneeded stuff
rm(myco, sheep2, test_issue, y2, y2_myco)

#-------------------------------------------#
# Bind the results to the sheep database ####
#-------------------------------------------#

# Make dates match up
sheep$Cap_Date <- as.Date(sheep$Cap_Date,
                          format = "%m/%d/%Y") %>%
  as.character()

myco2$Cap_Date <- as.character(myco2$Cap_Date)

# Get all the animals matching myco2
sheep2 <- sheep %>%
  semi_join(., myco2,
            by = c("EarTag" = "AnimalID",
                   "Cap_Date"))

# Remove all animals matching myco2
sheep <- sheep %>%
  anti_join(., myco2,
            by = c("EarTag" = "AnimalID",
                   "Cap_Date"))

# Check 
table(sheep2$EarTag %in% myco2$AnimalID)

# Now join myco data to sheep2
sheep2 <- left_join(sheep2,
                    myco2 %>% select(AnimalID, PCR_Result, ELISA_Result,
                                     Cap_Date, PermID),
                    by = c("EarTag" = "AnimalID",
                           "Cap_Date",
                           "AnimalID" = "PermID")) %>%
  select(Cap_Date, AnimalID, EarTag, PCR_Result, ELISA_Result)

# Bring the original sheep data back 
sheep <- read_csv("./Data/Capture Database updated 2023.csv") 

# Reformat the sheep2 date for the join
sheep2$Cap_Date <- as.Date(sheep2$Cap_Date) %>%
  format.Date(., "%m/%d/%Y")

# Join the sheep2 data 
sheep <- left_join(sheep,
                   sheep2, 
                   by = c("Ear Tag" = "EarTag",
                          "Capture Date" = "Cap_Date",
                          "Permanent Id" = "AnimalID"))

# Recode 
sheep$`Movi PCR` <- ifelse(!is.na(sheep$PCR_Result), sheep$PCR_Result, sheep$`Movi PCR`)
sheep$`Movi Sero` <- ifelse(!is.na(sheep$ELISA_Result), sheep$ELISA_Result, sheep$`Movi Sero`)

# Write that data
write_csv(sheep,
          "./Data/Capture Database updated 2023_CJP.csv")
