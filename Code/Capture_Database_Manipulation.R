# Capture Database Manipulation
# Colton Padilla
# 4/7/2023

#-----------------#
# Housekeeping ####
#-----------------#

# Load Packages
library(tidyverse)

# Load the capture data
sheep <- read_csv("./Data/Capture Database updated 2023.csv")

# Bring in Myco data
myco <- read_csv("./Data/BHS_Myco_Testing_4_7_2023.csv") %>%
  select(-c(`...12`, `...13`, `...14`))

# Check for duplicated rows 
table(duplicated(sheep))
dups <- sheep[(duplicated(sheep) | duplicated(sheep, fromLast = T)),]

# Remove a bunch of rows with all NAs
sheep <- sheep[!rowSums(is.na(sheep)) >= 43,]

# Now get a more simple database to work with
simple_sheep <- sheep %>%
  select(-c(Subspecies, Source, Herd, Method, 8:43)) %>%
  drop_na(`Permanent Id`, `Ear Tag`)

#-----------------------------------------#
# Place myco results into capture data ####
#-----------------------------------------#

# Remove random string from Eartag
simple_sheep$`Ear Tag` <- gsub("_right", "", simple_sheep$`Ear Tag`)
simple_sheep$`Ear Tag` <- gsub("_RIGHT", "", simple_sheep$`Ear Tag`)

# Check for duplicate ear tags
table(duplicated(simple_sheep$`Ear Tag`))

# Pull out those duplicates to check against the disease data
dup_sheep <- simple_sheep[duplicated(simple_sheep$`Ear Tag`), "Ear Tag"] %>%
  drop_na() %>%
  pull()

# Get Y2 specifically, two animals, same capture, same day
y2_myco <- myco %>% filter(AnimalID == "Y2")

# Remove the 3 duplicate test results until sorted out
myco2 <- myco[!(duplicated(myco$AnimalID) | duplicated(myco$AnimalID, fromLast = T)),]

# Join y2 back
myco2 <- rbind(myco2, y2_myco)
rm(y2_myco)

# Pull duplicate IDs from myco2
myco_dup <- myco2 %>% filter(AnimalID %in% dup_sheep)

# Remove duplicates from original myco data
myco2 <- myco2 %>% filter(!AnimalID %in% dup_sheep)

# Remove the original duplicates list
rm(dup_sheep)

# Pull the duplicates from the sheep data to sort it out
sheep_dup <- simple_sheep[simple_sheep$`Ear Tag` %in% myco_dup$AnimalID,] 

# Remove Y2, 2 individuals, same capture, same day
y2 <- sheep_dup %>% filter(`Ear Tag` == "Y2")
sheep_dup <- sheep_dup %>% filter(!`Ear Tag` == "Y2")

# Make columns dates for sorting out issues
sheep_dup$`Capture Date` <- as.Date(sheep_dup$`Capture Date`, 
                                    format = "%m/%d/%Y") 
myco_dup$Received_Date <- as.Date(myco_dup$Received_Date,
                                  format = "%m/%d/%Y")

# Now join the data together and then get the animal captured closest to test date
sheep_dup <- left_join(sheep_dup, 
                       myco_dup,
                       by = c(`Ear Tag` = "AnimalID")) %>%
  mutate(date_diff = abs(`Capture Date` - Received_Date)) %>%
  group_by(`Ear Tag`) %>%
  filter(date_diff == min(date_diff)) %>%
  select(-date_diff)

# Reformat the date
sheep_dup$`Capture Date` <- format.Date(sheep_dup$`Capture Date`,
                                        "%m/%d/%Y")

# Reformat simple_sheep date to fix compatability
simple_sheep$`Capture Date` <- as.Date(simple_sheep$`Capture Date`, format = "%m/%d/%Y") %>%
  format.Date(., "%m/%d/%Y")

# Now join Y2 
y2 <- left_join(y2,
                myco_dup, 
                by = c(`Permanent Id` = "PermID"))

# Now bind these together
sheep_dup <- rbind(sheep_dup,
                   y2) %>%
  select(-AnimalID)

# Join non duplicate data to the sheep data
simple_sheep <- left_join(simple_sheep, 
                   myco2,
                   by = c(`Ear Tag` = "AnimalID"))
  
# Now remove the rows from sheep and re-add them from the duplicate
simple_sheep <- simple_sheep %>%
  filter(!(`Permanent Id` %in% sheep_dup$`Permanent Id` &
             `Capture Date` %in% sheep_dup$`Capture Date`))

# Now rbind the data back together & remove unneeded rows
simple_sheep2 <- rbind(simple_sheep,
                      sheep_dup) %>%
  select(-c(`Ear Tag`, CaseID, Received_Date, Report_Date, ELISA_Level,
            Lab_Notes, Notes, VetPak, PermID)) %>%
  rename(Myco_PCR = PCR_Result,
         Myco_Elisa = ELISA_Result)

# Now join this back to original sheep data
sheep <- left_join(sheep,
                   simple_sheep2,
                   by = c("Capture Date", "Permanent Id"))
  
# Now fill in the normal columns
sheep$`Movi PCR` <- ifelse(!is.na(sheep$Myco_PCR), sheep$Myco_PCR, sheep$`Movi PCR`)
sheep$`Movi Sero` <- ifelse(!is.na(sheep$Myco_Elisa), sheep$Myco_Elisa, sheep$`Movi Sero`)

# Remove unneeded columns
sheep <- sheep %>% select(-c(Myco_PCR, Myco_Elisa))

# Write the database
write_csv(sheep,
          "./Data/Capture Database updated 2023_CJP.csv")

