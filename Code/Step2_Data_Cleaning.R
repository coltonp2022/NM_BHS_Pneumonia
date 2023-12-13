# Bighorn Sheep Mycoplasma Data Cleaning
# Colton Padilla
# 4/10/2023

# Last Updated: 4/10/2023

#-----------------#
# Housekeeping ####
#-----------------#

# Packages
library(tidyverse)
library(parasite)

# Bring in the database
myco <- read_csv("./Data/Clean/BHS_Myco_Data.csv")

# Hunter Harvest
harv <- read_csv("./Data/Hunter_Harvest_CR4282023.csv")

# Extra data
extra <- read_csv("./Data/Extra_Myco_Data.csv")

# If else
harv$Herd <- ifelse(harv$Herd == "WH_GH", "GH", harv$Herd)

# Pull myco in harv
myco_harv <- myco %>% filter(AnimalID %in% harv$AnimalID)
myco_non <-  myco %>% filter(!AnimalID %in% harv$AnimalID)

# Now overwrite age, herd, and sex of myco data in hunted herds
myco_harv$Age <- ifelse(myco_harv$AnimalID == harv$AnimalID, harv$Age, myco_harv$Age)
myco_harv$Herd <- ifelse(myco_harv$AnimalID == harv$AnimalID, harv$Herd, myco_harv$Herd)
myco_harv$Sex <- ifelse(myco_harv$AnimalID == harv$AnimalID, harv$Sex, myco_harv$Sex)

# Bind harv and non back together
myco <- rbind(myco_harv, myco_non)

#-----------------------------#
# Recode Sources and Herds ####
#-----------------------------#

# Check the unique names in these columns
unique(myco$Source)

# Recode
myco$Source <- ifelse(myco$Source == "CA", "CAB",
                      ifelse(myco$Source == "RR", "RRoc",
                             ifelse(myco$Source == "Fra", "FRA",
                                    ifelse(myco$Source == "na", NA, myco$Source))))
unique(myco$Source)

# Check uniques
unique(myco$Herd)

# Check
unique(myco$Herd)

# Recode
myco$Herd <- ifelse(myco$Herd == "Red River", "RRiv",
                    ifelse(myco$Herd == "San Andres", "SA",
                           ifelse(myco$Herd == "Ladron", "LA",
                                  ifelse(myco$Herd == "Gold Hill", "GH", 
                                         ifelse(myco$Herd == "FRA", "FC", myco$Herd)))))


#----------------------#
# Recode Sex column ####
#----------------------#

# Check
unique(myco$Sex)

# Recode
myco$Sex <- ifelse(!(myco$Sex == "M" | myco$Sex == "F"), NA, myco$Sex)

# Check
unique(myco$Sex)

#--------------------------#
# Recode the Age Column ####
#--------------------------#

# Remove extra characters
myco$Age <- myco$Age %>%
  gsub("\\+", "", .) %>%
  gsub("\\?", NA, .) %>%
  gsub("\\-Apr", "", .) 

# Check
unique(myco$Age)

# Recode lambs
myco$Age <- myco$Age %>%
  ifelse(grepl("lamb", .), "L", .) %>%
  ifelse(. == "Ad", "A", .) %>%
  ifelse(. == "adult", "A", .)

# Check
unique(myco$Age)

# Bin ages
Y <- c("1", "1.5", "2")
A <- c("2.5", "3", "3.5", "4", "4.5",
        "5", "5.5", "6", "6.5", "7", "7.5", "8",
        "8.5", "9", "9.5", "10", "10.5", "11")

# Now ifelse
myco$Age2 <- ifelse(myco$Age %in% Y, "Y", 
                    ifelse(myco$Age %in% A, "A", myco$Age))

# Check
unique(myco$Age2)

# Recode the last weird ones
myco$Age2 <- myco$Age2 %>%
  gsub("4/6", "A", .) %>%
  gsub("5/7", "A", .) %>%
  gsub("1.5/4", NA, .)

# Check
unique(myco$Age2)

#--------------------------------------#
# Get a database of unknown animals ####
#--------------------------------------#

# Pick animals where HERD is NA
sort <- myco %>% filter(is.na(Herd))

# Get hunter 
hunt <- sort[grepl("Hunter", sort$Notes),]
hunt1 <- sort[grepl("hunter", sort$Notes),]
hunt <- rbind(hunt, hunt1)

# Remove from sort
sort <- sort %>% filter(!AnimalID %in% hunt$AnimalID)

#------------------------------------------------#
# Create Extra Myco Data from unknown samples ####
#------------------------------------------------#

# Wrangle the data to equal the myco data
extra1 <- extra %>%
  group_by(Herd, Year, Test) %>%
  summarise(Result = c(rep(1, Positive),
                       rep(0, Total - Positive)))

pcr <- extra1 %>% 
  filter(Test == "PCR") %>%
  rename(PCR_Result = Result) %>%
  ungroup() %>%
  select(-Test)

elisa <- extra1 %>% 
  filter(Test == "ELISA" | Test == "PH") %>%
  rename(ELISA_Result = Result) %>%
  ungroup() %>%
  select(-Test)

# Use rbind.fill to bind them together
myco <- do.call(plyr::rbind.fill, list(myco, pcr, elisa))

#--------------------#
# Write databases ####
#--------------------#

# write myco
write_csv(myco,
          "./Data/Clean/BHS_Myco_Data_Clean.csv")

# write hunter
write_csv(hunt,
          "./Data/Hunter_Harvest.csv")

# Write other unknowns
write_csv(sort,
          "./Data/Other_Unknowns.csv")
