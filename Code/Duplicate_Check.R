# Check for duplicates in BHS Disease Data
# Colton Padilla
# 4/4/2023

# Last Ran: 4/7/2023

# Packages
library(tidyverse)

# Bring in the data set
myco <- read_csv("./Data/BHS_Myco_Testing_4_7_2023.csv")

# Duplicate rows check
table(duplicated(myco))

# Duplicate ID check
table(duplicated(myco$AnimalID))

# Pull duplicates to check in data sheet
myco[duplicated(myco$AnimalID) | duplicated(myco$AnimalID, fromLast = T),]
  # P1014 - Two Different results from different cases
  # P762 - Two Different results from the same case
  # R578 - Different ELISA Results from same case. Different Animal????

