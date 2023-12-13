# Cleaning Population Sizes
# Colton Padilla
# 7/11/2023

# Last Updated: 7/11/2023

#-----------------#
# Housekeeping ####
#-----------------#

# Packages
library(tidyverse)

# Data
pop <- read_csv("./Data/Pop_Sizes.csv")

#--------------------#
# Clean Pop sizes ####
#--------------------#

# Remove non numerics from column
pop$Year <- gsub("[^0-9.-]", "", pop$Date) %>%
  gsub("\\-", "", .) %>%
  gsub("\\.", "", .)

# Add 20 if the number is only two digits
pop <- pop %>%
  mutate(Year = if_else(str_detect(Year, "^.{2}$"), paste0("20", Year), Year))

# Now Move BH and LH to herd column
pop$Herd <- ifelse(grepl("BH", pop$Date), "BH",
                   ifelse(grepl("LH", pop$Date), "LH", pop$Herd))

# Now toss those strings and numbers
pop$Month <- gsub("[0-9]+", "", pop$Date) %>%
  gsub("BH", "", .) %>%
  gsub("LH", "", .) %>%
  gsub("\\-", "", .) %>%
  gsub("'", "", .) %>%
  gsub("\\_", "", .) %>%
  gsub("\\.", "", .) %>%
  gsub(" ", "", .)

# Add NAs
pop$Month <- ifelse(pop$Month == "", NA, pop$Month)

# Check
unique(pop$Month)

# Now recode months
pop$Month <- ifelse(pop$Month == "Jul", "JUL",
                    ifelse(pop$Month == "Aug", "AUG",
                           ifelse(pop$Month == "Mar", "MAR",
                                  ifelse(pop$Month == "Sep", "SEP", 
                                         ifelse(pop$Month == "July", "JUL",
                                                ifelse(pop$Month == "Oct", "OCT", pop$Month))))))

# Check
unique(pop$Month)

# Recode more
pop$Month <- ifelse(pop$Month == "Nov", "NOV",
                    ifelse(pop$Month == "Jun", "JUN",
                           ifelse(pop$Month == "Dec", "DEC",
                                  ifelse(pop$Month == "May", "MAY",
                                         ifelse(pop$Month == "Jan", "JAN",
                                                ifelse(pop$Month == "Sept", "SEP", pop$Month))))))

# Check
unique(pop$Month)

# One more
pop$Month <- ifelse(pop$Month == "June", "JUN",
                    ifelse(pop$Month == "Feb", "FEB",
                           ifelse(pop$Month == "JULY", "JUL",
                                  ifelse(pop$Month == "Apr", "APR",
                                         ifelse(pop$Month == "jun", "JUN",
                                                ifelse(pop$Month == "April", "APR", pop$Month))))))

# Check
unique(pop$Month)

# Average estimates in the same column
pop <- pop %>% 
  drop_na(Est) %>%
  rowwise() %>%
  mutate(Est = round(mean(as.numeric(unlist(strsplit(Est, "-"))))), 0) %>%
  select(-c(Date, `0`))

# Keep the second value
pop <- pop %>%
  group_by(Herd, Year) %>%
  filter(Est == max(Est))

# Write out the cleaned population sizes
write_csv(pop,
          "./Data/Clean/Population_Sizes.csv")
