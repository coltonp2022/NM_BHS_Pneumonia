# Cleaning Lamb Ewe Ratio Data
# Colton Padilla
# 5/15/2023

# Last Updated: 6/6/2023

#-----------------#
# Housekeeping ####
#-----------------#

# Packages
library(tidyverse)
library(openxlsx)
library(parasite)

# Load in the data
lamb <- read.xlsx("./Data/LambEweRatios.xlsx", sheet = 2)

# Load in extra lamb data found 7/1/2023
extra <- read_csv("./Data/Extra_Ratios.csv")

# Calculate extra lamb ewe ratios
extra$Lambs <- round((100 / extra$Ewes) * extra$Lambs)

# Remove lambs
extra <- extra %>% select(-Ewes)

# Load myco data
myco <- read_csv("./Data/Clean/BHS_Myco_Analysis_Data.csv")

# Build timeline table
time <- data.frame(
  Herd = c("SA", "FC", "CAB", "LA", "SAC", "RGG", "RRiv", "WH", "CUL"),
  Year = c(2007, 2017, 2019, 2017, 2018, 2020, 2022, 2022, 2022)
)

#-----------------------#
# Begin Manipulation ####
#-----------------------#

# Remove other numbers & characters
lamb$Lambs <- lamb$Lambs %>% 
  gsub(".*100:", "", .) %>%
  gsub("[^[:digit:]]", "", .)

# Remove years with no ratio
lamb <- lamb %>%
  drop_na(Lambs)

# Drop rows with nothing
lamb <- lamb[!lamb$Lambs == "", ]

# Now code BH and LH
lamb$Herd <- ifelse(grepl("BH", lamb$Year), "BH",
                    ifelse(grepl("LH", lamb$Year), "LH", lamb$Herd))

# Now make the rest the overall Hatchets
lamb$Herd <- ifelse(is.na(lamb$Herd), "HAT", lamb$Herd)

# Now remove those prefixes
lamb$Year <- lamb$Year %>% 
  gsub("BH", "", .) %>%
  gsub("LH", "", .)

# Now move months to other column
lamb$Month <- lamb$Year %>%
  gsub("[[:digit:]]+", "", .) %>%
  gsub("\\.", "", .) %>%
  gsub("\\_", "", .) %>%
  gsub("\\-", "", .) %>%
  gsub("\\'", "", .) %>%
  gsub(" x", "", .) %>%
  gsub(" y", "", .) %>%
  gsub(" c", "", .) %>%
  gsub(" d", "", .) %>%
  gsub(" h", "", .) %>%
  gsub(" g", "", .) %>%
  gsub("Sepm", "Sep", .) %>%
  gsub("Sept", "Sep", .) %>%
  gsub(" ", "", .) %>%
  gsub("\n", "", .) %>%
  toupper(.) %>%
  gsub("JULY", "JUL", .) %>%
  gsub("APRIL", "APR", .) %>%
  gsub("JUNE", "JUN", .)

# Check
unique(lamb$Month)

# Now get rid of characters in the year
lamb$Year <- lamb$Year %>%
  gsub("[^[:digit:]]", "", .)

# Now remove and then readd 20
lamb$Year <- lamb$Year %>%
  gsub("2014", "14", .) %>%
  gsub("2018", "18", .) %>%
  gsub("2019", "19", .) %>%
  gsub("2022", "22", .) %>%
  paste0("20", .)
unique(lamb$Year) # Check

# Make Year numeric
lamb$Year <- as.numeric(lamb$Year)
lamb$Lambs <- as.numeric(lamb$Lambs)

# Bind the extra data
lamb <- rbind(lamb, extra)

#---------------------------#
# Bring in the myco data ####
#---------------------------#

# Summarise the data by herd and year
myco <- myco %>%
  group_by(Herd, Year) %>%
  filter(!PCR_Result == 0.5) %>%
  summarise(Positives = sum(PCR_Result),
            Samp_Size = n())

# Join this data
lamb <- left_join(lamb, myco, by = c("Herd", "Year"))

# Code up years with PCR positives
lamb$Status <- ifelse(lamb$Positives > 0, "Post", NA)

# Get positive years
pos <- lamb %>% filter(Status == "Post") %>%
  distinct(Herd, Year)

# Rbind pos and time and then get distinct
time <- rbind(time, pos) %>%
  distinct()

# Make subsequent years post exposure
for(i in 1:nrow(time)){
  dat <- time[i, ]
  
  lamb$Status <- ifelse((lamb$Herd == dat$Herd &
                           between(lamb$Year,
                                   dat$Year,
                                   dat$Year + 500)),
                        "Post",
                        lamb$Status)
}

# Add Latir Positivity
lamb$Status <- ifelse((lamb$Herd == "LAT" & 
                         lamb$Year == 2022),
                      "Post",
                      lamb$Status)

# Make all others pre
lamb$Status <- ifelse(is.na(lamb$Status),
                      "Pre",
                      lamb$Status)

# Use the two sample ttest to check for patterns
boot2samp_t_test(lamb, 
                 "Lambs",
                 "Status",
                 "two.sided")

# Now create a 3 year post column
lamb$Status3 <- NA
for(i in 1:nrow(time)){
  dat <- time[i, ]
  
  lamb$Status3 <- ifelse((lamb$Herd == dat$Herd &
                           between(lamb$Year,
                                   dat$Year,
                                   dat$Year + 3)),
                        "Post",
                        lamb$Status3)
}

# Add latir positivity
lamb$Status3 <- ifelse(lamb$Herd == "LAT",
                       lamb$Status,
                       lamb$Status3)

# Make all others pre
lamb$Status3 <- ifelse(is.na(lamb$Status3),
                       "Pre",
                       lamb$Status3)

# Make a Latir Specific column to test for 2016 positivity
lamb$StatusLat <- ifelse((lamb$Herd == "LAT" &
                            between(lamb$Year,
                                    2016,
                                    2023)),
                         "Post", 
                         lamb$Status)

#---------------------------#
# Add in Rocky vs Desert ####
#---------------------------#

# Get unique herds
unique(lamb$Herd)

# Get rocky
rocky <- c("PEC", "WH", "SFR", "TC", "DC", "MZ", "JZ", "RGG", "LAT", "CUL")

# If else
lamb$Subspecies <- ifelse(lamb$Herd %in% rocky,
                          "Rocky",
                          "Desert")

#-----------------#
# Write the DF ####
#-----------------#

# write
write_csv(lamb,
          "./Data/Clean/Lamb-Ewe_Ratios_Clean.csv")
  


