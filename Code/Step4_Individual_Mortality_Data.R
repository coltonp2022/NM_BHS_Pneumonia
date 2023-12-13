# Add Individual Mortality Data
# Colton Padilla
# 5/15/2023

# Last Updated: 5/15/2023

#-----------------#
# Housekeeping ####
#-----------------#

# Packages
library(tidyverse)
library(openxlsx)
library(survival)
library(ggsurvfit)
library(gtsummary)

# Load the myco data
myco <- read_csv("./Data/Clean/BHS_Myco_Analysis_Data.csv")

# Load mortality data
mort <- read.xlsx("./Data/Bighorn GPS Collars March 2022_3_28_23.xlsx", sheet = 1)

#----------------------------#
# Begin Data Manipulation ####
#----------------------------#

# Get only needed columns in mort
mort <- mort %>%
  dplyr::select(BighornID, BeginDate, EndDate, EndCause, ReleaseArea, MortalityCause)

# Join these to the myco data
myco <- left_join(myco, mort %>% select(BighornID, EndDate, EndCause, MortalityCause), by = c("PermID" = "BighornID"))

# Remove harvested individuals
myco <- myco %>%
  filter(!grepl("harvest", MortalityCause),
         !grepl("Harvest", MortalityCause),
         !grepl("Capture", MortalityCause))

# Recode the date column
myco$EndDate <- convertToDate(myco$EndDate) %>%
  as.Date(.)

# Now Filter only to collared individuals
myco <- myco %>%
  filter(PermID %in% mort$BighornID) %>%
  drop_na(PermID) %>%
  select(PCR_Result, ELISA_Result, PermID, Herd, Wt_kg, Age, Sex,
         Cap_Date, EndDate, EndCause, MortalityCause) %>%
  mutate(Age = as.numeric(Age))

# Add sex for one individual
myco$Sex <- ifelse(is.na(myco$Sex), "F", myco$Sex)

# Create a database for population survival
pop_surv <- mort

# Now get a survival data set
surv <- myco %>%
  select(PermID,
         ELISA_Result,
         PCR_Result,
         Cap_Date,
         EndDate,
         EndCause,
         MortalityCause,
         Herd,
         Sex,
         Age) 

# Create a status column
surv$status <- ifelse((grepl("transmitting", surv$EndCause) | is.na(surv$EndCause) |
                         grepl("breaks", surv$EndCause) |
                         grepl("Re-collared", surv$EndCause)),
                      0,
                      1)

# Fix Ladron caputure dates
surv$Cap_Date <- gsub("0020", "2020", surv$Cap_Date) %>%
  as.Date(., format = "%Y-%m-%d")

# Make all end dates today if na
surv$EndDate <- ifelse(is.na(surv$EndDate),
                       "2023-03-28",
                       format.Date(surv$EndDate, format = "%Y-%m-%d")) %>%
  as.Date(., format = "%Y-%m-%d")

# Days survived
surv$Days <- as.character(surv$EndDate - surv$Cap_Date) %>%
  as.numeric()

# Sero
Sero <- surv %>% 
  filter(PCR_Result != 1,
         ELISA_Result != 0.5)

# PCR
PCR <- surv

#--------------------------------------------#
# Get individuals in positive populations ####
#--------------------------------------------#

# Fix Date columns
pop_surv$BeginDate <- convertToDate(pop_surv$BeginDate) %>%
  as.Date(.)
pop_surv$EndDate <- convertToDate(pop_surv$EndDate) %>%
  as.Date(.)

# Create a status column
pop_surv$status <- ifelse(!pop_surv$EndCause == "Mortality" | is.na(pop_surv$EndCause),
                      0,
                      1)

# Filter to Ladron and Fras
pop_surv <- pop_surv %>%
  filter(ReleaseArea %in% c("Ladron", "Fra Cristobal"))

# Make all end dates today if na
pop_surv$EndDate <- ifelse(is.na(pop_surv$EndDate),
                       "2023-03-28",
                       format.Date(pop_surv$EndDate, format = "%Y-%m-%d")) %>%
  as.Date(., format = "%Y-%m-%d")

# Days survived
pop_surv$Days <- as.character(pop_surv$EndDate - pop_surv$BeginDate) %>%
  as.numeric()

# Get a start year and end year
pop_surv$StartYear <- year(pop_surv$BeginDate)
pop_surv$EndYear <- year(pop_surv$EndDate)

# Now create a database of first positives
pop_disease <- data.frame(
  Herd = c("FC",  "LA"),
  Year = c(2017, 2017)
)

# Recode release area
pop_surv$ReleaseArea <- ifelse(pop_surv$ReleaseArea == "Ladron",
                               "LA",
                               "FC")

# Unique populations
LA <- 2017
FC <- 2017

# Initialize the column
pop_surv$PopStatus <- NA


# Get ladron data
Ladron <- pop_surv %>% filter(ReleaseArea == "LA")

# Get FC data 
Fras <- pop_surv %>% filter(ReleaseArea == "FC")

# Ladron popstatus
Ladron$PopStatus <- ifelse(between(Ladron$StartYear,
                                   2017,
                                   2023),
                           "Post",
                           ifelse(between(Ladron$EndYear,
                                          2017,
                                          2023),
                                  "Post", "Pre"))

# Ladron popstatus
Fras$PopStatus <- ifelse(between(Fras$StartYear,
                                   2017,
                                   2023),
                           "Post",
                           ifelse(between(Fras$EndYear,
                                          2017,
                                          2023),
                                  "Post", "Pre"))


#------------------------------------#
# Kaplan-Meier for seropositivity ####
#------------------------------------#

# Reset factor
Sero$ELISA_Result <- ifelse(Sero$ELISA_Result == 0,
                            "Negative",
                            "Positive")

# Run the kaplan-meier
n <- survfit(Surv(Days, status) ~ ELISA_Result, 
        data = Sero)

# Plot the kaplan-meier
survfit2(Surv(Days, status) ~ ELISA_Result, 
         data = Sero) %>%
  ggsurvfit() +
  labs(color = "Serology",
       y = "Survival probability",
       x = "Days after capture") +
  scale_color_manual(values = c("gray70", "#E5140A")) +
  add_censor_mark() +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))

# Save
ggsave("./TablesFigures/Serological_Survival.jpg",
       dpi = 600,
       height = 5,
       width = 6.5)

# Presentation
ggsave("./TablesFigures/Serological_Survival_Presentation.jpg",
       dpi = 600,
       height = 7,
       width = 7)

# Check ages
ggplot(Sero, 
       aes(x = Age,
           fill = ELISA_Result)) + 
  geom_density(alpha = 0.4) + 
  theme_bw() + 
  scale_fill_manual(values = c("gray70", "#E5140A"))


# Plot with CI
survfit2(Surv(Days, status) ~ ELISA_Result, 
         data = Sero) %>%
  ggsurvfit() +
  labs(linetype = "Seropositivity\nat Capture",
       y = "Survival Probability",
       x = "Days after Capture") +
  scale_fill_manual(values = c("#E69F00", "#56B4E9")) +
  scale_color_manual(values = c("#E69F00", "#56B4E9")) +
  add_confidence_interval()

# Get yearly survival estimates
summary(survfit(Surv(Days, status) ~ ELISA_Result, 
                data = Sero),
        times = 850)

options(scipen = 999)

# Differences
survdiff(Surv(Days, status) ~ ELISA_Result,
         data = Sero)

# Regression
coxph(Surv(Days, status) ~ ELISA_Result * Sex,
      data = Sero) %>% tbl_regression(exp = T)

# Run the K-M for PCR
survfit(Surv(Days, status) ~ PCR_Result, 
        data = PCR)

# Plot the kaplan-meier
survfit2(Surv(Days, status) ~ PCR_Result, 
         data = PCR) %>%
  tidy_survfit() %>%
  ggplot(aes(x = time,
             y = estimate,
             linetype = strata)) +
  geom_step(linewidth = 1) +
  labs(linetype = "Active Infection\nat Capture",
       y = "Survival Probability",
       x = "Days") 


# Regression
coxph(Surv(Days, status) ~ PCR_Result,
      data = PCR) %>% tbl_regression(exp = T)

#-----------------#
# Test by Year ####
#-----------------#

# Create a Year column
surv$Year <- year(surv$Cap_Date)

# Loop 
survfit2(Surv(Days, status) ~ ELISA_Result, 
         data = surv %>% filter(PCR_Result != 1,
                                  Year == 2017)) %>%
  tidy_survfit() %>%
  ggplot(aes(x = time,
               y = estimate,
               linetype = strata)) +
  geom_step(linewidth = 1) +
  labs(linetype = "Status\nat Capture",
       y = paste0("2022", "Survival Probability"),
       x = "Days")

# Add CIs
survfit2(Surv(Days, status) ~ ELISA_Result, 
         data = surv %>% filter(Year == 2022,
                                PCR_Result != 1)) %>%
  ggsurvfit() +
  labs(linetype = "Seropositivity\nat Capture",
       y = "Survival Probability",
       x = "Days after Capture") +
  scale_fill_manual(values = c("#E69F00", "#56B4E9")) +
  scale_color_manual(values = c("#E69F00", "#56B4E9")) + 
  add_confidence_interval()

#-----------------------------#
# Kaplan-Meier Sacramentos ####
#-----------------------------#

# Filter to the SACs
Sac <- Sero %>% filter(Herd == "SAC")

# Run the kaplan-meier
survfit(Surv(Days, status) ~ ELISA_Result, 
             data = Sac)

# Plot the kaplan-meier
survfit2(Surv(Days, status) ~ ELISA_Result, 
         data = Sac) %>%
  tidy_survfit() %>%
  ggplot(aes(x = time,
             y = estimate,
             color = strata)) +
  geom_step(linewidth = 1) +
  labs(color = "Serology",
       y = "Survival probability",
       x = "Days after capture") +
  scale_color_manual(values = c("gray70", "#E5140A")) +
  theme(legend.position = "bottom")

# Differences
survdiff(Surv(Days, status) ~ ELISA_Result,
         data = Sac)

# Save
ggsave("./TablesFigures/Sacramento_Survival.jpg",
       dpi = 600,
       height = 7,
       width = 7)

# Plot with CI
survfit2(Surv(Days, status) ~ ELISA_Result, 
         data = Sac) %>%
  ggsurvfit() +
  labs(linetype = "Seropositivity\nat Capture",
       y = "Survival Probability",
       x = "Days after Capture") +
  scale_fill_manual(values = c("#E69F00", "#56B4E9")) +
  scale_color_manual(values = c("#E69F00", "#56B4E9")) +
  add_confidence_interval()

# Get yearly survival estimates
summary(survfit(Surv(Days, status) ~ ELISA_Result, 
                data = Sac),
        times = 365.25 * 2)

# Differences
survdiff(Surv(Days, status) ~ ELISA_Result,
         data = Sac)


# View each population statuses
check <- Sero %>% 
  group_by(Herd, ELISA_Result)  %>%
  summarise(n = n(),
            sum = sum(status),
            perc_mort = sum / n) %>%
  filter(n() > 1) %>%
  mutate(perc_mort = if_else(perc_mort == 0, 0.01, perc_mort))

ggplot(check,
       aes(x = Herd,
           y = perc_mort,
           fill = ELISA_Result)) +  
  geom_bar(stat = "identity",
           position = "dodge")

Sero %>% filter(Herd == "RRiv") %>% count(ELISA_Result, EndCause) %>%
  group_by(ELISA_Result) %>% mutate(perc = n / sum(n))

Sero %>% count(MortalityCause, ELISA_Result)
Sero %>% count(ELISA_Result)
