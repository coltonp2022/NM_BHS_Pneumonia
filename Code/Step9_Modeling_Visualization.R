# Modeling Visualization
# Colton Padilla
# 8/3/2023

# Last Updated: 7/9/2024

# Checked model fit using DHARMa
# Created model averaged figures for desert models

#-----------------#
# Housekeeping ####
#-----------------#

# Packages
library(tidyverse)
library(ggeffects)
library(ggpubr)
library(flextable)
library(DHARMa)
library(performance)

# Remove everything
rm(list = ls())

# Load the models
load("./Data/Top_Models.RData")

# Check correlation between PreP and PostP temp
cor(desert$PrePart_Temp, desert$PostPart_Temp) # r = 0.74
cor(desert$PrePart_Precip, desert$PostPart_Precip) # r = 0.27
cor(desert$PrePart_Temp, desert$PostPart_Precip) # r = -0.34

# Fit
res <- simulateResiduals(mod_rocky)
plot(res) # Fits
res <- simulateResiduals(mod_des)
plot(res) # Fits

# Odds ratios
#emmeans::emmeans(mod_rocky, type = "response")

#------------------------------------------------------------#
# Make a scatterplot of values for desert temp and precip ####
#------------------------------------------------------------#

# Plot
ggplot(desert,
       aes(x = PrePart_Precip,
           y = PrePart_Temp)) +
  geom_point(aes(fill = Lambs),
             pch = 21,
             size = 4)+
  theme_bw() +
  scale_fill_binned(breaks = c(0, 25, 50),
                    type = "viridis")

#------------#
# Rockies ####
#------------#

# Create model summary table
mod_sum <- summary(mod_rocky)
coefs <- data.frame(mod_sum$coefficients$cond) %>%
  mutate_if(is.numeric, round, 3) %>%
  rownames_to_column(var = "Coefficient") %>%
  rename(SE = `Std..Error`,
         Z = z.value,
         `P-value` = `Pr...z..`) %>%
  mutate(`P-value` = as.character(`P-value`)) %>%
  mutate(`P-value` = case_when(
    `P-value` < 0.05 & `P-value` > 0.01 ~ paste0(`P-value`, " ", "*"),
    `P-value` <= 0.01 & `P-value` > 0.001 ~ paste0(`P-value`, " ", "**"),
    `P-value` < 0.001  ~ paste0(`P-value`, " ", "***"),
    .default = `P-value`
  ))
coefs

# Reset some names
coefs$Coefficient <- gsub("\\)", "", coefs$Coefficient) %>%
  gsub("\\(", "", .) %>%
  gsub("StatusPost", "Disease Status (Post-Exposure)", .) %>%
  gsub("Est", "Abundance", .) %>%
  gsub("PrevLamb_Precip", "PLS Precipitation", .) %>%
  gsub("Lamb_Temp_ClassCold", "PLS Temperature (Cold)", .) %>%
  gsub("PrePart_Precip", "PreP Precipitation", .) %>%
  gsub("\\:", " * ", .)
tab <- flextable(coefs)

# Save the table
save_as_docx(tab,
             path = "./TablesFigures/Rocky_Model_Summary.docx")

# Set a theme
theme_set(theme_bw())

# Status
dat <- ggpredict(mod_rocky,
                 terms = "Status") %>%
  mutate(x = if_else(x == "Pre", "Pre-Exposure", "Post-Exposure"),
         x = factor(x, levels = c("Pre-Exposure", "Post-Exposure"))) %>%
  drop_na()
status <- ggplot(dat,
       aes(x,
           predicted)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = conf.low,
                    ymax = conf.high),
                linewidth = 0.7,
                width = 0.2) +
  labs(y = "",
       x = "Population Infection Status")

# Get data from ggpredict
dat <- ggpredict(mod_rocky,
                 terms = "PrePart_Precip [all]") %>%
  as.data.frame()

dat <- dat %>%
  mutate(PrePart_Precip = (x * 91.26726) + 120.933)

dat %>%
  summarise(range(predicted))
preprecip <- ggplot(dat,
       aes(PrePart_Precip,
           predicted)) +
  geom_line(linewidth = 0.7) +
  geom_ribbon(aes(ymin = conf.low,
                  ymax = conf.high),
              alpha = 0.3) +
  labs(y = "Lambs per 100 ewes",
       x = "Preparturition Precipitation (mm)") +
  theme(axis.title.x = element_text(margin = margin(t = 15)),
        axis.title.y = element_text(margin = margin(r = 15)))

# Save
ggsave(preprecip,
       filename = "./TablesFigures/PrePart_Precip_Graph.jpg",
       width = 7,
       height = 7,
       dpi = 1200)

# Plot the interaction
# Get data from ggpredict
dat <- ggpredict(mod_rocky,
                 terms = c("Est[all]",
                           "Status")) %>%
  as.data.frame()
dat <- dat %>%
  mutate(Population_Est = (x * 110.0805) + 146.245,
         group = if_else(group == "Pre", "Pre-Exposure", "Post-Exposure"))
rocky <- rocky %>%
  mutate(Population_Est = (Est * 110.0805) + 146.245,
         group = if_else(Status == "Pre", "Pre-Exposure", "Post-Exposure")) 

dat$group <- factor(dat$group, levels = c("Pre-Exposure", "Post-Exposure"))
dens <- ggplot(dat,
       aes(Population_Est,
           predicted,
           fill = group))  +
  geom_line(aes(color = group,
                linetype = group),
            linewidth = 0.7) +
  #geom_point(data = rocky,
  #           aes(x = Population_Est,
  #               y = Lambs,
  #               color = group),
  #           inherit.aes = F) +
  geom_ribbon(aes(ymin = conf.low,
                  ymax = conf.high),
              alpha = 0.4) +
  labs(y = "Lambs per 100 ewes",
       x = "Population Estimate",
       fill = "Infection Status",
       color = "Infection Status",
       linetype = "Infection Status") +
  scale_fill_manual(values = c("gray70","#E5140A")) +
  scale_color_manual(values = c("gray70","#E5140A")) +
  theme(legend.box = "horizontal",
        legend.position = "bottom",
        axis.title.x = element_text(margin = margin(t = 15)),
        axis.title.y = element_text(margin = margin(r = 15)))


# Plot the other interaction
dat <- ggpredict(mod_rocky,
                 terms = c("PrevLamb_Precip[all]",
                           "Lamb_Temp_Class"),
                 ci.lvl = 0.95) %>%
  as.data.frame()
dat <- dat %>%
  mutate(PrevLamb_Precip = (x * 77.88366) + 167.6656)
dat %>%
  group_by(group) %>%
  summarise(Range = range(PrevLamb_Precip))
dat %>% 
  group_by(group) %>%
  summarise(range(predicted))

rocky <- rocky %>%
  mutate(Unscale_Temp = (PrevLamb_Temp * attr(rocky$PrevLamb_Temp, "scaled:scale")) + 
           attr(rocky$PrevLamb_Temp, "scaled:center"))

rocky %>%
  group_by(PrevLamb_Temp_Class) %>%
  summarise(Range = range(Unscale_Temp))

prevprecip <- ggplot(dat,
       aes(PrevLamb_Precip,
           predicted,
           fill = group))  +
  geom_line(aes(color = group,
                linetype = group),
            linewidth = 0.7) +
  geom_ribbon(aes(ymin = conf.low,
                  ymax = conf.high),
              alpha = 0.4) +
  labs(y = "",
       x = "Previous Lamb Season Precipitation (mm)",
       fill = "Temperature Class",
       color = "Temperature Class",
       linetype = "Temperature Class") +
  scale_fill_manual(values = c("grey70", "deepskyblue")) +
  scale_color_manual(values = c("grey70", "deepskyblue")) +
  theme(legend.box = "horizontal",
        legend.position = "bottom",
        axis.title.x = element_text(margin = margin(t = 15)))

# GGarrange
figure <- ggarrange(dens, prevprecip, labels = c("A", "B"))
figure

# Save the figure
ggsave(filename = "./TablesFigures/Rocky_Model_Graphs.jpg",
       height = 6,
       width = 8.5,
       dpi = 500)

# Save overall figures
ggsave(filename = "./TablesFigures/DensityGraph.jpg",
       dens,
       height = 6, width = 6.5, dpi = 1200)
ggsave(filename = "./TablesFigures/TempClassGraph.jpg",
       prevprecip,
       height = 6, width = 6.5, dpi = 1200)


#------------#
# Deserts ####
#------------#

# Plot model
dat <- ggpredict(mod_des,
                 terms = "PrePart_Precip[all]") %>%
  as.data.frame()
dat <- dat %>%
  mutate(PrePart_Precip = (x * 34.04618) + 52.63677)
desplot <- desert %>%
  mutate(PrePart_Precip = (PrePart_Precip * 34.04618) + 52.63677) %>%
  mutate(label = str_c(Herd, Year))
precip <- ggplot(dat,
       aes(PrePart_Precip,
           predicted)) +
  geom_point(data = desplot,
             aes(x = PrePart_Precip,
                 y = Lambs)) +
  geom_line(linewidth = 0.7) +
  geom_ribbon(aes(ymin = conf.low,
                  ymax = conf.high),
              alpha = 0.3) +
  labs(y = "Lambs per 100 ewes",
       x = "Preparturition Precipitation (mm)") +
  theme(legend.box = "horizontal",
        legend.position = "bottom",
        axis.title.x = element_text(margin = margin(t = 15)),
        axis.title.y = element_text(margin = margin(r = 15)))

ggsave(precip,
       filename = "./TablesFigures/Desert_PrePart_Precip_Graph.jpg",
       width = 7,
       height = 7,
       dpi = 1200)


# Add labels
precip + geom_label(data = desplot, aes(x = PrePart_Precip,
                                        y = Lambs, 
                                        label = label))


# Prepart temp
dat <- ggpredict(mod_des,
                 terms = "PrePart_Temp[all]")
dat <- dat %>%
  mutate(PrePart_Temp = (x * 1.470264) + 10.26064)
temp <- ggplot(dat,
                 aes(PrePart_Temp,
                     predicted)) +
  geom_line(linewidth = 0.7) +
  geom_ribbon(aes(ymin = conf.low,
                  ymax = conf.high),
              alpha = 0.3) +
  labs(y = "",
       x = "Preparturition Temperature (C)") +
  theme(legend.box = "horizontal",
        legend.position = "bottom",
        axis.title.x = element_text(margin = margin(t = 15)),
        axis.title.y = element_text(margin = margin(r = 15)))

ggsave(temp,
       filename = "./TablesFigures/Desert_PrePart_Temp_Graph.jpg",
       width = 7,
       height = 7,
       dpi = 1200)

# Arrange together
# GGarrange
figure <- ggarrange(precip, temp, labels = c("A", "B"))
figure

# Save the figure
ggsave(filename = "./TablesFigures/Desert_Model_Graphs.jpg",
       height = 6,
       width = 8.5,
       dpi = 1200)

# Plot
dat <- ggpredict(mod_des,
                 terms = "Status")

dat$group <- factor(dat$group, levels = c("Pre-Exposure", "Post-Exposure"))
dens <- ggplot(dat,
               aes(Population_Est,
                   predicted,
                   fill = group))  +
  geom_line(aes(color = group,
                linetype = group),
            linewidth = 0.7) +
  geom_ribbon(aes(ymin = conf.low,
                  ymax = conf.high),
              alpha = 0.4) +
  labs(y = "",
       x = "Population Estimate",
       fill = NULL,
       color = NULL,
       linetype = NULL) +
  scale_fill_manual(values = c("gray70","#E5140A")) +
  scale_color_manual(values = c("gray70","#E5140A")) +
  theme(legend.direction = "horizontal",
        legend.position = c(0.5, 0.05),
        legend.background = element_blank(),
        axis.title.x = element_text(margin = margin(t = 15)),
        axis.title.y = element_text(margin = margin(r = 15)))

fig <- (ggarrange(precip, dens))

ggsave(filename = "./TablesFigures/Desert_Model_Graphs.jpg",
       height = 6,
       width = 8.5,
       dpi = 1200)


#-------------------------------#
# Plot Desert precipitations ####
#-------------------------------#

# unscale
unscale <- desert %>%
  mutate(across(11:18, ~ (. * attr(., "scaled:scale")) + 
                  attr(., "scaled:center")))

# Reshape
plotdat <- reshape2::melt(unscale, id.vars = c("Herd", "Lambs"),
               measure.vars = c("PrePart_Precip", "PostPart_Precip",
                                "Recovery_Precip", "PrevLamb_Precip"),
               vars = "Precip")

# Plot
ggplot(plotdat,
       aes(y = Lambs,
           x = value)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ variable, scales = "free")

#----------------------------------#
# Model averaged plots manually ####
#----------------------------------#

# Expand a grid 
grid <- expand.grid(
  Status = unique(desert$Status),
  Est = mean(desert$Est),
  PrePart_Temp = c(min(desert$PrePart_Temp), max(desert$PrePart_Temp)),
  PrePart_Precip = mean(desert$PrePart_Precip),
  Lion_REM = unique(desert$Lion_REM)
) %>%
  mutate(Status = if_else(Status == "Pre", 0, 1))

# Overwrite the full subset with the conditional subset
coefs <- avg_desert$coefficients[1,]

# Now write out the equation for mean
pred <- numeric()
for(i in 1:nrow(grid)){
  pred[i] <- coefs[1] + coefs[2] * grid$Status[i] + coefs[3] * grid$PrePart_Temp[i] +
    coefs[4] * grid$PrePart_Precip[i] + coefs[5] * (grid$PrePart_Precip[i] ^ 2) +
    coefs[6] * grid$Lion_REM[i] + coefs[7] * grid$Est[i] + coefs[8] * (grid$Est[i] * grid$Status[i])
}

# Bind it back and calculate
grid$pred <- exp(pred)
grid <- grid %>%
  group_by(PrePart_Temp) %>%
  summarise(pred = mean(pred))

# Expand a grid 
grid <- expand.grid(
  Status = unique(desert$Status),
  Est = mean(desert$Est),
  PrePart_Temp = mean(desert$PrePart_Temp),
  PrePart_Precip = c(min(desert$PrePart_Precip), max(desert$PrePart_Precip)),
  Lion_REM = unique(desert$Lion_REM)
) %>%
  mutate(Status = if_else(Status == "Pre", 0, 1))

# Overwrite the full subset with the conditional subset
coefs <- avg_desert$coefficients[1,]

# Now write out the equation for mean
pred <- numeric()
for(i in 1:nrow(grid)){
  pred[i] <- coefs[1] + coefs[2] * grid$Status[i] + coefs[3] * grid$PrePart_Temp[i] +
    coefs[4] * grid$PrePart_Precip[i] + coefs[5] * (grid$PrePart_Precip[i] ^ 2) +
    coefs[6] * grid$Lion_REM[i] + coefs[7] * grid$Est[i] + coefs[8] * (grid$Est[i] * grid$Status[i])
}

# Bind it back and calculate
grid$pred <- exp(pred)
grid <- grid %>%
  group_by(PrePart_Precip) %>%
  summarise(pred = mean(pred))
