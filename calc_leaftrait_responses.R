# Code for modelling idealized leaf trait responses to environmental changes as used in:
# Temporal constraints on leaf-level trait plasticity for next-generation land surface models
# Odé at al. (yr), in preparation

#### variable key ####

# P-model variables:
#tg_c: acclimated temperature (degrees Celsius)
#vpdo: vapor pressure deficit at sea level (kPa)
#cao: atmospheric CO2 at sea level (umol mol-1)
#paro: photosynthetically active radiation at sea level (µmol m-2 s-1)
#gs: stomatal conductance to water vapour
#Ci: intercellular CO2 concentration (ppm)
#Chi: χ(optimal) (unitless)
#Vcmax: maximum carboxylation rate (µmol·m-2·s-1)
#ALEAF: net assimilation rate (µmol·m-2·s-1)

# Photosyn variables:
#Ci: intercellular CO2 concentration (ppm)
#VPD: vapour pressure deficit (kPa)
#Tleaf: leaf temperature (degrees Celsius)
#Ca: atmospheric CO2 concentration (ppm)
#PPFD: photosynthetic photon flux density (µmol m-2·s-1)
#GS: stomatal conductance to water vapour 

# Load libraries
library(R.utils)
library(plantecophys)
library(ggplot2)
library(tidyverse)

# Load necessary functions using path to own working directory
sourceDirectory("functions", modifiedOnly = FALSE)
source("calc_optimal_vcmax.R")

#### CO2 scenario ####
#### Part 1: framework scenario ####

# Choose environmental start- and end variables #
cao_start<-400 
cao_end<- 800 
vpdo_start <- 1
vpdo_end <- 1
paro_start <- 800 
paro_end <- 800
t_start<-25
t_end<-25

# Calculate start (optimal) trait values#
rpmodel_step0<-data.frame()
rpmodel_step0<-calc_optimal_vcmax(cao=cao_start, par=paro_start, vpdo=vpdo_start, tg_c=t_start) #calculate start values
vcmax_start<-rpmodel_step0$vcmax[1] #Store start Vcmax value
chi_start<-rpmodel_step0$chi[1] #Store start chi value

# Calculate end (optimal) trait values #
step_end_rpmodel <-calc_optimal_vcmax(cao=cao_end, par=paro_end, vpdo=vpdo_end, tg_c=t_end) #calculate end values
vcmax_end<-step_end_rpmodel$vcmax[1] #Store end Vcmax value
chi_end<-step_end_rpmodel$chi[1] #Store end chi value

# Step 0: Start conditions
step0 <- Photosyn(Ca=cao_start, PPFD=paro_start, Tleaf=t_start, VPD=vpdo_start, Vcmax=vcmax_start)
ci_start<-chi_start*step0$Ca
step0 <- Photosyn(Ci=ci_start, Ca=cao_start, PPFD=paro_start, Tleaf=t_start, VPD=vpdo_start, Vcmax=vcmax_start, gsmodel="BBOpti")
gs_start<-step0$GS

# Step 0-1: Instantaneous response 
step1 <- Photosyn(GS=gs_start, Ca=cao_end, Vcmax=vcmax_start, VPD=vpdo_end, PPFD=paro_end, Tleaf=t_end, gsmodel="BBOpti")

# Step 1-2: Stomatal aperture response
step2 <- Photosyn(cao_end, Vcmax=vcmax_start, VPD=vpdo_end, PPFD=paro_end, Tleaf=t_end, gsmodel="BBOpti")

# Step 2-3: Biochemical acclimation response
step3 <- Photosyn(Ca=cao_end, VPD=vpdo_end, PPFD=paro_end, Vcmax=vcmax_end, Tleaf=t_end, gsmodel="BBOpti")
ci_end<-chi_end*step3$Ca
step3 <- Photosyn(Ci=ci_end, Ca=cao_end, VPD=vpdo_end, PPFD=paro_end, Vcmax=vcmax_end, Tleaf=t_end, gsmodel="BBOpti")

# Step 4: Developmental response (same as step 3, but we will adjust gsmax in line 89)
step4 <- step3

# Combine all steps in one dataframe
combined <- rbind(step0, step1, step2, step3, step4)

# Add Chi and Vcmax columns to dataframe
combined$Chi <- c(chi_start, rep(chi_end, 4))
combined$Vcmax <- c(vcmax_start, vcmax_start, vcmax_start, vcmax_end, vcmax_end)

# Add gsmax values to dataframe
gsmax_start <- step0$GS / 0.25
gsmax_end <- step3$GS / 0.25
combined$gsmax <- c(gsmax_start, gsmax_start, gsmax_start, gsmax_start, gsmax_end)

# Calculate x- and y-axes
combined$xaxis <- combined$GS / combined$gsmax
combined$yaxis <- (combined$Ci / c(cao_start, rep(cao_end, 4))) / combined$Chi

# Label columns as steps and reorder
combined$Step <- c("0", "1", "2", "3", "4")
combined <- combined[, c("Step", "Chi","Ci","GS","Vcmax","gsmax", "xaxis", "yaxis", names(combined)[!names(combined) %in% c("Step", "GS", "Ci", "Chi", "Vcmax", "gsmax", "xaxis", "yaxis")])]

# Plot the scenario #
ggplot(combined, aes(x=xaxis, y=yaxis, label=Step)) +
  geom_point() +
  geom_path(linetype='solid', linewidth=1.5) +
  theme_classic(base_size = 20) +
  xlab("gs/gsmax") + ylab("ci:ca/chi") +
  ggtitle("Scenario title") +
  geom_hline(yintercept=1, color='Blue', linewidth=1, linetype='dashed') +
  geom_vline(xintercept=0.25, color='Blue', linewidth=1, linetype='dashed') +
  ylim(0.90,1.2) + xlim(0.12,0.3) +
  geom_text(hjust=-0.2, vjust=-0.5, show.legend = FALSE)  # Add this line to show labels

#### Part 2: timelapse of normalized leaf trait values ####

# Normalize trait values to the maximum value per trait #
combined <- combined %>%
  mutate(GS_normalized = GS / max(GS, na.rm = TRUE))
combined <- combined %>%
  mutate(ALEAF_normalized = ALEAF / max(ALEAF, na.rm = TRUE))
combined <- combined %>%
  mutate(Ci_normalized = Ci / max(Ci, na.rm = TRUE))
combined <- combined %>%
  mutate(Vcmax_normalized = Vcmax / max(Vcmax, na.rm = TRUE))

# Create the data frame with normalized values and pivot to correct format #
timelapse_dataframe <- combined %>%
  select(Step, GS_normalized, ALEAF_normalized, Ci_normalized, Vcmax_normalized) %>%
  pivot_longer(cols = c(GS_normalized, ALEAF_normalized, Ci_normalized, Vcmax_normalized),
               names_to = "Variable",
               values_to = "Normalized_value") %>%
  mutate(Variable = sub("_normalized", "", Variable))

# Create a plot for the timelapse grouping the values per trait value # 
timelapse_dataframe<-ggplot(timelapse_dataframe, aes(x=Step, y=Normalized_value, color=Variable)) +
  geom_line(linewidth = 1, aes(group=Variable)) +
  geom_point(size=3)+
  theme_bw()+
  ggtitle("Timelapse ")+
  theme_classic(base_size = 15)

print(timelapse_dataframe)

#### VPD scenario ####

#### Part 1: framework scenario ####

# Choose environmental start- and end variables #
cao_start<-400 
cao_end<- 400 
vpdo_start <- 1
vpdo_end <- 2
paro_start <- 800 
paro_end <- 800
t_start<-25
t_end<-25

# Calculate start (optimal) trait values#
rpmodel_step0<-data.frame()
rpmodel_step0<-calc_optimal_vcmax(cao=cao_start, par=paro_start, vpdo=vpdo_start, tg_c=t_start) #calculate start values
vcmax_start<-rpmodel_step0$vcmax[1] #Store start Vcmax value
chi_start<-rpmodel_step0$chi[1] #Store start chi value

# Calculate end (optimal) trait values #
step_end_rpmodel <-calc_optimal_vcmax(cao=cao_end, par=paro_end, vpdo=vpdo_end, tg_c=t_end) #calculate end values
vcmax_end<-step_end_rpmodel$vcmax[1] #Store end Vcmax value
chi_end<-step_end_rpmodel$chi[1] #Store end chi value

# Step 0: Start conditions
step0 <- Photosyn(Ca=cao_start, PPFD=paro_start, Tleaf=t_start, VPD=vpdo_start, Vcmax=vcmax_start)
ci_start<-chi_start*step0$Ca
step0 <- Photosyn(Ci=ci_start, Ca=cao_start, PPFD=paro_start, Tleaf=t_start, VPD=vpdo_start, Vcmax=vcmax_start, gsmodel="BBOpti")
gs_start<-step0$GS

# Step 0-1: Instantaneous response 
step1 <- Photosyn(GS=gs_start, Ca=cao_end, Vcmax=vcmax_start, VPD=vpdo_end, PPFD=paro_end, Tleaf=t_end, gsmodel="BBOpti")

# Step 1-2: Stomatal aperture response
step2 <- Photosyn(cao_end, Vcmax=vcmax_start, VPD=vpdo_end, PPFD=paro_end, Tleaf=t_end, gsmodel="BBOpti")

# Step 2-3: Biochemical acclimation response
step3 <- Photosyn(Ca=cao_end, VPD=vpdo_end, PPFD=paro_end, Vcmax=vcmax_end, Tleaf=t_end, gsmodel="BBOpti")
ci_end<-chi_end*step3$Ca
step3 <- Photosyn(Ci=ci_end, Ca=cao_end, VPD=vpdo_end, PPFD=paro_end, Vcmax=vcmax_end, Tleaf=t_end, gsmodel="BBOpti")

# Step 4: Developmental response (same as step 3, but we will adjust gsmax in line 89)
step4 <- step3

# Combine all steps in one dataframe
combined <- rbind(step0, step1, step2, step3, step4)

# Add Chi and Vcmax columns to dataframe
combined$Chi <- c(chi_start, rep(chi_end, 4))
combined$Vcmax <- c(vcmax_start, vcmax_start, vcmax_start, vcmax_end, vcmax_end)

# Add gsmax values to dataframe
gsmax_start <- step0$GS / 0.25
gsmax_end <- step3$GS / 0.25
combined$gsmax <- c(gsmax_start, gsmax_start, gsmax_start, gsmax_start, gsmax_end)

# Calculate x- and y-axes
combined$xaxis <- combined$GS / combined$gsmax
combined$yaxis <- (combined$Ci / c(cao_start, rep(cao_end, 4))) / combined$Chi

# Label columns as steps and reorder
combined$Step <- c("0", "1", "2", "3", "4")
combined <- combined[, c("Step", "Chi","Ci","GS","Vcmax","gsmax", "xaxis", "yaxis", names(combined)[!names(combined) %in% c("Step", "GS", "Ci", "Chi", "Vcmax", "gsmax", "xaxis", "yaxis")])]

# Plot the scenario #
ggplot(combined, aes(x=xaxis, y=yaxis, label=Step)) +
  geom_point() +
  geom_path(linetype='solid', linewidth=1.5) +
  theme_classic(base_size = 20) +
  xlab("gs/gsmax") + ylab("ci:ca/chi") +
  ggtitle("Scenario title") +
  geom_hline(yintercept=1, color='Blue', linewidth=1, linetype='dashed') +
  geom_vline(xintercept=0.25, color='Blue', linewidth=1, linetype='dashed') +
  ylim(0.90,1.2) + xlim(0.12,0.3) +
  geom_text(hjust=-0.2, vjust=-0.5, show.legend = FALSE)  # Add this line to show labels

#### Part 2: timelapse of normalized leaf trait values ####

# Normalize trait values to the maximum value per trait #
combined <- combined %>%
  mutate(GS_normalized = GS / max(GS, na.rm = TRUE))
combined <- combined %>%
  mutate(ALEAF_normalized = ALEAF / max(ALEAF, na.rm = TRUE))
combined <- combined %>%
  mutate(Ci_normalized = Ci / max(Ci, na.rm = TRUE))
combined <- combined %>%
  mutate(Vcmax_normalized = Vcmax / max(Vcmax, na.rm = TRUE))

# Create the data frame with normalized values and pivot to correct format #
timelapse_dataframe <- combined %>%
  select(Step, GS_normalized, ALEAF_normalized, Ci_normalized, Vcmax_normalized) %>%
  pivot_longer(cols = c(GS_normalized, ALEAF_normalized, Ci_normalized, Vcmax_normalized),
               names_to = "Variable",
               values_to = "Normalized_value") %>%
  mutate(Variable = sub("_normalized", "", Variable))

# Create a plot for the timelapse grouping the values per trait value # 
timelapse_dataframe<-ggplot(timelapse_dataframe, aes(x=Step, y=Normalized_value, color=Variable)) +
  geom_line(linewidth = 1, aes(group=Variable)) +
  geom_point(size=3)+
  theme_bw()+
  ggtitle("Timelapse ")+
  theme_classic(base_size = 15)

print(timelapse_dataframe)
