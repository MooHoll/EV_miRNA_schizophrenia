#####################################################
# Mary's miRNA data analysis: normalise by batch
#####################################################

# ---------------------------------------------------
# Admin 
# ---------------------------------------------------

# Set directory, where the data file is
setwd("Z:/My Documents/PROJECTS/Mary_miRNAs")

# Load packages to do the stuff we need to
library(readxl)
library(goeveg)
library(ggplot2)
library(data.table)
library(dplyr)
library(reshape2)
library(doBy)
library(tidyr)
library(car)
library(fitdistrplus)

# ---------------------------------------------------
# Check variation in data by groups (miR26a_126 normalised)
# ---------------------------------------------------

# Read in the raw data
mary_data <- read_excel("data_norm_by_miR26a_and_126.xlsx")

# Change dataframe from wide to long (easier format to work with)
long_data <- melt(setDT(mary_data), id.vars = c("Condition","Sex", "Sample_ID", "Batch", "Norm_by"), variable.name = "miRNA")
long_data$miRNA <- as.character(long_data$miRNA)

# Calculate the coefficient of variation within the data per miRNA per group
cv_data <- long_data %>% 
  group_by(Condition, Sex, Batch, miRNA) %>%
  summarise(
    "coeff_var" = cv(value, na.rm = T))

# Remove rows with NA and make a new column to help plotting
cv_data <- cv_data[complete.cases(cv_data),]
cv_data$group <- paste0(cv_data$Condition, cv_data$Sex)

# Boxplot of coefficient of variations
ggplot(cv_data, aes(x=group, y=coeff_var, fill=Batch))+
  geom_boxplot()

# Have a look at distribution for stats
hist(cv_data$coeff_var)
qqPlot(cv_data$coeff_var)

m = lm(coeff_var ~ Condition*Sex*Batch, data=cv_data)
Anova(m, type=3, test.statistic="F")

shapiro.test(residuals(m)) # Not normally distributed

descdist(cv_data$coeff_var, discrete = FALSE)
look <- fitdist(cv_data$coeff_var, "gamma")
plot(look) # Yes, is gamma

# Run stat to see if the variation is indeef different between groups
m = glm(coeff_var ~ Condition*Sex*Batch, data=cv_data, family=Gamma)
Anova(m, type=3)

# No interactions so run without the interaction effect
m = glm(coeff_var ~ Condition+Sex+Batch, data=cv_data, family=Gamma)
Anova(m, type=3)

# Yes, there is a significant difference in the variation of the miRNA data
# based on the sex and batch


# ---------------------------------------------------
# Check variation in data by groups (miR16_let7a normalised)
# ---------------------------------------------------

# Read in the raw data
mary_data <- read_excel("data_norm_by_miR16_and_let7a.xlsx")

# Change dataframe from wide to long (easier format to work with)
long_data <- melt(setDT(mary_data), id.vars = c("Condition","Sex", "Sample_ID", "Batch", "Norm_by"), variable.name = "miRNA")
long_data$miRNA <- as.character(long_data$miRNA)


# Calculate the coefficient of variation within the data per miRNA per group
cv_data <- long_data %>% 
  group_by(Condition, Sex, Batch, miRNA) %>%
  summarise(
    "coeff_var" = cv(value, na.rm = T))

# Remove rows with NA and make a new column to help plotting
cv_data <- cv_data[complete.cases(cv_data),]
cv_data$group <- paste0(cv_data$Condition, cv_data$Sex)

# Boxplot of coefficient of variations
ggplot(cv_data, aes(x=group, y=coeff_var, fill=Batch))+
  geom_boxplot()

# Have a look at distribution for stats
hist(cv_data$coeff_var)
qqPlot(cv_data$coeff_var)

m = lm(coeff_var ~ Condition*Sex*Batch, data=cv_data)
Anova(m, type=3, test.statistic="F")

shapiro.test(residuals(m)) # Not normally distributed

descdist(cv_data$coeff_var, discrete = FALSE)

look <- fitdist(cv_data$coeff_var, "gamma")
plot(look) # Yes, is gamma

# Run stat to see if the variation is indeef different between groups
m = glm(coeff_var ~ Condition*Sex*Batch, data=cv_data, family=Gamma)
Anova(m, type=3)

# No interactions so run without the interaction effect
m = glm(coeff_var ~ Condition+Sex+Batch, data=cv_data, family=Gamma)
Anova(m, type=3)

# Significant effect of batch

# Decision point: we go with the second normalised dataset and I normalise by batch
# following Flav's logic

# ---------------------------------------------------
# Normalise by batch - using HC Females
# ---------------------------------------------------

# Read in the raw data
mary_data <- read_excel("data_norm_by_miR16_and_let7a.xlsx")

# Change dataframe from wide to long (easier format to work with)
long_data <- melt(setDT(mary_data), id.vars = c("Condition","Sex", "Sample_ID", "Batch", "Norm_by"), variable.name = "miRNA")
long_data$miRNA <- as.character(long_data$miRNA)

# ---------------------
# Batch 1
# ---------------------

batch_1 <- long_data[long_data$Batch=="Batch_1",]

# Get the control group for batch 1
control_data <- batch_1[(batch_1$Sex=="F" & batch_1$Condition=="HC"),]

# Get mean value per miRNA
mean_per_miRNA <- control_data %>% 
  group_by(miRNA) %>%
  summarise(
    "mean" = mean(value, na.rm = T))

# Remove NAs
mean_per_miRNA <- mean_per_miRNA[complete.cases(mean_per_miRNA),]

# Normalise by mean
all_data <- merge(batch_1, mean_per_miRNA, by = "miRNA", all=T)
all_data$normalised <- as.numeric(all_data$value / all_data$mean)

# Eyeball
mean(all_data$normalised[(all_data$Sex=="F" & all_data$Condition=="HC")], na.rm=T)
ggplot(all_data, aes(x = Condition, y = normalised, fill = Sex))+
  geom_boxplot()

# ---------------------
# Batch 2
# ---------------------

batch_2 <- long_data[long_data$Batch=="Batch_2",]

# Get the control group for batch 1
control_data <- batch_2[(batch_2$Sex=="F" & batch_2$Condition=="HC"),]

# Get mean value per miRNA
mean_per_miRNA <- control_data %>% 
  group_by(miRNA) %>%
  summarise(
    "mean" = mean(value, na.rm = T))

# Remove NAs
mean_per_miRNA <- mean_per_miRNA[complete.cases(mean_per_miRNA),]

# Normalise by mean
all_data_2 <- merge(batch_2, mean_per_miRNA, by = "miRNA", all=T)
all_data_2$normalised <- as.numeric(all_data_2$value / all_data_2$mean)

# Eyeball
mean(all_data_2$normalised[(all_data_2$Sex=="F" & all_data_2$Condition=="HC")], na.rm=T)
ggplot(all_data_2, aes(x = Condition, y = normalised, fill = Sex))+
  geom_boxplot()

# ---------------------
# Put batches back together into one dataframe
# ---------------------

both <- rbind(all_data, all_data_2) # same number of rows as original data
write.csv(both, file = "normalised_data.txt", row.names = F, quote = F)

# Double check the batch effect has gone as above

# Calculate the coefficient of variation within the data per miRNA per group
cv_data <- both %>% 
  group_by(Condition, Sex, Batch, miRNA) %>%
  summarise(
    "coeff_var" = cv(normalised, na.rm = T))

# Remove rows with NA and make a new column to help plotting
cv_data <- cv_data[complete.cases(cv_data),]
cv_data$group <- paste0(cv_data$Condition, cv_data$Sex)

# Boxplot of coefficient of variations
ggplot(cv_data, aes(x=group, y=coeff_var, fill=Batch))+
  geom_boxplot()

# Have a look at distribution for stats
hist(cv_data$coeff_var)
qqPlot(cv_data$coeff_var)

m = lm(coeff_var ~ Condition*Sex*Batch, data=cv_data)
Anova(m, type=3, test.statistic="F")

shapiro.test(residuals(m)) # Not normally distributed

descdist(cv_data$coeff_var, discrete = FALSE)

look <- fitdist(cv_data$coeff_var, "gamma")
plot(look) # Yes, is gamma

# Run stat to see if the variation is indeef different between groups
m = glm(coeff_var ~ Condition*Sex*Batch, data=cv_data, family=Gamma)
Anova(m, type=3)

# No interactions so run without the interaction effect
m = glm(coeff_var ~ Condition+Sex+Batch, data=cv_data, family=Gamma)
Anova(m, type=3)


# ---------------------------------------------------
# Normalise by batch - using HC Males
# ---------------------------------------------------

# Read in the raw data
mary_data <- read_excel("data_norm_by_miR16_and_let7a.xlsx")

# Change dataframe from wide to long (easier format to work with)
long_data <- melt(setDT(mary_data), id.vars = c("Condition","Sex", "Sample_ID", "Batch", "Norm_by"), variable.name = "miRNA")
long_data$miRNA <- as.character(long_data$miRNA)

# ---------------------
# Batch 1
# ---------------------

batch_1 <- long_data[long_data$Batch=="Batch_1",]

# Get the control group for batch 1
control_data <- batch_1[(batch_1$Sex=="M" & batch_1$Condition=="HC"),]

# Get mean value per miRNA
mean_per_miRNA <- control_data %>% 
  group_by(miRNA) %>%
  summarise(
    "mean" = mean(value, na.rm = T))

# Remove NAs
mean_per_miRNA <- mean_per_miRNA[complete.cases(mean_per_miRNA),]

# Normalise by mean
all_data <- merge(batch_1, mean_per_miRNA, by = "miRNA", all=T)
all_data$normalised <- as.numeric(all_data$value / all_data$mean)

# Eyeball
mean(all_data$normalised[(all_data$Sex=="M" & all_data$Condition=="HC")], na.rm=T)
ggplot(all_data, aes(x = Condition, y = normalised, fill = Sex))+
  geom_boxplot()

# ---------------------
# Batch 2
# ---------------------

batch_2 <- long_data[long_data$Batch=="Batch_2",]

# Get the control group for batch 1
control_data <- batch_2[(batch_2$Sex=="M" & batch_2$Condition=="HC"),]

# Get mean value per miRNA
mean_per_miRNA <- control_data %>% 
  group_by(miRNA) %>%
  summarise(
    "mean" = mean(value, na.rm = T))

# Remove NAs
mean_per_miRNA <- mean_per_miRNA[complete.cases(mean_per_miRNA),]

# Normalise by mean
all_data_2 <- merge(batch_2, mean_per_miRNA, by = "miRNA", all=T)
all_data_2$normalised <- as.numeric(all_data_2$value / all_data_2$mean)

# Eyeball
mean(all_data_2$normalised[(all_data_2$Sex=="M" & all_data_2$Condition=="HC")], na.rm=T)
ggplot(all_data_2, aes(x = Condition, y = normalised, fill = Sex))+
  geom_boxplot()

# ---------------------
# Put batches back together into one dataframe
# ---------------------

both <- rbind(all_data, all_data_2) # same number of rows as original data
#write.csv(both, file = "normalised_data.txt", row.names = F, quote = F)

# Double check the batch effect has gone as above

# Calculate the coefficient of variation within the data per miRNA per group
cv_data <- both %>% 
  group_by(Condition, Sex, Batch, miRNA) %>%
  summarise(
    "coeff_var" = cv(normalised, na.rm = T))

# Remove rows with NA and make a new column to help plotting
cv_data <- cv_data[complete.cases(cv_data),]
cv_data$group <- paste0(cv_data$Condition, cv_data$Sex)

# Boxplot of coefficient of variations
ggplot(cv_data, aes(x=group, y=coeff_var, fill=Batch))+
  geom_boxplot()

# Have a look at distribution for stats
hist(cv_data$coeff_var)
qqPlot(cv_data$coeff_var)

m = lm(coeff_var ~ Condition*Sex*Batch, data=cv_data)
Anova(m, type=3, test.statistic="F")

shapiro.test(residuals(m)) # Not normally distributed

descdist(cv_data$coeff_var, discrete = FALSE)

look <- fitdist(cv_data$coeff_var, "gamma")
plot(look) # Yes, is gamma

# Run stat to see if the variation is indeef different between groups
m = glm(coeff_var ~ Condition*Sex*Batch, data=cv_data, family=Gamma)
Anova(m, type=3)

# No interactions so run without the interaction effect
m = glm(coeff_var ~ Condition+Sex+Batch, data=cv_data, family=Gamma)
Anova(m, type=3)

