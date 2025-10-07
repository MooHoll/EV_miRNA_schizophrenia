#####################################################
# Mary's miRNA data analysis: both batches of data
#####################################################

# ---------------------------------------------------
# Admin 
# ---------------------------------------------------

# Set directory, where the data file is
setwd("Z:/My Documents/PROJECTS/Mary_miRNAs")

# Load packages to do the stuff we need to
library(readr)
library(car)
library(ggplot2)
library(multcomp)
library(emmeans)
library(lme4)
library(lmerTest)
library(data.table)
library(dplyr)
library(reshape2)
library(doBy)
library(tidyr)
library(fitdistrplus)
library(ggbreak) 
library(patchwork)

# Read in the raw data
mary_data <- read_csv("normalised_data.txt")

# Everything needs to be a factor (this is just an R communication thing)
mary_data$Condition <- as.factor(mary_data$Condition)
mary_data$Sex <- as.factor(mary_data$Sex)

# ---------------------------------------------------
# Clean up data
# ---------------------------------------------------

# Remove rows with missing data
mary_data <- mary_data[complete.cases(mary_data),]

# Remove the samples with low isolation efficiency
mary_data <- mary_data[!(mary_data$Sample_ID=="1627"
                         | mary_data$Sample_ID=="S1030"),]

# Remove miRNAs which don't have at least 3 replicates for each condition/sex
mary_data$both <- paste0(mary_data$Condition, mary_data$Sex)
length(unique(mary_data$miRNA)) # 78

sample_number <- summarise(group_by(mary_data,miRNA,both),count =n())
sample_number <- spread(sample_number, key = both, value = count)
sample_number <- sample_number[(sample_number$HCF >=3 & sample_number$HCM >=3 &
                                  sample_number$SZF >=3 & sample_number$SZM >=3 ),]

mary_data <- mary_data[mary_data$miRNA %in% sample_number$miRNA,]
length(unique(mary_data$miRNA)) # 63

# ---------------------------------------------------
# Test run with just one miRNA
# ---------------------------------------------------

# Pull out data for just one as a test run
test <- mary_data[mary_data$miRNA=="hsa-let-7a-5p",]

# Take a look at distribution, it's not looking super 'normal'
hist(test$normalised)
qqPlot(test$normalised)

# Run linear model with interaction (basically an anova)
m = lm(normalised ~ Condition*Sex, data=test)
Anova(m, type=3, test.statistic="F")

# Test the normality formally - it's not normal
shapiro.test(residuals(m))

# What is the distribution? 
descdist(test$normalised, discrete = FALSE)

# Looks like exponential for some which is also analysed as Gamma, lets compare to normal, yup it's gamma
look <- fitdist(test$normalised, "gamma")
plot(look)
look <- fitdist(test$normalised, "norm")
plot(look)

# Run a general linear model instead (basically an anova which allows gamma)

# First need to account for marginal stuff
contrasts(test$Condition) <-"contr.sum"
contrasts(test$Sex) <-"contr.sum"

# Actual model, nothing is significant
m = glm(normalised ~ Condition*Sex, data=test, family=Gamma)
Anova(m, type=3)

# Quick eyeball of the data, yup checks out, doesn't look significant
ggplot(test, aes(x=Condition, y=normalised, fill=Sex))+
  geom_boxplot()

# Do a posthoc test anyway to check the code is working
summary(emmeans(m, pairwise ~ Condition*Sex, adjust="tukey", mode="linear.predictor", type="Score"))

# ---------------------------------------------------
# Everything above looks good, now lets look at the data distribution for all
# ---------------------------------------------------

# Have a look at the qqplots
ggplot(mary_data, aes(sample = normalised)) +
  stat_qq() +
  stat_qq_line() +
  facet_wrap( ~ miRNA, scales="free") # Mixed bag, some normal, some gamma

# Run the shapiro test for all
mary_data$miRNA <- as.factor(mary_data$miRNA)
normaility_check <- do.call(rbind,
                            lapply(levels(mary_data$miRNA),
                                   function (miRNA_name) { 
                                     dfsubs <- mary_data[mary_data$miRNA==miRNA_name,]
                                     
                                     fit <- lm(normalised ~ Condition*Sex, data=dfsubs)
                                     output <- shapiro.test(residuals(fit))
                                     
                                     df=data.frame(output[2], miRNA_name)
                                     df
                                   }
                            ) )

normaility_check$normal <- "Yes"
normaility_check$normal[normaility_check$p.value < 0.05] <- "No"
table(normaility_check$normal) # 4 normal, 59 gamma

# Put miRNA column back to character as the factor will mess up later
mary_data$miRNA <- as.character(mary_data$miRNA)

# Make two dataframes, one for the normal, one for the gamma
normal_miRNAs <- mary_data[mary_data$miRNA %in% 
                                   normaility_check$miRNA_name[normaility_check$normal =="Yes"],]
length(unique(normal_miRNAs$miRNA)) # 4

ggplot(normal_miRNAs, aes(sample = normalised)) +
  stat_qq() +
  stat_qq_line() +
  facet_wrap( ~ miRNA, scales="free")

gamma_miRNAs <- mary_data[mary_data$miRNA %in% 
                                  normaility_check$miRNA_name[normaility_check$normal =="No"],]
length(unique(gamma_miRNAs$miRNA)) # 59

ggplot(gamma_miRNAs, aes(sample = normalised)) +
  stat_qq() +
  stat_qq_line() +
  facet_wrap( ~ miRNA, scales="free")

# ---------------------------------------------------
# Models for normal and gamma data sets
# ---------------------------------------------------

# Run linear model (basically an anova) for all normally distributed miRNAs
normal_miRNAs$miRNA <- as.factor(normal_miRNAs$miRNA)
contrasts(normal_miRNAs$Condition) <-"contr.sum"
contrasts(normal_miRNAs$Sex) <-"contr.sum"

model_output_normal <- do.call(rbind,
                               lapply(levels(normal_miRNAs$miRNA),
                                      function (miRNA_name) { 
                                        dfsubs <- normal_miRNAs[normal_miRNAs$miRNA==miRNA_name,]
                                        
                                        fit <- lm(normalised ~ Condition*Sex, data=dfsubs)
                                        
                                        model_output <- data.frame(anova(fit))
                                        
                                        df <- data.frame(miRNA_name,model_output)
                                        colnames(df)=c("miRNA","Df","Sum.Sq","Mean.Sq","F.value","p")
                                        df
                                      }
                               ) )

# Tidy up
model_output_normal <- model_output_normal[complete.cases(model_output_normal),]
model_output_normal$comparison <- row.names(model_output_normal)
row.names(model_output_normal) <- NULL
model_output_normal <- model_output_normal[,c(1,5,6,7)]
model_output_normal$comparison <- gsub('[[:digit:]]+', '', model_output_normal$comparison)


# Run general linear model (basically an anova) for all gamma distributed miRNAs
gamma_miRNAs$miRNA <- as.factor(gamma_miRNAs$miRNA)
contrasts(gamma_miRNAs$Condition) <-"contr.sum"
contrasts(gamma_miRNAs$Sex) <-"contr.sum"

model_output_gamma <- do.call(rbind,
                              lapply(levels(gamma_miRNAs$miRNA),
                                     function (miRNA_name) { 
                                       dfsubs <- gamma_miRNAs[gamma_miRNAs$miRNA==miRNA_name,]
                                       
                                       fit = glm(normalised ~ Condition*Sex, data=dfsubs, family=Gamma)
                                       
                                       
                                       model_output <- data.frame(Anova(fit, type=3))
                                       
                                       df <- data.frame(miRNA_name,model_output)
                                       colnames(df)=c("miRNA","Chi.Sq","Df","p")
                                       df
                                     }
                              ) )

# Tidy up
model_output_gamma$comparison <- row.names(model_output_gamma)
row.names(model_output_gamma) <- NULL
model_output_gamma <- model_output_gamma[,-3]
model_output_gamma$comparison <- gsub('[[:digit:]]+', '', model_output_gamma$comparison)


# Put both dataframes together
all_model_output <- rbindlist(list(model_output_normal, model_output_gamma), fill=T)

# ---------------------------------------------------
# Posthoc test the significant ones for interaction effects
# ---------------------------------------------------

# Pull out significant miRNAs for interaction only
significant <- all_model_output[(all_model_output$p < 0.05 & 
                                   all_model_output$comparison=="Condition:Sex"),] 
length(unique(significant$miRNA)) # 7 miRNAs have something going on
write.table(significant, file="sig_interaction_miRNAs.txt", sep='\t', quote = F,
            col.names = T, row.names = F)


# Do manually so can eyeball each one

# Significant normal distribution:

#hsa-miR-143-3p (SZ-F diff to other 3)

# Change the name in the quotes below and run the chunk to get the posthoc p-values
miRNA <- "hsa-miR-143-3p"
data <- mary_data[mary_data$miRNA == miRNA,]
contrasts(data$Condition) <-"contr.sum"
contrasts(data$Sex) <-"contr.sum"
m = lm(normalised ~ Condition*Sex, data=data)
anova(m)
summary(emmeans(m, pairwise ~ Condition*Sex, adjust="tukey", mode="linear.predictor", type="Score"))


# Significant gamma distribution:

#hsa-miR-145-5p (interaction: but all lost in posthoc)
#hsa-miR-146a-5p (interaction: but all lost in posthoc)
#hsa-miR-17-5p (everything sig! posthoc: HC-F diff to other 3, found in first analysis)
#hsa-miR-19b-3p (interaction: but all lost in posthoc)
#hsa-miR-224-5p (interaction: but all lost in posthoc)
#hsa-miR-200c-3p (interaction: but all lost in posthoc)

# Change the name in the quotes below and run the chunk to get the posthoc p-values
miRNA <- "hsa-miR-17-5p"
data <- mary_data[mary_data$miRNA == miRNA,]
contrasts(data$Condition) <-"contr.sum"
contrasts(data$Sex) <-"contr.sum"
m = glm(normalised ~ Condition*Sex, data=data, family=Gamma)
Anova(m, type=3)
summary(emmeans(m, pairwise ~ Condition*Sex, adjust="tukey", mode="linear.predictor", type="Score"))

# ---------------------------------------------------
# Plot the significant interaction miRNAs
# ---------------------------------------------------

# Define function for confidence intervals
condifence_intervals <- function(x) {
  m <- mean(x)
  sd1 <- sd(x)
  n <- length(x)
  error <- qnorm(0.95)*sd1/sqrt(n)
  ymin <- m-error
  ymax <- m+error
  return(c(y=m,ymin=ymin,ymax=ymax))
}

# Make categories for adding the dots and mean
mary_data$combined <- paste0(mary_data$Condition,"_",mary_data$Sex)


# Pull out and plot interaction significant

#hsa-miR-143-3p (SZ-F diff to other 3)
#hsa-miR-17-5p (everything sig! posthoc: HC-F diff to other 3, found in first analysis)

int_sig <- as.data.frame(c("hsa-miR-143-3p","hsa-miR-17-5p"))
colnames(int_sig) <- "miRNA"
int_sig_raw <- mary_data[mary_data$miRNA %in%
                           int_sig$miRNA,]

ggplot(int_sig_raw, aes(y = normalised, x= combined, fill=Sex)) +
  geom_violin()+
  geom_point()+
  stat_summary(fun.data=condifence_intervals, color="red",size=0.5)+
  theme_bw()+
  scale_fill_manual(breaks=c("F","M"), values=c("seagreen2","dodgerblue2"))+
  facet_wrap( ~ miRNA, scales="free")+
  xlab("Condition")+
  ylab("Normalised Expression Value")+
  ggtitle("Significant interaction between condition and sex")+
  scale_x_discrete(limits=c("HC_F","HC_M","SZ_F","SZ_M"),
                   labels=c("HC","HC", "SZ","SZ"))


# ---------------------------------------------------
# Re-run the models without the interaction effect
# ---------------------------------------------------

# Remove ones which have interaction
length(unique(normal_miRNAs$miRNA)) #4
normal_miRNAs <- normal_miRNAs[!normal_miRNAs$miRNA=="hsa-miR-143-3p",]
length(unique(normal_miRNAs$miRNA)) #3
normal_miRNAs <- droplevels(normal_miRNAs)

# Run linear model (basically an anova) for all normally distributed miRNAs
contrasts(normal_miRNAs$Condition) <-"contr.sum"
contrasts(normal_miRNAs$Sex) <-"contr.sum"

model_output_normal <- do.call(rbind,
                               lapply(levels(normal_miRNAs$miRNA),
                                      function (miRNA_name) { 
                                        dfsubs <- normal_miRNAs[normal_miRNAs$miRNA==miRNA_name,]
                                        
                                        fit <- lm(normalised ~ Condition+Sex, data=dfsubs)
                                        
                                        model_output <- data.frame(anova(fit))
                                        
                                        df <- data.frame(miRNA_name,model_output)
                                        colnames(df)=c("miRNA","Df","Sum.Sq","Mean.Sq","F.value","p")
                                        df
                                      }
                               ) )

# Tidy up
model_output_normal <- model_output_normal[complete.cases(model_output_normal),]
model_output_normal$comparison <- row.names(model_output_normal)
row.names(model_output_normal) <- NULL
model_output_normal <- model_output_normal[,c(1,5,6,7)]
model_output_normal$comparison <- gsub('[[:digit:]]+', '', model_output_normal$comparison)


# Remove ones which have interaction
length(unique(gamma_miRNAs$miRNA)) #59
gamma_miRNAs <- gamma_miRNAs[!(gamma_miRNAs$miRNA=="hsa-miR-145-5p" |
                                 gamma_miRNAs$miRNA=="hsa-miR-146a-5p" |
                                 gamma_miRNAs$miRNA=="hsa-miR-17-5p" |
                                 gamma_miRNAs$miRNA=="hsa-miR-19b-3p" |
                                 gamma_miRNAs$miRNA=="hsa-miR-224-5p" |
                                 gamma_miRNAs$miRNA=="hsa-miR-200c-3p"),]
length(unique(gamma_miRNAs$miRNA)) #53
gamma_miRNAs <- droplevels(gamma_miRNAs)

# Run general linear model (basically an anova) for all gamma distributed miRNAs
contrasts(gamma_miRNAs$Condition) <-"contr.sum"
contrasts(gamma_miRNAs$Sex) <-"contr.sum"

model_output_gamma <- do.call(rbind,
                              lapply(levels(gamma_miRNAs$miRNA),
                                     function (miRNA_name) { 
                                       dfsubs <- gamma_miRNAs[gamma_miRNAs$miRNA==miRNA_name,]
                                       
                                       fit = glm(normalised ~ Condition+Sex, data=dfsubs, family=Gamma)
                                       
                                       
                                       model_output <- data.frame(Anova(fit, type=3))
                                       
                                       df <- data.frame(miRNA_name,model_output)
                                       colnames(df)=c("miRNA","Chi.Sq","Df","p")
                                       df
                                     }
                              ) )

# Tidy up
model_output_gamma$comparison <- row.names(model_output_gamma)
row.names(model_output_gamma) <- NULL
model_output_gamma <- model_output_gamma[,-3]
model_output_gamma$comparison <- gsub('[[:digit:]]+', '', model_output_gamma$comparison)


# Put both dataframes together
all_model_output <- rbindlist(list(model_output_normal, model_output_gamma), fill=T)

# ---------------------------------------------------
# Plot the significant condition miRNAs
# ---------------------------------------------------

# Pull out significant miRNAs for interaction only
significant <- all_model_output[(all_model_output$p < 0.05 & 
                                   all_model_output$comparison=="Condition"),] 
length(unique(significant$miRNA)) # 3 miRNAs have something going on
write.table(significant, file="sig_condition_miRNAs.txt", sep='\t', quote = F,
            col.names = T, row.names = F)


# Pull out and plot just condition significant

#hsa-miR-103a-3p...91
#hsa-miR-200b-3p 
#hsa-miR-30e-5p 

cond_sig <- as.data.frame(c("hsa-miR-103a-3p...91","hsa-miR-200b-3p",
                            "hsa-miR-30e-5p"))
colnames(cond_sig) <- "miRNA"
cond_sig_raw <- mary_data[mary_data$miRNA %in%
                            cond_sig$miRNA,]

ggplot(cond_sig_raw, aes(y = normalised, x= combined, fill=Sex)) +
  geom_violin()+
  geom_point()+
  stat_summary(fun.data=condifence_intervals, color="red",size=0.5)+
  theme_bw()+
  scale_fill_manual(breaks=c("F","M"), values=c("seagreen2","dodgerblue2"))+
  facet_wrap( ~ miRNA, scales="free")+
  xlab("Condition")+
  ylab("Normalised Expression Value")+
  ggtitle("Significant difference between conditions")+
  scale_x_discrete(limits=c("HC_F","HC_M","SZ_F","SZ_M"),
                   labels=c("HC","HC", "SZ","SZ"))


# ---------------------------------------------------
# Plot the significant sex miRNAs
# ---------------------------------------------------

# Pull out significant miRNAs for interaction only
significant <- all_model_output[(all_model_output$p < 0.05 & 
                                   all_model_output$comparison=="Sex"),] 
length(unique(significant$miRNA)) # 1 miRNAs have something going on
write.table(significant, file="sig_sex_miRNAs.txt", sep='\t', quote = F,
            col.names = T, row.names = F)

# Pull out and plot sex significant

#hsa-miR-423-5p

sex_sig <- as.data.frame(c("hsa-miR-423-5p"))
colnames(sex_sig) <- "miRNA"
sex_sig_raw <- mary_data[mary_data$miRNA %in%
                            sex_sig$miRNA,]

ggplot(sex_sig_raw, aes(y = normalised, x= combined, fill=Sex)) +
  geom_violin()+
  geom_point()+
  stat_summary(fun.data=condifence_intervals, color="red",size=0.5)+
  theme_bw()+
  scale_fill_manual(breaks=c("F","M"), values=c("seagreen2","dodgerblue2"))+
  facet_wrap( ~ miRNA, scales="free")+
  xlab("Condition")+
  ylab("Normalised Expression Value")+
  ggtitle("Significant difference between sexes")+
  scale_x_discrete(limits=c("HC_F","HC_M","SZ_F","SZ_M"),
                   labels=c("HC","HC", "SZ","SZ"))

