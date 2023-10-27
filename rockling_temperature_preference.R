# Bass Dye, created: 30/05/2023, last revision: 05/06/2023
# Code to 1) read in rockling experimental  data from csv file 
# 2) Statistically analyze the data:
  # test for non random use of chamber, 
  # ANOVA - test for differences in compartment usage observational phases,
  # rank (most preferred to least preferred) temperature compartments
  # statistical test rankings.
# Code based off of Aebischer et al. 1993 and Schram et al. 2013 to analyse data from the preference chamber experiments.
# Fish acclimated to 16°C were exposed to a temperature range of 10-13-16-19-22°C. 
  # dependent variables (measured): count of individual fish in each preference chamber compartment 
  # independent variables (observation phases): non-gradient, acute temperature gradient, final temperature gradient

#	Step 1: read in the raw data: counts of numbers of fish per zone in the preference chamber, with column headings:
#	experiment_group: Indicator for experiment (experimental group): 1,2,...,22
#	acclimation_temperature: Acclimation temperature of an experiment: the temperature of all 
#				zones within the preference chamber during the acclimation phase
# compartment: Indicator for compartment within the preference chamber: 1,2,...,8
# temperature: Temperature of a compartment in the acute and final temperature observation phases
# phase: Phase of the experiment indicated by roman numerals I (non-gradient),II (acute temp pref.) or III (final temp pref.)
#	      T0, T5,...,T60: Count location of fish at the start of the phase (T0), at 5 minutes from the start of the phase (T5),….,at 60 minutes from the start of 
#				the phase (T60).

# House keeping; clean environment 
rm(list=ls())

# Install and load packages ####
if(!require(fs)) install.packages("fs")
if(!require(tidyverse))install.packages("tidyverse")
if(!require(here))install.packages("here")
if(!require(ggpubr))install.packages("ggpubr")

# Load packages
library("fs")
library("tidyverse")
library("here")
library("ggpubr")

# Read in data
dat <- read.table('rockling_temp_preference_data.csv', header=TRUE, sep=",");
#View(dat)

# Step 2: transform data of raw counts per compartment to proportions per temperature zone 
# 2a: sum of raw counts over all counting times
dat$Sum <- rowSums(dat[, c("T0","T5","T10","T15","T20","T25","T30","T35","T40","T45","T50","T55","T60")]); 

# 2b: create a dataset ‘aggdat’ with aggregate counts per experiment, acclimation temperature, phase, and compartment
aggdat <- aggregate(dat$Sum, list(as.factor(dat$experiment),
                                  as.factor(dat$acclimation_temperature),
                                  as.factor(dat$phase),
                                  as.factor(dat$compartment)),
                    sum);

# 2c: create a dataset ‘aggdatsum’ with aggregate counts per experiment, acclimation, temperature and phase 
aggdatsum <- aggregate(dat$Sum,
                       list(as.factor(dat$experiment),
                            as.factor(dat$acclimation_temperature),
                            as.factor(dat$phase)),
                       sum); 

# Create a dataset ‘expdat’ with each row representing a combination of experiment and phase
expdat <- data.frame(aggdatsum,
                     z1=rep(NA, nrow(aggdatsum)),
                     z2=rep(NA, nrow(aggdatsum)),
                     z3=rep(NA, nrow(aggdatsum)), 
                     z4=rep(NA, nrow(aggdatsum)),
                     z5=rep(NA, nrow(aggdatsum)),
                     z6=rep(NA, nrow(aggdatsum)),
                     z7=rep(NA, nrow(aggdatsum)), 
                     z8=rep(NA,nrow(aggdatsum)));
for ( ii in 1:nrow(aggdatsum) ) {
  ttt <- aggdat[aggdat$Group.1==aggdatsum$Group.1[ii] & aggdat$Group.2==aggdatsum$Group.2[ii] & aggdat$Group.3==aggdatsum$Group.3[ii], ];
  expdat[ii, 5:12] <- ttt$x/expdat$x[ii];
};

# Change column names
colnames(expdat)[c(1:4)] <- c("experiment", "acclimation_temp", "phase", "count")

# 2d: create data set with proportions per temperature zone, by summing up the proportions of
# usage of compartments with the same temperature zone
# C1: 22 degrees: compartment 4
# C2: 19 degrees: compartments 3 and 5
# C3: 16 degrees: compartments 2 and 6
# C4: 13 degrees: compartments 1 and 7
# C5: 10 degrees: compartment 8
prop_dat <- expdat[, 1:4];
prop_dat$C1 <- expdat$z4;
prop_dat$C2 <- expdat$z3 + expdat$z5;
prop_dat$C3 <- expdat$z2 + expdat$z6;
prop_dat$C4 <- expdat$z1 + expdat$z7;
prop_dat$C5 <- expdat$z8;

# Zero inflated data
#Plot number of zeros in data to show zero inflation
#library(dplyr)
#temp2 <- expdat[,5:12]
#temp2 <- data.frame(value=unlist(temp2, use.names = FALSE))
#count_zero <- temp2[temp2 == 0] 
#percent_zero <- length(count_zero) / nrow(temp2)
#h <- hist(temp2$value, main = "Zero inflated data: 56% = 0 ", xlab = "Percent compartment use",
#breaks = c(0,seq(0.01,1.0,.01)))
#h

# Find and replace zero numerator or denominator with small value
  # As a zero numerator or denominator in the log-ratio transformation is invalid, a small positive value, 
  # less than the smallest recorded nonzero proportion, should be substituted (Aebischer et al. 1993)
  # replacing zeros in original data set with values less than the smallest recorded nonzero proportion
  # "An order of magnitude is generally sufficient" (Aebischer et al. 1993)
  # Find smallest value and replace with much smaller value in prop_data 
temp <- prop_dat[, 5:9]
replace_value <- min(temp[temp > 0])  # smallest value in our data set is 0.02564103
# Replace 0 with small value
prop_dat[prop_dat == 0] <- replace_value / 1000; # 1000

# Step 3: compute log-ratios ####
# An appropriate statistical framework for handling compositional data is to replace the observed proportions 
# With a set of ratios by choosing the observed fraction of one particular component as the denominator by which all other fractions are divided
# The results of the analyses do not depend upon the choice of denominator.
# In this way, the unit sum constraint on the composition is broken (Aitchison (1986) The statistical analysis of compositional data. London: Chapman and Hall. 416 p.)

# A version of the central limit theorem exists stating why random variation in log-ratios can often be assumed to be normally distributed (Aitchison, 1986)
# Take the compartment C3 (21*C temperature compartment) as the denominator
prop_dat$LR1 <- log(prop_dat$C1/prop_dat$C3);
prop_dat$LR2 <- log(prop_dat$C2/prop_dat$C3);
prop_dat$LR3 <- log(prop_dat$C4/prop_dat$C3);
prop_dat$LR4 <- log(prop_dat$C5/prop_dat$C3);

# Step 4: offset log-ratios of usage with log-ratios of availability of temperature zones ####
# I.e. when the temperature gradient is present within the chamber, 1/8 of the chamber is 15ºC while 1/4 is 18ºC.
# E.g. prop_dat$LRA1 <- prop_dat$LR1 - log(0.5) is d = yu - ya (pairwise habitat differences = utilized habitat - available habitat; Aebischer et al. 1993)
prop_dat$LRA1 <- prop_dat$LR1 - log(0.5); # 0.5 = (1/8)/(1/4)
prop_dat$LRA2 <- prop_dat$LR2 - log(1); # 1 = (1/4)/(1/4)
prop_dat$LRA3 <- prop_dat$LR3 - log(1);
prop_dat$LRA4 <- prop_dat$LR4 - log(0.5);

# Step 5: assess evidence of non-random usage (i.e. preference) during the three observational phases: ####
# 5a: create a matrix ‘CM’ of log-ratios of the non-gradient phase (I) or temperature gradient phases (II and III); substitute phase below accordingly
BAdat <- subset(prop_dat, phase == "I") 
BALR1_non <- BAdat$LRA1
BALR2_non <- BAdat$LRA2 
BALR3_non <- BAdat$LRA3
BALR4_non <- BAdat$LRA4

CM_non <- rbind(BALR1_non, BALR2_non, BALR3_non, BALR4_non);

# 5b: compute the generalized likelihood ratio statistic which compares a model with no preference
# against a model with preference (Aebischer, 1993).
# H1MAT - reduced model (m2); matrix of raw sums of squares and cross products calculated from the matrix ‘CM’ of log-ratios
# H2MAT - general model (m1); matrix of mean-corrected sums of squares and cross-products calculated from the matrix ‘CM’ of log-ratios
H1MAT <- H2MAT <- matrix(data= NA, nrow = 4, ncol = 4); 
for (i in 1:4 ) {
  for (j in 1:4) {
    H1MAT[i, j] <- sum(CM_non[i, ]  * CM_non[j, ]);
    H2MAT[i, j] <- sum((CM_non[i, ] - mean(CM_non[i, ])) * (CM_non[j, ] - mean(CM_non[j, ])));
  };
};

# Test for general preference (deviation from random use)
# -N*ln($) $=|r1|/|r2| ~ ln(det(H2MAT)/det(H1MAT)) ~ ln(general model/reduced model) - Aebischer et al. 1993
TEST <- -ncol(CM_non) * log(det(H2MAT) / det(H1MAT)); TEST;

# By default, pchisq() gives the proportion of the distribution to the left of the value.
# To get the proportion more extreme than your difference, you can specify lower.tail = FALSE or subtract the result from 1.
# Used to compute cumulative chi square density for a vector of elements
# If significant (p<0.05), evidence for non-random use (null - no preference;  alternative - preference)
1 - pchisq(TEST, df = 4);

rm(list=setdiff(ls(), c("prop_dat"))) # remove variables

# 5c: The matrix of log-ratios may also be analysed using MANOVA, for example to assess
# evidence of treatment effects:
# MANOVA tests for temperature preference experiments

# Test for difference between phases ####
# Select which experimental phase you would like to test (non-gradient - (I), acute - (II), final - (III))
data_phase <- as_tibble(prop_dat) %>% 
  filter(phase != "I") 

# Conduct MANOVA although multivariate normality assumption is not met
#summary(manova(cbind(LRA1, LRA2, LRA3, LRA4) ~ treatment +  Error(factor(subject) / treatment), data = data_phase))
summary(manova(cbind(LRA1, LRA2, LRA3, LRA4) ~ phase, data = data_phase))

# Conduct Permutation test so we can ignore MANOVA assumptions 
# Produce the same values by setting the seed
set.seed(1234)
# Define number of permutations
nperm <- 5000
# Allocate storage (i.e. vector with NAs to be replaced by permutation results)
# Here we want 5000 + 1 for our observed stat
res <- rep(NA, nperm + 1)

for(i in seq(along = res[-1])) {
  # First create matrix containing log-ratios of usage (corrected for the availability of temperature zones)
  # from the selected phase on line 10
  data_to_shuffle <- cbind(data_phase$LRA1, data_phase$LRA2, data_phase$LRA3, data_phase$LRA4) 
  # Then shuffle (randomly rearrange) the matrix for use in the MANOVA
  shuffled <- data_to_shuffle[sample(1:nrow(data_to_shuffle)), ]
  # Here we replace the ith value in the vector res with the shuffled data MANOVA result
  tmp <- summary(manova(shuffled ~ data_phase$phase), test = "Wilks")
  # Get Wilks value (test statistic) from MANOVA summary
  res[i] <- tmp[["stats"]][1, 2]
  # Remove variables for next iteration
  # rm("tmp", "shuffled")
}
# Now we append the observed stat onto the end of the result vector res
# We also store this in 'obs' for convenience
tmp <- summary(manova(cbind(LRA1, LRA2, LRA3, LRA4) ~ phase, data = data_phase),
               test = "Wilks")
res[nperm+1] <- obs <- tmp[["stats"]][1, 2]

# This is the permutation p-value - the proportion of the nperm
# permutations + 1 that are greater than or equal to the observed stat 'obs'
# perm_p_value <- sum(res <= obs) / (nperm+1); perm_p_value
# obs_p_value <- tmp[["Error: factor(subject):treatment"]][["stats"]][1,6]; obs_p_value
perm_p_value <- 2 * (min(sum(res <= obs), sum(res >= obs)) / (nperm+1)); perm_p_value # 2 tailed
obs_p_value <- tmp[["stats"]][1, 6]; obs_p_value

# Plot histogram of permutation and observed (blue line) results 
hist(res)
abline(v = obs, col = "blue")

# Clean up environment
rm(list=setdiff(ls(), "prop_dat"))

# Step 6: Once overall evidence of non-random use has been assessed, the next step is to rank compartments from least to most preferred ####
# This can be done by computing a cross-table with pairwise differences between matching log-ratios of usage 
# and availability of temperature zones, and counting the number of times a particular temperature zone
# has been observed to be preferred over other temperature zones (e.g. see table 1 in Aebischer et al. 1993)

# Subset for non-gradient, acute, or final proportions; choose (I,II,III) accordingly 
exp_data <- subset(prop_dat, phase == "II")

# Proportion of available space per compartment and temperature 
c1 <- 1/8 # C1: 33 degrees: compartment 4
c2 <- 2/8 # C2: 30 degrees: compartments 3 and 5
c3 <- 2/8 # C3: 27 degrees: compartments 2 and 6
c4 <- 2/8 # C4: 24 degrees: compartments 1 and 7
c5 <- 1/8 # C5: 21 degrees: compartment 8

# Create empty list to store matrices
all_mat <- list()
for (tt in 1:nrow(exp_data)) {
  all_mat[[tt]] <- matrix(data = NA, nrow = nrow(exp_data), ncol = 5); 
}

# Fill in matrix for each experimental group (n = 5); (table 1 in Aebischer et al. 1993)
for (i in 1:nrow(exp_data)) {
  # Create blank matrix to temporarily hold data for each iteration
  rank_mat <- matrix(data = NA, nrow = 5, ncol = 5); 
  # Fill in matrix 
  # Column 1
  rank_mat[2,1] <- ifelse(log(exp_data$C2[i]/exp_data$C1[i]) == 0, 0, log(exp_data$C3[i]/exp_data$C1[i]) - log(c3/c1))
  rank_mat[3,1] <- ifelse(log(exp_data$C3[i]/exp_data$C1[i]) == 0, 0, log(exp_data$C3[i]/exp_data$C1[i]) - log(c3/c1))
  rank_mat[4,1] <- ifelse(log(exp_data$C4[i]/exp_data$C1[i]) == 0, 0, log(exp_data$C4[i]/exp_data$C1[i]) - log(c4/c1))
  rank_mat[5,1] <- ifelse(log(exp_data$C5[i]/exp_data$C1[i]) == 0, 0, log(exp_data$C5[i]/exp_data$C1[i]) - log(c5/c1))
  
  # Column 2
  rank_mat[1,2] <- ifelse(log(exp_data$C1[i]/exp_data$C2[i]) == 0, 0, log(exp_data$C1[i]/exp_data$C2[i]) - log(c1/c2))
  rank_mat[3,2] <- ifelse(log(exp_data$C3[i]/exp_data$C2[i]) == 0, 0, log(exp_data$C3[i]/exp_data$C2[i]) - log(c3/c2))
  rank_mat[4,2] <- ifelse(log(exp_data$C4[i]/exp_data$C2[i]) == 0, 0, log(exp_data$C4[i]/exp_data$C2[i]) - log(c4/c2))
  rank_mat[5,2] <- ifelse(log(exp_data$C5[i]/exp_data$C2[i]) == 0, 0, log(exp_data$C5[i]/exp_data$C2[i]) - log(c5/c2))
  
  # Column 3
  rank_mat[1,3] <- ifelse(log(exp_data$C1[i]/exp_data$C3[i]) == 0, 0, log(exp_data$C1[i]/exp_data$C3[i]) - log(c1/c3))
  rank_mat[2,3] <- ifelse(log(exp_data$C2[i]/exp_data$C3[i]) == 0, 0, log(exp_data$C2[i]/exp_data$C3[i]) - log(c2/c3))
  rank_mat[4,3] <- ifelse(log(exp_data$C4[i]/exp_data$C3[i]) == 0, 0, log(exp_data$C4[i]/exp_data$C3[i]) - log(c4/c3))
  rank_mat[5,3] <- ifelse(log(exp_data$C5[i]/exp_data$C3[i]) == 0, 0, log(exp_data$C5[i]/exp_data$C3[i]) - log(c5/c3))
  
  # Column 4
  rank_mat[1,4] <- ifelse(log(exp_data$C1[i]/exp_data$C4[i]) == 0, 0, log(exp_data$C1[i]/exp_data$C4[i]) - log(c1/c4))
  rank_mat[2,4] <- ifelse(log(exp_data$C2[i]/exp_data$C4[i]) == 0, 0, log(exp_data$C2[i]/exp_data$C4[i]) - log(c2/c4))
  rank_mat[3,4] <- ifelse(log(exp_data$C3[i]/exp_data$C4[i]) == 0, 0, log(exp_data$C3[i]/exp_data$C4[i]) - log(c3/c4))
  rank_mat[5,4] <- ifelse(log(exp_data$C5[i]/exp_data$C4[i]) == 0, 0, log(exp_data$C5[i]/exp_data$C4[i]) - log(c5/c4))
  
  # Column 5
  rank_mat[1,5] <- ifelse(log(exp_data$C1[i]/exp_data$C5[i]) == 0, 0, log(exp_data$C1[i]/exp_data$C5[i]) - log(c1/c5))
  rank_mat[2,5] <- ifelse(log(exp_data$C2[i]/exp_data$C5[i]) == 0, 0, log(exp_data$C2[i]/exp_data$C5[i]) - log(c2/c5))
  rank_mat[3,5] <- ifelse(log(exp_data$C3[i]/exp_data$C5[i]) == 0, 0, log(exp_data$C3[i]/exp_data$C5[i]) - log(c3/c5))
  rank_mat[4,5] <- ifelse(log(exp_data$C4[i]/exp_data$C5[i]) == 0, 0, log(exp_data$C4[i]/exp_data$C5[i]) - log(c4/c5))
  
  all_mat[[i]] <- rank_mat
  rm(rank_mat) # Remove variable for next iteration
}

# Calculate element wise mean, sd, se
# Make a 3D array from list of matrices for calculation purposes
arr <- array(unlist(all_mat), c(5, 5, nrow(exp_data)))
# Calculate mean, sd, se of third dimension (i.e. mean of all experiments)
rank_mat_mean <- apply(arr, 1:2, mean)
rank_mat_sd <- apply(arr, 1:2, sd)
rank_mat_var <- apply(arr, 1:2, var)
rank_mat_se <-  rank_mat_sd / sqrt(nrow(arr))
rank_mat_test_stat <- abs(rank_mat_mean / rank_mat_se)

# Sum the number of positive values for each row to rank the habitats in increasing order of preference (i.e. greatest number = most preferred)
# Add empty column to matrix to store data
rank_mat_mean <- cbind(rank_mat_mean, NA)
rank_mat_mean[is.na(rank_mat_mean)] <- 0 # Change NAs to 0 (changing the diagonal NAs to zero doesn't change calculation)
for (ii in 1:nrow(rank_mat_mean)) {
  tt <- rank_mat_mean[ii, ] > 0
  rank_mat_mean[ii, 6] <- length(tt[tt] == TRUE)
}
# Temperature ranking matrix
row.names(rank_mat_mean) <- c("22*C", "19*C", "16*C", "13*C", "10*C")
colnames(rank_mat_mean) <- c("22*C", "19*C", "16*C", "13*C", "10*C", "Rankings")
rank_mat_mean

# Test for normality: Shapiro-Wilk's method
# p value > 0.05 implies distribution of data is not significantly different from normal distribution (i.e. normally distributed).
test_norm_mean <- c(rank_mat_mean[1, c(2:5)], rank_mat_mean[2, c(3:5)], rank_mat_mean[3, c(4:5)], rank_mat_mean[4, c(5)])
test_norm_se <- c(rank_mat_se[1, c(2:5)], rank_mat_se[2, c(3:5)], rank_mat_se[3, c(4:5)], rank_mat_se[4, c(5)])

shapiro.test(test_norm_mean)
shapiro.test(test_norm_se) 

# Test significance between compartments ####
# Subset of non-gradient, acute, or final proportions; choose (I,II,III) accordingly
#rank_data <- subset(prop_dat, phase == "I" & treatment == "5") 
rank_data <- exp_data # Same subset as line 255 to reduce any potential mistakes

# Clear space in environment
rm(list=setdiff(ls(), c("rank_data", "prop_dat", "rank_mat_mean")))

# Proportion of available space per compartment and temperature 
c1 <- 1/8 # C1: 33 degrees: compartment 4
c2 <- 2/8 # C2: 30 degrees: compartments 3 and 5
c3 <- 2/8 # C3: 27 degrees: compartments 2 and 6
c4 <- 2/8 # C4: 24 degrees: compartments 1 and 7
c5 <- 1/8 # C5: 21 degrees: compartment 8

# Temperature compartments comparison
compare_seq <- c("compare_22_19", "compare_22_16", "compare_22_13", "compare_22_10", "compare_19_16",
                 "compare_19_13", "compare_19_10", "compare_16_13", "compare_16_10", "compare_13_10")

# Create empty data frame to store matrices
OG_mat <- list()
# Create empty data frame to store p values
OG_p_vals <- data.frame(p_values = rep(NA, 1, 10))
# Rename rows corresponding to temperature compartment comparisons
row.names(OG_p_vals) <- compare_seq

# Run loop to calculate test statistic and resulting p value for each temperature compartment comparison and experimental phase
for (ttt in 1:length(compare_seq)) {
  # Counter
  test <- compare_seq[ttt] 
  for (xxx in 1:nrow(rank_data)) { 
    # Create empty data frame to store loop value
    OG_mat[[xxx]] <- matrix(data=NA, nrow=1, ncol=1); 
    for (i in 1:nrow(rank_data)) {
      # Create blank matrix to temporarily hold data for each iteration
      OG_rank_mat <- matrix(data=NA, nrow=1, ncol=1); 
      
      # Fill in matrix 
      if (test == "compare_22_19") {
        print("comparing 22*C to 19*C")
        OG_rank_mat[1,1] <- ifelse(log(rank_data$C1[i]/rank_data$C2[i]) == 0, 0, log(rank_data$C1[i]/rank_data$C2[i]) - log(c1/c2))
      } else if (test == "compare_22_16") {
        print("comparing 22*C to 16*C")
        OG_rank_mat[1,1] <- ifelse(log(rank_data$C1[i]/rank_data$C3[i]) == 0, 0, log(rank_data$C1[i]/rank_data$C3[i]) - log(c1/c3))
      } else if (test == "compare_22_13") {
        print("comparing 22*C to 13*C")
        OG_rank_mat[1,1] <- ifelse(log(rank_data$C1[i]/rank_data$C4[i]) == 0, 0, log(rank_data$C1[i]/rank_data$C4[i]) - log(c1/c4))
      } else if (test == "compare_22_10") {
        print("comparing 22*C to 10*C")
        OG_rank_mat[1,1] <- ifelse(log(rank_data$C1[i]/rank_data$C5[i]) == 0, 0, log(rank_data$C1[i]/rank_data$C5[i]) - log(c1/c5))
      } else if (test == "compare_19_16") {
        print("comparing 19*C to 16*C")
        OG_rank_mat[1,1] <- ifelse(log(rank_data$C2[i]/rank_data$C3[i]) == 0, 0, log(rank_data$C2[i]/rank_data$C3[i]) - log(c2/c3))
      } else if (test == "compare_19_13") {
        print("comparing 19*C to 13*C")
        OG_rank_mat[1,1] <- ifelse(log(rank_data$C2[i]/rank_data$C4[i]) == 0, 0, log(rank_data$C2[i]/rank_data$C4[i]) - log(c2/c4))
      } else if (test == "compare_19_10") {
        print("comparing 19*C to 10*C")
        OG_rank_mat[1,1] <- ifelse(log(rank_data$C2[i]/rank_data$C5[i]) == 0, 0, log(rank_data$C2[i]/rank_data$C5[i]) - log(c2/c5))
      } else if (test == "compare_16_13") {
        print("comparing 16*C to 13*C")
        OG_rank_mat[1,1] <- ifelse(log(rank_data$C3[i]/rank_data$C4[i]) == 0, 0, log(rank_data$C3[i]/rank_data$C4[i]) - log(c3/c4))
      } else if (test == "compare_16_10") {
        print("comparing 16*C to 10*C")
        OG_rank_mat[1,1] <- ifelse(log(rank_data$C3[i]/rank_data$C5[i]) == 0, 0, log(rank_data$C3[i]/rank_data$C5[i]) - log(c3/c5))
      } else if (test == "compare_13_10") {
        print("comparing 13*C to 10*C")
        OG_rank_mat[1,1] <- ifelse(log(rank_data$C4[i]/rank_data$C5[i]) == 0, 0, log(rank_data$C4[i]/rank_data$C5[i]) - log(c4/c5))
      }
      OG_mat[[i]] <- OG_rank_mat
      rm(OG_rank_mat) # Remove variable for next iteration
      
      # Calculate element wise mean, sd, se
      # Make a 3D array from list of matrices for calculation purposes
      OG_arr <- array(unlist(OG_mat), c(1, 1, nrow(rank_data)))
      # Get mean, sd, se of third dimension (i.e. mean of all experiments)
      OG_mean <- apply(OG_arr, 1:2, mean)
      OG_sd <- apply(OG_arr, 1:2, sd)
      OG_var <- apply(OG_arr, 1:2, var)
      OG_se <- OG_sd / sqrt(nrow(rank_data))
      OG_test_stat <- abs(OG_mean / OG_se)
      
      # Calculate P value from test statistic (T score = mean / se)
      OG_p_vals[ttt, 1] <- 2 * pt(q = OG_test_stat, df = 4, lower.tail = FALSE)
      rm(OG_test_stat)
    } 
  }
}
# P values for chosen phase and salinity treatment compartment comparison
OG_p_vals <- format(OG_p_vals, scientific = FALSE); OG_p_vals
rank_mat_mean 

# Clear space in environment
rm(list=setdiff(ls(), c("prop_dat")))
