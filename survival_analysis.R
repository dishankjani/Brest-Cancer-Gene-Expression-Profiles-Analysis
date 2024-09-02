install.packages("survival")
install.packages("survminer")

# Load the packages
library(survival)
library(survminer)
library("dplyr")
library(ggplot2)
library(tidyverse)

metabric_data <- read.csv("METABRIC_RNA_Mutation.csv")
head(metabric_data)

#step 1: Cleaning data
# Count missing values in each column
colSums(is.na(metabric_data))
metabric_data <- na.omit(metabric_data)  # Removes rows with any NA values

threshold <- 0.5  # Set your threshold (e.g., remove columns with more than 50% missing data)
metabric_data <- metabric_data[, colMeans(is.na(metabric_data)) < threshold]

# Assume 'time' is the survival time, 'status' is the event indicator,
# and 'group' is a grouping variable such as a treatment or mutation status.
time <- metabric_data$overall_survival_months
status <- metabric_data$overall_survival
group <- metabric_data$type_of_breast_surgery

# Create a survival object
surv_object <- Surv(time, status)

# If you have a grouping variable, stratify the analysis
#we assign surv_object as 1 if there is ust one group, but here there is multiple groups.
fit <- survfit(surv_object ~ group, data = metabric_data)

summary(metabric_data$overall_survival)

# Basic Kaplan-Meier plot
ggsurvplot(fit, 
           data = metabric_data, 
           pval = TRUE,               # Show p-value
           conf.int = TRUE,           # Show confidence intervals
           risk.table = TRUE,         # Show risk table
           surv.median.line = "hv",   # Add median survival line
           ggtheme = theme_minimal(), # Set theme
           palette = "Dark2")         # Set color palette


#we compare survival between two groups using the log-rank test (also known as Mantel-Cox test). 
#This test is the most common hypothesis test to compare survival between two groups.
survdiff(Surv(time, status) ~ group,
         data = metabric_data
)


#summarize the age at which diagnosis is done
#This analysis helps to better undertand at what age an average person is diagnosed from brest cancer
age_of_diagnosis_summary <- metabric_data %>%
  summarise(n = length(age_at_diagnosis),
            Mean = mean(age_at_diagnosis, na.rm = TRUE),
            Median = median(age_at_diagnosis, na.rm = TRUE),
            SD = sd(age_at_diagnosis, na.rm = TRUE),
            Min = min(age_at_diagnosis, na.rm = TRUE),
            Max = max(age_at_diagnosis, na.rm = TRUE)) 
print(age_of_diagnosis_summary)

#identify average age of people and their overall survival
Table_1 <- ddply(metabric_data, "overall_survival", 
                 summarise, 
                 n = length(age_at_diagnosis), 
                 Mean.age_at_diagnosis = mean(age_at_diagnosis), 
                 Stdev.age_at_diagnosis = sd(age_at_diagnosis),
                 Min = min(age_at_diagnosis),
                 Max = max(age_at_diagnosis))
Table_1


# Perform hierarchical clustering
dist_matrix <- dist(metabric_data$mutation_count)  # Calculate distance matrix
cluster_fit <- hclust(dist_matrix, method = "ward.D2")  # Perform hierarchical clustering

# Cut the dendrogram to form clusters
data$cluster <- cutree(cluster_fit, k = 3)  # Choose the number of clusters (e.g., 3)

# Plot the dendrogram
plot(cluster_fit, labels = FALSE, hang = -1)
rect.hclust(cluster_fit, k = 3, border = "red")  # Highlight clusters with red borders
