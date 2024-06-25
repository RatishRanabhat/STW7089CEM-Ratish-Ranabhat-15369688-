setwd = ("C:\\Users\\ranar\\OneDrive\\Desktop\\R")
#install.packages("readxl")
#install.packages("ggplot2")
#install.packages("dplyr")
#install.packages("tidyr")
library(readxl) # library for reading excel files 
library(dplyr)   # library for data manipulation
library(ggplot2) # library for plotting 
library(tidyr)   # library for data reshaping 

# Loading the Excel file format Data set of gene expression  
gene_data <- read_excel("C:\\Users\\ranar\\OneDrive\\Desktop\\R\\gene_data.csv_1718251756549.xlsx")

# Converting values of genes into numeric if not 
gene_data$x1 <- as.numeric(gene_data$x1)
gene_data$x2 <- as.numeric(gene_data$x2)
gene_data$x3 <- as.numeric(gene_data$x3)
gene_data$x4 <- as.numeric(gene_data$x4)
gene_data$x5 <- as.numeric(gene_data$x5)
str(gene_data)

# 1. Time series plots

# Plot each gene against time individual

#  For x1 against time 
ggplot(gene_data, aes(x = `Time (min)`)) +
  geom_line(aes(y = x1), color = "blue")+
  labs(title = "Gene Expression over Time",
       x = "Time (minutes)", y = "x1") +
  theme_minimal()

#  For x2 against time 
ggplot(gene_data, aes(x = `Time (min)`)) +
  geom_line(aes(y = x2), color = "red")+
  labs(title = "Gene Expression over Time",
       x = "Time (minutes)", y = "x2") +
  theme_minimal()

#  For x3 against time 
ggplot(gene_data, aes(x = `Time (min)`)) +
  geom_line(aes(y = x3), color = "green")+
  labs(title = "Gene Expression over Time",
       x = "Time (minutes)", y = "x3") +
  theme_minimal()

#  For x4 against time 
ggplot(gene_data, aes(x = `Time (min)`)) +
  geom_line(aes(y = x4), color = "purple")+
  labs(title = "Gene Expression over Time",
       x = "Time (minutes)", y = "x4") +
  theme_minimal()

#  For x5 against time 
ggplot(gene_data, aes(x = `Time (min)`)) +
  geom_line(aes(y = x5), color = "orange")+
  labs(title = "Gene Expression over Time",
       x = "Time (minutes)", y = "x5") +
  theme_minimal()


# Combined Time series Plot for each gene against time in Same graph
ggplot(gene_data, aes(x = `Time (min)`)) +
  geom_line(aes(y = x1, color = "Gene 1"), linewidth = 0.25) +
  geom_line(aes(y = x2, color = "Gene 2"), linewidth = 0.25) +
  geom_line(aes(y = x3, color = "Gene 3"), linewidth = 0.25) +
  geom_line(aes(y = x4, color = "Gene 4"), linewidth = 0.25) +
  geom_line(aes(y = x5, color = "Gene 5"), linewidth = 0.25) +
  labs(title = "Gene Expression over Time",
       x = "Time (minutes)", y = "Expression") +
  scale_color_manual(values = c("blue", "red", "green", "purple", "orange"),
                     labels = c("Gene 1", "Gene 2", "Gene 3", "Gene 4", "Gene 5")) +
  theme_minimal()


# Reshaping data gene_data from wide format to long format 
gene_data_long <- gene_data %>%
  pivot_longer(cols = -`Time (min)`, names_to = "Gene", values_to = "Expression")
gene_data_long

# Plotting Facet Time series Plot for all genes
ggplot(gene_data_long, aes(x = `Time (min)`, y = Expression, color = Gene, group = Gene)) +
  geom_line() +
  labs(x = "Time (minutes)", y = "Expression", title = "Gene Expression Time Series") +
  theme_classic() +
  facet_wrap(~ Gene, scales = "free")+
  theme(panel.background = element_rect(fill = "transparent"))

##############


# 2. Distribution for each gene

# Plotting histograms for gene x1
ggplot(gene_data, aes(x = x1)) +
  geom_histogram(bins = 30, color = "black", fill = "blue", alpha = 0.6) +
  labs(title = "Distribution of Gene x1 Expression", x = "Expression (x1)", y = "Frequency")

# Plotting histograms for gene x2
ggplot(gene_data, aes(x = x2)) +
  geom_histogram(bins = 30, color = "black", fill = "red", alpha = 0.6) +
  labs(title = "Distribution of Gene x2 Expression", x = "Expression (x2)", y = "Frequency")

# Plotting histograms for gene x3
ggplot(gene_data, aes(x = x3)) +
  geom_histogram(bins = 30, color = "black", fill = "green", alpha = 0.6) +
  labs(title = "Distribution of Gene x3 Expression", x = "Expression (x3)", y = "Frequency")

# Plotting histograms for gene x4
ggplot(gene_data, aes(x = x4)) +
  geom_histogram(bins = 30, color = "black", fill = "purple", alpha = 0.6) +
  labs(title = "Distribution of Gene x4 Expression", x = "Expression (x4)", y = "Frequency")

# Plotting histograms for gene x5
ggplot(gene_data, aes(x = x5)) +
  geom_histogram(bins = 30, color = "black", fill = "orange", alpha = 0.6) +
  labs(title = "Distribution of Gene x5 Expression", x = "Expression (x5)", y = "Frequency")

# Combine histograms into one plot with facets
ggplot(gene_data_long, aes(x = Expression, fill = Gene)) +
  geom_histogram(bins = 30, color = "black", alpha = 0.6) +
  labs(x = "Expression", y = "Frequency", title = "Distribution of Gene Expression") +
  facet_wrap(~ Gene, scales = "free") +
  theme_classic() +
  theme(panel.background = element_rect(fill = "transparent"))

# Reshape data to long format for easier plotting
gene_data_lo <- pivot_longer(gene_data, cols = starts_with("x"), names_to = "Genes", values_to = "Expressions")

# plotting Density plots for each gene's expression
ggplot(gene_data_lo, aes(x = Expressions, fill = Genes)) +
  geom_density(color = "black", alpha = 0.6) +
  labs(x = "Expressions", y = "Density", title = "Density Plot of Gene Expressions") +
  facet_wrap(~ Genes, scales = "free") +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "transparent"))

##################################

#3. Correlation and scatter plots
# Compute correlations between genes x1 to x5
library(ggcorrplot)
cor_matrix <- cor(gene_data[, c("x1", "x2", "x3", "x4", "x5")])

# Plot scatter plots for pairs of genes
pairs(gene_data[, c("x1", "x2", "x3", "x4", "x5")])

# Print correlation matrix
print(cor_matrix)

# Scatter plots between x1 and other genes
# Scatter plot: x1 vs x2
ggplot(gene_data, aes(x = x1, y = x2)) +
  geom_point(color = "red") + 
  labs(x = "x1 Expression", y = "x2 Expression", title = "Scatter Plot: x1 vs x2")

# Scatter plot: x1 vs x3
ggplot(gene_data, aes(x = x1, y = x3)) +
  geom_point(color = "blue") + 
  labs(x = "x1 Expression", y = "x3 Expression", title = "Scatter Plot: x1 vs x3")

# Scatter plot: x1 vs x4
ggplot(gene_data, aes(x = x1, y = x4)) +
  geom_point(color = "green") + 
  labs(x = "x1 Expression", y = "x4 Expression", title = "Scatter Plot: x1 vs x4")

# Scatter plot: x1 vs x5
ggplot(gene_data, aes(x = x1, y = x5)) +
  geom_point(color = "purple") + 
  labs(x = "x1 Expression", y = "x5 Expression", title = "Scatter Plot: x1 vs x5")

# Scatter plots between x2 and other genes
# Scatter plot: x2 vs x1
ggplot(gene_data, aes(x = x2, y = x1)) +
  geom_point(color = "red") + 
  labs(x = "x2 Expression", y = "x1 Expression", title = "Scatter Plot: x2 vs x1")

# Scatter plot: x2 vs x3
ggplot(gene_data, aes(x = x2, y = x3)) +
  geom_point(color = "blue") + 
  labs(x = "x2 Expression", y = "x3 Expression", title = "Scatter Plot: x2 vs x3")

# Scatter plot: x2 vs x4
ggplot(gene_data, aes(x = x2, y = x4)) +
  geom_point(color = "green") + 
  labs(x = "x2 Expression", y = "x4 Expression", title = "Scatter Plot: x2 vs x4")

# Scatter plot: x2 vs x5
ggplot(gene_data, aes(x = x2, y = x5)) +
  geom_point(color = "purple") + 
  labs(x = "x2 Expression", y = "x5 Expression", title = "Scatter Plot: x2 vs x5")

# Scatter plots between x3 and other genes
# Scatter plot: x3 vs x1
ggplot(gene_data, aes(x = x3, y = x1)) +
  geom_point(color = "red") + 
  labs(x = "x3 Expression", y = "x1 Expression", title = "Scatter Plot: x3 vs x1")

# Scatter plot: x3 vs x2
ggplot(gene_data, aes(x = x3, y = x2)) +
  geom_point(color = "blue") + 
  labs(x = "x3 Expression", y = "x2 Expression", title = "Scatter Plot: x3 vs x2")

# Scatter plot: x3 vs x4
ggplot(gene_data, aes(x = x3, y = x4)) +
  geom_point(color = "green") + 
  labs(x = "x3 Expression", y = "x4 Expression", title = "Scatter Plot: x3 vs x4")

# Scatter plot: x3 vs x5
ggplot(gene_data, aes(x = x3, y = x5)) +
  geom_point(color = "purple") + 
  labs(x = "x3 Expression", y = "x5 Expression", title = "Scatter Plot: x3 vs x5")

# Scatter plots between x4 and other genes
# Scatter plot: x4 vs x1
ggplot(gene_data, aes(x = x4, y = x1)) +
  geom_point(color = "red") + 
  labs(x = "x4 Expression", y = "x1 Expression", title = "Scatter Plot: x4 vs x1")

# Scatter plot: x4 vs x2
ggplot(gene_data, aes(x = x4, y = x2)) +
  geom_point(color = "blue") + 
  labs(x = "x4 Expression", y = "x2 Expression", title = "Scatter Plot: x4 vs x2")

# Scatter plot: x4 vs x3
ggplot(gene_data, aes(x = x4, y = x3)) +
  geom_point(color = "green") + 
  labs(x = "x4 Expression", y = "x3 Expression", title = "Scatter Plot: x4 vs x3")

# Scatter plot: x4 vs x5
ggplot(gene_data, aes(x = x4, y = x5)) +
  geom_point(color = "purple") + 
  labs(x = "x4 Expression", y = "x5 Expression", title = "Scatter Plot: x4 vs x5")

# Scatter plots between x5 and other genes
# Scatter plot: x5 vs x1
ggplot(gene_data, aes(x = x5, y = x1)) +
  geom_point(color = "red") + 
  labs(x = "x5 Expression", y = "x1 Expression", title = "Scatter Plot: x5 vs x1")

# Scatter plot: x5 vs x2
ggplot(gene_data, aes(x = x5, y = x2)) +
  geom_point(color = "blue") + 
  labs(x = "x5 Expression", y = "x2 Expression", title = "Scatter Plot: x5 vs x2")

# Scatter plot: x5 vs x3
ggplot(gene_data, aes(x = x5, y = x3)) +
  geom_point(color = "green") + 
  labs(x = "x5 Expression", y = "x3 Expression", title = "Scatter Plot: x5 vs x3")

# Scatter plot: x5 vs x4
ggplot(gene_data, aes(x = x5, y = x4)) +
  geom_point(color = "purple") + 
  labs(x = "x5 Expression", y = "x4 Expression", title = "Scatter Plot: x5 vs x4")

# Plot correlation matrix as a circle plot
ggcorrplot(cor_matrix, method = "circle", lab = TRUE,
           title = "Pairwise Pearson Correlation Between Genes")
          
