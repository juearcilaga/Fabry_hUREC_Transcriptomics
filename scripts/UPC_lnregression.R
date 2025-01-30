# Script for Visualizing Urine Protein Creatinine (UPC) Over Time
# Author: Dr. Juliana E. Galvis
# Date: 18/Dec/2024
# Description: This script processes urine protein creatinine ratio data, visualizes trends over time, 
# and applies linear regression analysis to different treatment phases.
#
# Inputs:
# - A dataset containing dates and urine protein creatinine (UPC) values.
#
# Outputs:
# - A scatter plot with event markers and regression lines, saved as 'UPCR.svg'.

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Define the dates of measurements
dates <- as.Date(c("11-05-2022", "11-05-2022", "13-06-2022", 
                   "27-06-2022", "11-07-2022", "27-07-2022", "10-08-2022",
                   "20-10-2022", "05-01-2023", "12-04-2023", "12-04-2023",
                   "03-07-2023", "03-07-2023", "04-01-2024", "24-07-2023", "08-10-2023",
                   "26-10-2023", "18-12-2023", "04-04-2024", "10-04-2024",
                   "30-04-2024", "03-05-2024", "11-06-2024", "10-07-2024",
                   "22-07-2024", "02-09-2024", "06-09-2024", "26-09-2024",
                   "28-11-2024", "02-12-2024"), format = "%d-%m-%Y")

# Define urine UPC values
urine_UPC <- c("","","","",
               "380","","","",
               "","","226","",
               "249","", "","",
               "","185","249","",
               "","180","","",
               "","","","243",
               "","")

# Create a data frame
data <- data.frame(date = dates, Urine_UPC = urine_UPC)

# Filter out empty values
data <- data %>% filter(Urine_UPC != "")

# Scatter plot of Urine UPC over time
ggplot(data, aes(x = date, y = Urine_UPC)) +
  geom_point(color = "black", size = 1) +  # Add points
  scale_x_date(
    date_breaks = "2 months",  # Breaks every 2 months
    date_labels = "%m-%Y"      # Format as Year-Month
  ) +
  labs(title = "Scatter Plot of Urine UPC Over Time", 
       x = "Dates", 
       y = "Urine UPC Values") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1)  # Rotate x-axis labels
  )

# Define event phases with start and end dates
event_ranges <- data.frame(
  Event = c("Pre-treatment", "Chaperone", "ERT"),
  Start = as.Date(c("6-05-2022", "13-01-2023", "03-08-2024"), format = "%d-%m-%Y"),
  End = as.Date(c("12-01-2023", "31-07-2024", "02-12-2024"), format = "%d-%m-%Y")
)

# Filter data for each event range
range_data_Pretreatment <- data %>% filter(date >= event_ranges$Start[1] & date <= event_ranges$End[1])
range_data_Chaperone <- data %>% filter(date >= event_ranges$Start[2] & date <= event_ranges$End[2])
range_data_ERT <- data %>% filter(date >= event_ranges$Start[3] & date <= event_ranges$End[3])

# Assign event categories
data$Event <- "Other"
data[data$date %in% range_data_Pretreatment$date, ]$Event <- "Pre-treatment"
data[data$date %in% range_data_Chaperone$date, ]$Event <- "Chaperone"
data[data$date %in% range_data_ERT$date, ]$Event <- "ERT"

data$Event <- as.character.factor(data$Event)

# Scatter plot with event colors
eUPC_plot <- ggplot(data, aes(x = date, y = Urine_UPC, colour = Event)) +
  geom_point(color = "black", size = 1) +  # Add points
  scale_x_date(
    date_breaks = "2 months",  # Breaks every 2 months
    date_labels = "%m-%Y"      # Format as Year-Month
  ) +
  labs(title = "Scatter Plot of Urine UPC Over Time", 
       x = "Dates", 
       y = "Urine UPC Values") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1)  # Rotate x-axis labels
  ) +
  geom_vline(xintercept = event_ranges$Start, color = "red", linetype = "dashed")

# Convert date column to numeric for regression analysis
data$date <- as.Date(data$date)

eUPC_plot <- ggplot(data, aes(x = date, y = as.numeric(Urine_UPC), colour = Event)) +
  geom_point(color = "black", size = 1) +  # Add points
  scale_x_date(
    date_breaks = "2 months",  # Breaks every 2 months
    date_labels = "%m-%Y"      # Format as Year-Month
  ) +
  labs(title = "", 
       x = "Dates", 
       y = "Urine Protein Creatinine Ratio") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 11),  # Rotate x-axis labels
    axis.text.y = element_text(size = 11)
  ) +
  geom_vline(xintercept = event_ranges$Start, color = "black", linetype = "dashed") +
  geom_smooth(method = "lm", aes(group = Event), color = "black", se = TRUE)   # Regression lines

# Linear regression analysis for each event group
slopes <- data %>%
  group_by(Event) %>%
  do(model = lm(Urine_UPC ~ as.numeric(date), data = .)) %>%
  summarise(
    slope = coef(model)[2],  # Extract slope (rate of change)
    intercept = coef(model)[1],  # Extract intercept
    Event = Event
  )

# Print the regression slopes
print(slopes)

# Save the plot
ggsave(file = "UPCR.svg", plot = eUPC_plot, width = 10, height = 6)