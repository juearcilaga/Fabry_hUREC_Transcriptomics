# Script for Visualizing Estimated Glomerular Filtration Rate (eGFR) Over Time
# Author: [Your Name]
# Date: [Today's Date]
# Description: This script processes eGFR data, visualizes trends over time, 
# and applies linear regression analysis to different treatment phases.
#
# Inputs:
# - A dataset containing dates and eGFR values.
#
# Outputs:
# - A scatter plot with event markers and regression lines, saved as 'eGFR.svg'.

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Define the dates of measurements
dates <- as.Date(c("6-05-2022","11-05-2022","11-05-2022","13-06-2022",
                   "27-06-2022","11-07-2022","27-07-2022","10-08-2022",
                   "20-10-2022","05-01-2023","12-04-2023","12-04-2023",
                   "03-07-2023","03-07-2023","04-01-2024", "24-07-2023","08-10-2023",
                   "26-10-2023","18-12-2023","04-04-2024","10-04-2024",
                   "30-04-2024","03-05-2024","11-06-2024","10-07-2024",
                   "22-07-2024","02-09-2024","06-09-2024","26-09-2024",
                   "28-11-2024","02-12-2024"), format = "%d-%m-%Y")

# Define eGFR values
gfr <- c("36","47","49","50",
         "51","48","45","48",
         "41","39","39","",
         "36","", "41","37",
         "34","39","36","31",
         "37","31","30","31",
         "27","29","31","29",
         "28","24", "26")

# Create a data frame
data <- data.frame(date = dates, GFR = gfr)
data <- data %>% filter(GFR != "")

# Define event phases with start and end dates
event_ranges <- data.frame(
  Event = c("Pre-treatment", "Chaperone", "ERT"),
  Start = as.Date(c("6-05-2022", "13-01-2023", "03-08-2024"), format = "%d-%m-%Y"),
  End = as.Date(c("12-01-2023", "31-07-2024", "02-12-2024"), format = "%d-%m-%Y")
)

# Assign event categories
data$Event <- "Other"
data$Event[data$date >= event_ranges$Start[1] & data$date <= event_ranges$End[1]] <- "Pre-treatment"
data$Event[data$date >= event_ranges$Start[2] & data$date <= event_ranges$End[2]] <- "Chaperone"
data$Event[data$date >= event_ranges$Start[3] & data$date <= event_ranges$End[3]] <- "ERT"

data$Event <- as.factor(data$Event)

# Scatter plot with event colors
eGFR_plot <- ggplot(data, aes(x = date, y = as.numeric(GFR), colour = Event)) +
  geom_point(size = 1) +  # Add points
  scale_x_date(
    date_breaks = "2 months",  # Breaks every 2 months
    date_labels = "%m-%Y"      # Format as Year-Month
  ) +
  labs(title = "", 
       x = "Dates", 
       y = "eGFR (mL/min/1.73m^2)") +
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
  do(model = lm(as.numeric(GFR) ~ as.numeric(date), data = .)) %>%
  summarise(
    slope = coef(model)[2],  # Extract slope (rate of change)
    intercept = coef(model)[1],  # Extract intercept
    Event = Event
  )

# Print the regression slopes
print(slopes)

# Save the plot
ggsave(file = "eGFR.svg", plot = eGFR_plot, width = 10, height = 6)
