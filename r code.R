#========================================================
# PROJECT: Crime Trends in Indian States (2001–2012)
# FILE: 01_data_cleaning.R
# AUTHOR: Anushka Das
# PURPOSE: Import + reshape + create panel dataset
#========================================================
# PROJECT PIPELINE
# Phase 1 : Data Import & Panel Construction  ✅ DONE
# Phase 2 : Descriptive Statistics & Visualisation  ← NOW
# Phase 3 : Panel Econometrics (FE/RE models)
# Phase 4 : Robustness Checks
# Phase 5 : Publication Figures & Tables
#========================================================

#========================================================

# Clear environment
rm(list = ls())

# Set working directory
setwd("C:/Users/HP/Desktop/CRIME2-India-Project")

# Load libraries
library(tidyverse)

#--------------------------------------------------------
# Import raw crime dataset
#--------------------------------------------------------

crime_raw <- read_csv("01_data/CRIME2.csv")

# Preview data
head(crime_raw)
glimpse(crime_raw)


# Clean column names
names(crime_raw)


crime_raw <- crime_raw %>%
  rename(
    state_ut = 1,
    crime_head = 2
  )


crime_long <- crime_raw %>%
  pivot_longer(
    cols = `2001`:`2012`,
    names_to = "year",
    values_to = "crime_count"
  )

head(crime_long)

crime_long <- crime_long %>%
  mutate(
    year = as.numeric(year),
    crime_count = as.numeric(crime_count)
  )

crime_long <- crime_long %>%
  mutate(state_id = as.numeric(factor(state_ut)))


write_csv(crime_long, "01_data/crime_panel.csv")

#========================================================
# PHASE 2 : DESCRIPTIVE ANALYSIS
#========================================================
crime <- read_csv("01_data/crime_panel.csv")

crime %>%
  summarise(
    total_observations = n(),
    number_states = n_distinct(state_ut),
    number_crime_types = n_distinct(crime_head),
    mean_crime = mean(crime_count, na.rm = TRUE),
    median_crime = median(crime_count, na.rm = TRUE),
    sd_crime = sd(crime_count, na.rm = TRUE),
    max_crime = max(crime_count, na.rm = TRUE)
  )
summary_stats <- crime_long %>%
  group_by(crime_head) %>%
  summarise(
    mean_crime = mean(crime_count, na.rm = TRUE),
    sd_crime   = sd(crime_count, na.rm = TRUE),
    min_crime  = min(crime_count, na.rm = TRUE),
    max_crime  = max(crime_count, na.rm = TRUE)
  )

write.csv(summary_stats,
          "03_output/tables/summary_stats.csv",
          row.names = FALSE)


crime_state <- crime %>%
  group_by(state_ut) %>%
  summarise(total_crime = sum(crime_count))

ggplot(crime_state, aes(reorder(state_ut, total_crime), total_crime)) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Total Crime by State (2001–2012)",
    x = "State",
    y = "Total Crime"
  ) +
  theme_minimal()

ggsave("03_output/figures/state_comparison.png", width = 8, height = 6)

library(scales)

crime_trend <- crime_long %>%
  group_by(year, crime_head) %>%
  summarise(total = sum(crime_count, na.rm = TRUE))

p1 <- ggplot(crime_trend,
             aes(x = year, y = total, colour = crime_head)) +
  geom_line(linewidth = 1.2) +
  scale_y_log10(labels = comma) +   # ⭐ KEY FIX (log scale)
  labs(title = "Crime Trends by Crime Category (Log Scale)",
       subtitle = "India, 2001–2012",
       x = "Year",
       y = "Total Crime (log scale)",
       colour = "Crime Type") +
  theme_minimal()

ggsave("03_output/figures/crime_trend_india.png",
       plot = p1, width = 9, height = 6)

#========================================================
# PHASE 3 : PANEL ECONOMETRICS
# STEP 3.1 – Create modelling dataset + panel structure
#========================================================

library(plm)

# 1️⃣ Keep only TOTAL crimes against children
crime_panel <- crime_long %>%
  filter(crime_head == "TOTAL CRIMES AGAINST CHILDREN")

# 2️⃣ Check dataset
unique(crime_panel$crime_head)
dim(crime_panel)

# 3️⃣ Quick sanity trend plot

crime_panel %>%
  group_by(year) %>%
  summarise(total = sum(crime_count, na.rm = TRUE)) %>%
  ggplot(aes(year, total)) +
  geom_line(linewidth = 1)
  geom_point(size = 2) +
  labs(title = "Total Crimes Against Children in India") +
  theme_minimal()
  library(scales)
  
  # National trend of crimes against children
  p_trend <- crime_panel %>%
    group_by(year) %>%
    summarise(total = sum(crime_count, na.rm = TRUE)) %>%
    ggplot(aes(x = year, y = total)) +
    geom_line(color = "#1f77b4", linewidth = 1.3) +
    geom_point(color = "#1f77b4", size = 3) +
    scale_y_continuous(labels = comma) +
    labs(
      title = "Trend in Crimes Against Children in India (2001–2012)",
      subtitle = "Total reported cases across all states",
      x = "Year",
      y = "Number of Cases"
    ) +
    theme_minimal(base_size = 13)
  
  # Save figure
  ggsave("03_output/figures/national_child_crime_trend.png",
         plot = p_trend,
         width = 8, height = 5)
  
  
  library(plm)
  library(lmtest)
  library(sandwich)
  library(stargazer)
  
  # Keep only total crimes
  crime_panel <- crime_long %>%
    filter(crime_head == "TOTAL CRIMES AGAINST CHILDREN")
  
  # Convert to panel format
  pdata <- pdata.frame(crime_panel, index = c("state_ut", "year"))
  
  # Log transformation (very important)
  pdata$log_crime <- log(pdata$crime_count + 1)
  
  
  pooled_ols <- plm(log_crime ~ year,
                    data = pdata,
                    model = "pooling")
  
  summary(pooled_ols)
  
  sink("03_output/tables/pooled_ols.txt")
  summary(pooled_ols)
  sink()
  
  fe_model <- plm(log_crime ~ year,
                  data = pdata,
                  model = "within")
  
  summary(fe_model)
  
  sink("03_output/tables/fixed_effects.txt")
  summary(fe_model)
  sink()
  re_model <- plm(log_crime ~ year,
                  data = pdata,
                  model = "random")
  
  summary(re_model)
  
  sink("03_output/tables/random_effects.txt")
  summary(re_model)
  sink()
  re_model <- plm(log_crime ~ year,
                  data = pdata,
                  model = "random")
  
  summary(re_model)
  
  sink("03_output/tables/random_effects.txt")
  summary(re_model)
  sink()
  
  hausman <- phtest(fe_model, re_model)
  hausman
  
  sink("03_output/tables/hausman_test.txt")
  print(hausman)
  sink()
  
  
  stargazer(pooled_ols, fe_model, re_model,
            type = "text",
            out = "03_output/tables/panel_models.txt",
            title = "Panel Regression Results")
  
  
  
  