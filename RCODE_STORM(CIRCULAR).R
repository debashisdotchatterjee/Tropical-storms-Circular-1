# -------------------------
# Load Required Libraries
# -------------------------
library(dplyr)
library(ggplot2)
library(circular)
library(geosphere)
library(evd)
library(tidyr)
library(purrr)
library(viridis)
library(broom)
library(ismev)

# -------------------------
# Load and Prepare Storms Dataset
# -------------------------
data(storms)

storms <- storms %>%
  mutate(
    lat_rad = lat * pi / 180,
    lon_rad = long * pi / 180,
    x = cos(lat_rad) * cos(lon_rad),
    y = cos(lat_rad) * sin(lon_rad),
    z = sin(lat_rad),
    storm_id = paste(name, year, sep = "_"),
    time = as.POSIXct(paste(year, month, day, hour, sep = "-"), format = "%Y-%m-%d-%H")
  ) %>%
  arrange(storm_id, time)

# Compute headings, turning angles, and step lengths
storms <- storms %>%
  group_by(storm_id) %>%
  mutate(
    lat_prev = lag(lat),
    lon_prev = lag(long),
    heading_deg = ifelse(is.na(lat_prev), NA, bearing(cbind(lon_prev, lat_prev), cbind(long, lat))),
    heading_rad = heading_deg * pi / 180,
    turning_rad = heading_rad - lag(heading_rad),
    turning_rad = ifelse(is.na(turning_rad), NA, (turning_rad + pi) %% (2 * pi) - pi),
    r_prev_x = lag(x),
    r_prev_y = lag(y),
    r_prev_z = lag(z),
    dot = x * r_prev_x + y * r_prev_y + z * r_prev_z,
    step_length = ifelse(is.na(dot), NA, acos(pmin(pmax(dot, -1), 1)))
  ) %>%
  ungroup()

# Assign regimes
storms <- storms %>%
  mutate(regime = case_when(
    status == "tropical depression" ~ 1,
    status == "tropical storm" ~ 2,
    status == "hurricane" ~ 3,
    TRUE ~ 4
  ))

# -------------------------
# Diagnostic
# -------------------------
turnings <- na.omit(storms$turning_rad)
print("Diagnostic: Turning Angles Summary")
print(summary(turnings))
print(paste("Number of valid turning angles:", length(turnings)))
print(paste("Unique turning angles:", length(unique(turnings))))
print(paste("Range of turning angles:", paste(range(turnings, na.rm = TRUE), collapse = " to ")))

# -------------------------
# Circular Statistics Helper
# -------------------------
compute_circular_stats <- function(data, min_n = 10) {
  if (length(data) < min_n || all(is.na(data)) || length(unique(na.omit(data))) < 2) {
    return(list(mean = NA, rho = NA, var = NA, error = TRUE))
  }
  circ_data <- circular(data, type = "angles", units = "radians")
  tryCatch({
    list(mean = mean(circ_data), rho = rho.circular(circ_data), var = var(circ_data), error = FALSE)
  }, error = function(e) {
    message("Error in circular stats: ", e$message)
    list(mean = NA, rho = NA, var = NA, error = TRUE)
  })
}

# -------------------------
# Overall Circular Summaries
# -------------------------
heading_stats <- compute_circular_stats(na.omit(storms$heading_rad))
turning_stats <- compute_circular_stats(turnings)

summary_table_headings <- data.frame(
  Statistic = c("Mean Direction (rad)", "Resultant Length", "Circular Variance"),
  Value = c(heading_stats$mean, heading_stats$rho, heading_stats$var)
)
print("Table 1: Directional Summary Statistics for Headings")
print(summary_table_headings)

summary_table_turnings <- data.frame(
  Statistic = c("Mean Direction (rad)", "Resultant Length", "Circular Variance"),
  Value = c(turning_stats$mean, turning_stats$rho, turning_stats$var)
)
print("Table 2: Directional Summary Statistics for Turning Angles")
print(summary_table_turnings)

# -------------------------
# Regime-specific Circular Summaries
# -------------------------
regime_summaries_headings <- storms %>%
  filter(!is.na(heading_rad)) %>%
  group_by(regime) %>%
  summarise(stats = list(compute_circular_stats(heading_rad)),
            Mean_Direction = map_dbl(stats, ~ifelse(.x$error, NA, .x$mean)),
            Resultant_Length = map_dbl(stats, ~ifelse(.x$error, NA, .x$rho)),
            Circular_Variance = map_dbl(stats, ~ifelse(.x$error, NA, .x$var))) %>%
  select(-stats) %>%
  filter(!is.na(Mean_Direction))
print("Table 3: Regime-Specific Directional Summaries for Headings")
print(regime_summaries_headings)

# -------------------------
# Plots: Headings, Turnings, Step Lengths
# -------------------------
# Rose diagram headings
if(length(na.omit(storms$heading_deg))>0){
  ggplot(data.frame(heading_deg = na.omit(storms$heading_deg)), aes(x = heading_deg)) +
    geom_histogram(bins = 36, fill = "steelblue", color = "black", alpha = 0.8) +
    coord_polar(start = 0) +
    theme_minimal() +
    labs(title = "Rose Diagram of Storm Headings", x = "Direction (deg)", y = "Count")
}

# Rose diagram turning angles
if(length(turnings)>0){
  ggplot(data.frame(turning_deg = turnings*180/pi), aes(x = turning_deg)) +
    geom_histogram(bins = 36, fill = "darkgreen", color = "black", alpha = 0.8) +
    coord_polar(start = 0) +
    theme_minimal() +
    labs(title = "Rose Diagram of Turning Angles", x = "Turning (deg)", y = "Count")
}

# Step lengths histogram by regime
if(nrow(storms %>% filter(!is.na(step_length)))>0){
  ggplot(storms %>% filter(!is.na(step_length)), aes(x = step_length, fill = factor(regime))) +
    geom_histogram(bins = 50, position = "dodge") +
    scale_fill_viridis_d() +
    theme_minimal() +
    labs(title = "Distribution of Step Lengths by Regime", x = "Geodesic Distance (rad)", y = "Count")
}

# Trajectories
if(nrow(storms)>0){
  ggplot(storms, aes(x = long, y = lat, group = storm_id, color = wind)) +
    geom_path(alpha = 0.7) +
    scale_color_viridis_c(option = "plasma") +
    theme_minimal() +
    labs(title = "Storm Trajectories Colored by Wind Speed", color = "Wind (knots)")
}

# -------------------------
# Fit Von Mises to Headings
# -------------------------
vm_fit_head <- tryCatch({
  if(length(na.omit(storms$heading_rad))>10){
    mle.vonmises(circular(na.omit(storms$heading_rad)))
  } else NULL
}, error = function(e) NULL)

if(!is.null(vm_fit_head)){
  vm_table_head <- data.frame(
    Parameter = c("Mu (mean direction)", "Kappa (concentration)"),
    Estimate = c(vm_fit_head$mu, vm_fit_head$kappa),
    SE = c(vm_fit_head$se.mu, vm_fit_head$se.kappa)
  )
  print("Table 4: Von Mises Fit for Headings")
  print(vm_table_head)
  
  # Plot density
  theta_seq <- seq(-pi, pi, length.out = 1000)
  vm_density <- dvonmises(theta_seq, vm_fit_head$mu, vm_fit_head$kappa)
  df_vm <- data.frame(theta = theta_seq, density = vm_density)
  ggplot() +
    geom_histogram(data = data.frame(theta = na.omit(storms$heading_rad)),
                   aes(x = theta, y = ..density..), bins = 36,
                   fill = "lightblue", alpha = 0.5) +
    geom_line(data = df_vm, aes(x = theta, y = density), color = "red", size = 1.2) +
    coord_polar() +
    theme_minimal() +
    labs(title = "Fitted Von Mises Density for Headings")
}

# -------------------------
# Fit Wrapped Cauchy to Turning Angles (Fixed)
# -------------------------
wc_table_turn <- data.frame(
  Parameter = c("Mu (location)", "Rho (concentration)"),
  Estimate = NA,
  SE = NA
)

wc_fit_turn <- tryCatch({
  if(length(turnings)>10 & length(unique(turnings))>2 & !all(turnings==0)){
    fit <- mle.wrappedcauchy(circular(turnings))
    if(!is.null(fit$mu) & !is.null(fit$rho) & !is.na(fit$mu) & !is.na(fit$rho)){
      wc_table_turn$Estimate <- c(fit$mu, fit$rho)
      wc_table_turn$SE <- c(fit$se.mu, fit$se.rho)
    }
    fit
  } else NULL
}, error=function(e) NULL)

print("Table 5: Wrapped Cauchy Fit for Turning Angles")
print(wc_table_turn)

# Plot Wrapped Cauchy density if fit exists
if(!is.null(wc_fit_turn)){
  theta_seq <- seq(-pi, pi, length.out = 1000)
  wc_density <- dwrappedcauchy(theta_seq, wc_fit_turn$mu, wc_fit_turn$rho)
  df_wc <- data.frame(theta = theta_seq, density = wc_density)
  ggplot() +
    geom_histogram(data = data.frame(theta = turnings), aes(x = theta, y = ..density..),
                   bins = 36, fill = "lightgreen", alpha = 0.5) +
    geom_line(data = df_wc, aes(x = theta, y = density), color = "purple", size = 1.2) +
    coord_polar() +
    theme_minimal() +
    labs(title = "Fitted Wrapped Cauchy Density for Turning Angles")
}

# -------------------------
# Wind Regression by Regime
# -------------------------
wind_models <- storms %>%
  filter(!is.na(heading_rad)) %>%
  group_by(regime) %>%
  nest() %>%
  mutate(model = map(data, ~lm(wind ~ cos(heading_rad) + sin(heading_rad) +
                                 sin(2*pi*month/12) + cos(2*pi*month/12), data=.)),
         tidy = map(model, broom::tidy)) %>%
  unnest(tidy)
print("Table 6: Wind Intensity Model Coefficients by Regime")
print(wind_models)

# Wind vs Heading plot
ggplot(storms %>% filter(!is.na(heading_deg)), aes(x = heading_deg, y = wind, color = factor(regime))) +
  geom_point(alpha = 0.5) +
  coord_polar() +
  scale_color_viridis_d() +
  theme_minimal() +
  labs(title = "Wind Speed vs Heading by Regime", color = "Regime")

# Wind by directional bins
ggplot(storms, aes(x=cut(heading_deg, breaks=seq(0,360,45)), y=wind, fill=factor(regime))) +
  geom_boxplot() +
  scale_fill_viridis_d() +
  theme_minimal() +
  labs(title="Wind by Directional Bins and Regime", x="Heading Bin (deg)", y="Wind (knots)")

# -------------------------
# Top 10 Wind Exceedances Table
# -------------------------
u <- quantile(storms$wind, 0.9, na.rm=TRUE)
exceed <- storms %>% filter(wind>u & !is.na(heading_rad)) %>% mutate(excess = wind-u)
top_exceed <- exceed %>% arrange(desc(excess)) %>% head(10)
print("Table 11: Top 10 Wind Exceedances")
print(top_exceed)

# -------------------------
# Correlation Table
# -------------------------
corr_matrix <- storms %>% select(step_length, turning_rad, wind) %>% na.omit() %>% cor()
print("Table 12: Correlation Between Step Length, Turning Angle, Wind")
print(corr_matrix)

#  Step length vs turning
ggplot(storms %>% filter(!is.na(step_length) & !is.na(turning_rad)),
       aes(x=step_length, y=turning_rad, color=factor(regime))) +
  geom_point(alpha=0.6) + scale_color_viridis_d() + theme_minimal() +
  labs(title=" Step Length vs Turning by Regime")

#  Step length vs wind
ggplot(storms %>% filter(!is.na(step_length) & !is.na(wind)),
       aes(x=step_length, y=wind, color=factor(regime))) +
  geom_point(alpha=0.6) + scale_color_viridis_d() + theme_minimal() +
  labs(title="Step Length vs Wind by Regime")

# Table 7: GPD fit for extreme winds
u <- quantile(storms$wind,0.9,na.rm=TRUE)
exceed <- storms %>% filter(wind>u & !is.na(heading_rad)) %>% mutate(excess=wind-u)
if(nrow(exceed)>10){
  gpd_fit <- tryCatch({gpd.fit(exceed$excess, threshold=0, show=FALSE)}, error=function(e) NULL)
  if(!is.null(gpd_fit)){
    table7 <- data.frame(Parameter=c("Scale","Shape"), Estimate=gpd_fit$mle)
  } else {
    table7 <- data.frame(Parameter=c("Scale","Shape"), Estimate=NA)
  }
} else {
  table7 <- data.frame(Parameter=c("Scale","Shape"), Estimate=NA)
}
print("Table 7: GPD Fit")
print(table7)

# Table 10: Regime-wise correlation
table10 <- storms %>% filter(!is.na(step_length) & !is.na(turning_rad) & !is.na(wind)) %>%
  group_by(regime) %>% summarise(corr_step_turn=cor(step_length,turning_rad),
                                 corr_step_wind=cor(step_length,wind),
                                 corr_turn_wind=cor(turning_rad,wind))
table10
