library(dlm)
library(ggplot2)
library(gplots)
library(forecast)
library(legion)
library(vars)

setwd("/Users/siddharthhaveliwala/Documents/Spring 2024/Courses/CDA 500 - Time Series Analysis/Presentation")
normalized.db <- read.csv('normalized_data.csv')
original.db <- read.csv('original_data.csv')

num_rows <- nrow(original.db)
date_sequence <- seq(0, length.out = num_rows)

# Add the time column to your dataframe
original.db$time <- date_sequence

# Plot fluorescence
plot(original.db$time, original.db$fluorescence, type = "l", xlab = "Time", ylab = "Fluorescence", main = "Fluorescence Time Series")

# Plot Flare_time
plot(original.db$time, original.db$Flare_time, type = "l", xlab = "Time", ylab = "Flare_time", main = "Flare_time Time Series")

# Plot XRAY
plot(original.db$time, original.db$XRAY, type = "l", xlab = "Time", ylab = "XRAY", main = "XRAY Time Series")

# Plot ACF for fluorescence
acf(original.db$fluorescence, main = "ACF for Fluorescence")

# Plot ACF for Flare_time
acf(original.db$Flare_time, main = "ACF for Flare_time")

# Plot ACF for XRAY
acf(original.db$XRAY, main = "ACF for XRAY")

# Apply exponential smoothing to each column
smoothed_exp_fluorescence <- HoltWinters(original.db$fluorescence, beta = FALSE, gamma = FALSE, alpha=0.3)
smoothed_exp_Flare_time <- HoltWinters(original.db$Flare_time, beta = FALSE, gamma = FALSE, alpha=0.3)
smoothed_exp_XRAY <- HoltWinters(original.db$XRAY, beta = FALSE, gamma = FALSE, alpha=0.3)

# Extract the smoothed values
smoothed_values_fluorescence <- fitted(smoothed_exp_fluorescence)
smoothed_values_Flare_time <- fitted(smoothed_exp_Flare_time)
smoothed_values_XRAY <- fitted(smoothed_exp_XRAY)

# Extract the time vector
time_vector <- original.db$time

# Trim the time vector to match the length of smoothed values vector
time_vector <- time_vector[-length(time_vector)]

# Plot the original and exponential smoothed data
par(mfrow = c(3, 1))  # Arrange plots in a 3x1 grid

# Fluorescence
plot(time_vector, original.db$fluorescence[-length(original.db$fluorescence)], type = "l", main = "Original vs Exponential Smoothed Fluorescence Data", xlab = "Time", ylab = "Fluorescence", col='gray')
lines(time_vector, smoothed_values_fluorescence[ , 2], col = "red")

# Flare_time
plot(time_vector, original.db$Flare_time[-length(original.db$Flare_time)], type = "l", main = "Original vs Exponential Smoothed Flare_time Data", xlab = "Time", ylab = "Flare_time", col='gray')
lines(time_vector, smoothed_values_Flare_time[ , 2], col = "red")

# XRAY
plot(time_vector, original.db$XRAY[-length(original.db$XRAY)], type = "l", main = "Original vs Exponential Smoothed XRAY Data", xlab = "Time", ylab = "XRAY", col='gray')
lines(time_vector, smoothed_values_XRAY[ , 2], col = "red")

original.db <- original.db[-1, ]

original.db$fluorescence <- smoothed_values_fluorescence[ , 2]
original.db$Flare_time <- smoothed_values_Flare_time[ , 2]
original.db$XRAY <- smoothed_values_XRAY[ , 2]
original.db <- subset(original.db, select = -c(Date))


# Create a correlation matrix
correlation_matrix <- cor(original.db)

heatmap.2(correlation_matrix, 
          col = colorRampPalette(c("blue", "white", "red"))(100),  # Color gradient
          symm = TRUE,  # Show symmetrically
          margins = c(10, 10),  # Adjust margins
          main = "Correlation Heatmap",  # Title
          key = TRUE,  # Show color key
          keysize = 1.5,  # Size of color key
          trace = "none",  # Turn off trace lines
          density.info = "none",  # Turn off density plot
          cellnote = correlation_matrix,  # Add correlation values
          notecol = "black",  # Color of annotation text
          notecex = 0.8)  # Size of annotation text

# Create data frames for each target variable
flare_time_pred_df <- data.frame(variable = names(original.db))
fluorescence_pred_df <- data.frame(variable = names(original.db))
XRAY_pred_df <- data.frame(variable = names(original.db))

# Compute the correlation matrix
correlation_matrix <- cor(original.db)

# Get variables with correlation above +-0.05 for Flare_time
variables_above_threshold <- names(which(abs(correlation_matrix["Flare_time", ]) > 0.05))

# Extract data for Flare_time_pred_df
flare_time_pred_df <- original.db[, c("time", variables_above_threshold)]

# Repeat the process for fluorescence and XRAY

# Get variables with correlation above +-0.05 for Fluorescence
variables_above_threshold <- names(which(abs(correlation_matrix["fluorescence", ]) > 0.05))

# Extract data for Fluorescence_pred_df
fluorescence_pred_df <- original.db[, c("time", variables_above_threshold)]

# Get variables with correlation above +-0.05 for XRAY
variables_above_threshold <- names(which(abs(correlation_matrix["XRAY", ]) > 0.05))

# Extract data for XRAY_pred_df
XRAY_pred_df <- original.db[, c("time", variables_above_threshold)]

XRAY_pred_df <- subset(XRAY_pred_df, select = -c(time, fluorescence))
fluorescence_pred_df <- subset(fluorescence_pred_df, select = -c(time, XRAY, Flare_time))
flare_time_pred_df <- subset(flare_time_pred_df, select = -c(time, fluorescence))


fluorescence_ts <- ts(fluorescence_pred_df$fluorescence)
other_variables_ts <- lapply(fluorescence_pred_df[-length(fluorescence_pred_df)], ts)
# Compute lagged correlation between fluorescence and other variables
lagged_correlations <- lapply(other_variables_ts, function(x) {
  ccf_result <- ccf(fluorescence_ts, x, lag.max = 10)
  if (length(ccf_result$acf) < 21) {
    warning(paste("Insufficient observations for variable:", deparse(substitute(x))))
    return(rep(NA, 21))
  } else {
    return(ccf_result$acf)
  }
})

# Define lag values
lags <- seq(-10, 10)

# Plot lagged correlations
matplot(lags, sapply(lagged_correlations, function(x) x), type = "l", xlab = "Lag", ylab = "Correlation", main = "Lagged Correlation with Fluorescence")

# Add legend
legend("topright", legend = names(lagged_correlations), col = 1:length(lagged_correlations), lty = 1:length(lagged_correlations), cex = 0.8)

# Convert dataframe to time series
flare_time_ts <- ts(flare_time_pred_df$Flare_time)
other_variables_ts <- lapply(flare_time_pred_df[-length(flare_time_pred_df)], ts)

# Compute lagged correlation between flare_time and other variables
lagged_correlations <- lapply(other_variables_ts, function(x) {
  ccf_result <- ccf(flare_time_ts, x, lag.max = 10)
  if (length(ccf_result$acf) < 21) {
    warning(paste("Insufficient observations for variable:", deparse(substitute(x))))
    return(rep(NA, 21))
  } else {
    return(ccf_result$acf)
  }
})

# Define lag values
lags <- seq(-10, 10)

# Plot lagged correlations
matplot(lags, sapply(lagged_correlations, function(x) x), type = "l", xlab = "Lag", ylab = "Correlation", main = "Lagged Correlation with Flare_time")

# Add legend
legend("topright", legend = names(lagged_correlations), col = 1:length(lagged_correlations), lty = 1:length(lagged_correlations), cex = 0.8)

# Convert dataframe to time series
XRAY_ts <- ts(XRAY_pred_df$XRAY)
other_variables_ts <- lapply(XRAY_pred_df[-length(XRAY_pred_df)], ts)

# Compute lagged correlation between XRAY and other variables
lagged_correlations <- lapply(other_variables_ts, function(x) {
  ccf_result <- ccf(XRAY_ts, x, lag.max = 10)
  if (length(ccf_result$acf) < 21) {
    warning(paste("Insufficient observations for variable:", deparse(substitute(x))))
    return(rep(NA, 21))
  } else {
    return(ccf_result$acf)
  }
})

# Define lag values
lags <- seq(-10, 10)

# Plot lagged correlations
matplot(lags, sapply(lagged_correlations, function(x) x), type = "l", xlab = "Lag", ylab = "Correlation", main = "Lagged Correlation with XRAY")

# Add legend
legend("topright", legend = names(lagged_correlations), col = 1:length(lagged_correlations), lty = 1:length(lagged_correlations), cex = 0.8)

# Function to perform min-max normalization
min_max_normalization <- function(df) {
  normalized_df <- as.data.frame(lapply(df, function(x) {
    (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  }))
  return(normalized_df)
}


# Function to fit ARIMA models and show diagnostics
# Function to fit ARIMA models and find the best model based on RMSE
fit_best_arima_and_diagnostics <- function(dataframe, variable_name) {
  
  predictor_variables <- dataframe[-length(dataframe)]
  predictor_variables <- predictor_variables[-length(predictor_variables)]
  
  predictor_variables <- min_max_normalization(predictor_variables)
  model <- auto.arima(dataframe[[variable_name]], xreg = as.matrix(predictor_variables))
    
  
  # Show model diagnostics for the best model
  cat("ARIMA Model Diagnostics for", variable_name, ":\n")
  checkresiduals(model)
  print(summary(model))
  # Return the best model
  return(model)
}

# Apply ARIMA models and show diagnostics for flare_time_pred_df
flare_time_arima_model <- fit_best_arima_and_diagnostics(flare_time_pred_df, "Flare_time")

# Apply ARIMA models and show diagnostics for fluorescence_pred_df
fluorescence_arima_model <- fit_best_arima_and_diagnostics(fluorescence_pred_df, "fluorescence")

# Apply ARIMA models and show diagnostics for XRAY_pred_df
XRAY_arima_model <- fit_best_arima_and_diagnostics(XRAY_pred_df, "XRAY")



fit_best_ets_and_diagnostics <- function(dataframe, variable_name) {
  # Extract predictor variables and perform min-max normalization
  dataframe[-ncol(dataframe)] <- min_max_normalization(dataframe[-ncol(dataframe)])
  
  # Fit ETS model
  model <- VAR(dataframe)
  
  print(summary(model))
  
  # Return the best model
  return(model)
}

# Apply VAR models and show diagnostics for flare_time_pred_df
flare_time_model <- fit_best_ets_and_diagnostics(flare_time_pred_df, "Flare_time")

# Apply VAR models and show diagnostics for fluorescence_pred_df
fluorescence_model <- fit_best_ets_and_diagnostics(fluorescence_pred_df, "fluorescence")

# Apply VAR models and show diagnostics for XRAY_pred_df
XRAY_model <- fit_best_ets_and_diagnostics(XRAY_pred_df, "XRAY")


