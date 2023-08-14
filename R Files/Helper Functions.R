#load necessary libraries
library(raster)

#function to create a raster stack of yearly averages
yearly_average_raster_data <- function(file_subset)
{
  years <- seq(2000, 2019, by = 1)
  
  #create raster stack for each year... then get average for each month
  raster_stack <- list()
  mean_raster_stack <- list()
  for (i in 1:length(years))
  {
    raster_stack[[i]] <- stack(file_subset[grep(paste(years[i], ".asc", sep = ""), file_subset)])
    mean_raster_stack[[i]] <- calc(raster_stack[[i]], fun = mean)
  }
  #put each year's average into a stack
  stack <- stack(mean_raster_stack)
  return(stack)
}

#takes a row vector of extracted BNAM values and returns an ordered and coherent time series
order_BNAM_values <- function(BNAM_df_row)
{
  timeseries <- c()
  years <- seq(2000, 2019, by = 1)
  for (i in 1:length(years))
  {
    for (j in 1:12)
    {
      col <- grep(paste(month.abb[j],years[i], sep = "_"), colnames(BNAM_df_row))
      timeseries <- c(timeseries, BNAM_df_row[,col])
    }
  }
  return(timeseries)
}

#returns vector denoting the months between 2000 and 2019
BNAM_month_year_order <- function()
{
  monthyears <- c()
  years <- seq(2000, 2019, by = 1)
  for (i in 1:length(years))
  {
    for (j in 1:12)
    {
      monthyears <- c(monthyears, paste(01, month.abb[j],years[i], sep = " "))
    }
  }
  return(monthyears)
}

#function to create a raster stack of monthly averages for 2000 to 2019 (and standard deviations)###
monthly_average_stack <- function(file_subset)
{
  #split up by month and create a raster stack for each month... then get average for each month
  raster_stack <- list()
  mean_raster_stack <- list()
  sd_raster_stack <- list()
  for (i in 1:12)
  {
    raster_stack[[i]] <- stack(file_subset[grep(month.abb[i], file_subset)][11:30])
    mean_raster_stack[[i]] <- calc(raster_stack[[i]], fun = mean)
    sd_raster_stack[[i]] <- calc(raster_stack[[i]], fun = sd)
  }
  #put each month's average into a stack
  stack_mean <- stack(mean_raster_stack)
  stack_sd <- stack(sd_raster_stack)
  return(list(stack_mean, stack_sd))
}

###function that returns the average min, max, and range for a set of BNAM files###
average_min_max_range_calc <- function(file_subset, years)
{
  #create lists containing min, max, and range for each year
  min_raster_list <- list()
  max_raster_list <- list()
  range_raster_list <- list()
  for (i in 1:length(years))
  {
    min_raster_list[[i]] <- calc(stack(file_subset[grep(paste(years[i], ".asc", 
                                                              sep = ""), file_subset)]), fun = min)
    max_raster_list[[i]] <- calc(stack(file_subset[grep(paste(years[i], ".asc", 
                                                              sep = ""), file_subset)]), fun = max)
    range_raster_list[[i]] <- max_raster_list[[i]] - min_raster_list[[i]]
  }
  #turn these lists into stacks and get the average for each
  min_stack <- calc(stack(min_raster_list), fun = mean)
  max_stack <- calc(stack(max_raster_list), fun = mean)
  range_stack <- calc(stack(range_raster_list), fun = mean)
  stack <- list(min_stack, max_stack, range_stack)
  return(stack)
}

#function to match BNAM values to the survey data for a particular BNAM raster 
#stack (using month, year, lat, lon)... returns a vector of the matching values
match_BNAM <- function(surveydata, BNAM_raster_stack, longitude, latitude)
{
  #extract the BNAM values, stick into a data frame
  BNAM_df <- as.data.frame(extract(BNAM_raster_stack, cbind(surveydata[longitude], surveydata[latitude])))
  
  #get the month and year combined to use in matching
  month_year_combined <- paste(surveydata$month_name, surveydata$year, sep = "_")
  
  #go through each row of RV survey data, find BNAM match (same month, year, coordinates) if there is one
  BNAM_values <-c()
  for (i in 1:nrow(surveydata))
  {
    #if no match, set as NA
    if(length(grep(month_year_combined[i], names(BNAM_df))) == 0)
    {
      BNAM_values <- c(BNAM_values, NA)
    }
    #if there is a match, grab the value
    else
    {
      BNAM_values <- c(BNAM_values, BNAM_df[i, grep(month_year_combined[i], names(BNAM_df))]) 
    }
  }
  #return the vector of matching values
  return(BNAM_values)
}

#function used by Figures 5 and 6... makes a partial effect plot for a GLMM covariate 
#(with 95% confidence bounds derived using Breheny and Burchett's contrast approach (2017))
#... let a and b be the index for the covariate's terms, and let original_vector be the 
#covariate's original training vector
effects_plot <- function(model, a, b = NULL, original_vector,..., line_color = "black", 
                         ci.color = "light blue", starve = F)
{
  #values to predict across
  x <- seq(min(original_vector), max(original_vector), length.out = 1000)
  #if not a quadratic
  if(is.null(b)==TRUE)
  {
    if(starve == T)
    {
      y <- model@TMB_out@sdr$par.fixed[a]*x
      se <- x*sqrt(model@tracing@parameter_covariance[a,a])       
    }
    else
    {
      y <- model$sd_report$par.fixed[a]*x
      se <- x*sqrt(vcov(model)[a,a])  
    }
  }
  #if quadratic
  else
  {
    if(starve == T)
    {
      y <- model@TMB_out@sdr$par.fixed[a]*x+model@TMB_out@sdr$par.fixed[b]*x^2
      se <- sqrt(x^2*model@tracing@parameter_covariance[a,a]+
                   x^4*model@tracing@parameter_covariance[b,b]+2*x^3*model@tracing@parameter_covariance[a,b]) 
    }
    else
    {
      y <- model$sd_report$par.fixed[a]*x+model$sd_report$par.fixed[b]*x^2
      se <- sqrt(x^2*vcov(model)[a,a]+x^4*vcov(model)[b,b]+2*x^3*vcov(model)[a,b]) 
    }
  }
  #make plot, add shading for CIs, rug, line
  plot(x, y,..., cex = 0.05, col = line_color)
  polygon(c(x, rev(x)), c(y-1.96*se, rev(y+1.96*se)),
          col = ci.color, lty = 0)
  lines(x,y, col = line_color, ...)
  rug(original_vector, ...)
}

#function to turn GMRF or starve predictions/SEs into a raster
raster_create <- function(prediction_df, prediction_column, year_start)
{
  year_stack <- stack()
  for (i in 1:(2020-year_start))
  {
    prediction_df_new <- prediction_df[which(prediction_df$year==(i+year_start-1)),]
    year_raster <- rasterFromXYZ(cbind(prediction_df_new$UTMX, prediction_df_new$UTMY, 
                                       prediction_df_new[,prediction_column]))
    year_stack <- stack(year_stack, year_raster)
  }
  return(year_stack)
}

#function to turn GAM predictions/SEs into a raster
raster_create_gam <- function(prediction_vector, data_df, year_start)
{
  year_stack <- stack()
  for (i in 1:(2020-year_start))
  {
    data_df_new <- data_df[which(data_df$year==(i+year_start-1)),]
    prediction_vector_new <- prediction_vector[which(data_df$year==(i+year_start-1))]
    year_raster <- rasterFromXYZ(cbind(data_df_new$UTMX, data_df_new$UTMY, prediction_vector_new))
    year_stack <- stack(year_stack, year_raster)
  }
  return(year_stack)
}