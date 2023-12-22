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

#used by Figures 2 and 3
#plots partial effects for all three models at once (i.e., overlapping)...
#let a and be the index for the variable's terms - in the starve and sdmTMB case, 
#this must match with the covariate's row numbers in the model variance-covariance matrix...
#let original_vector be the variable's original training vector
#let sig = a vector containing the names of the models (i.e., "starve", "sdmTMB", and/or "mgcv") 
#for which the term is significant (to add a star)
effects_plot_combined <- function(mgcv_model, starve_model, sdmTMB_model, mgcv_a, starve_a, starve_b = NULL, 
                                  sdmTMB_a, sdmTMB_b = NULL, original_vector, ylab = "Partial Effect", cex.main = 2,
                                  cex.axis = 1.6, cex.lab = 1.6, lwd = 3, ... ,
                                  line_colors = c("blue", "black","red"), sigstar_colors = line_colors,
                                  ci.colors = c("light blue", 
                                                adjustcolor("black", alpha.f=0.17), 
                                                adjustcolor("red", alpha.f=0.17)),
                                  sig = c("starve","sdmTMB","mgcv"))
{
  #mgcv model plot
  plot(mgcv_model, select = mgcv_a, shade = TRUE,  se = 1.96, lwd = lwd, rug = T, col = line_colors[1], 
       shade.col = ci.colors[1], ..., ylab = ylab, cex.main = cex.main, cex.axis = cex.axis, cex.lab = cex.lab)
  
  # Add rug plot to the x-axis
  rug(original_vector, side = 1, lwd = lwd)
  
  #values to predict across for other two models
  x <- seq(min(original_vector), max(original_vector), length.out = 1000)
  #if quadratic isn't allowed
  if(is.null(starve_b)&&is.null(sdmTMB_b))
  {
    #calculate starve model values using Breheny and Burchett's contrast approach (2017)
    y_starve <- starve_model@TMB_out@sdr$par.fixed[starve_a]*x
    se_starve <- x*sqrt(starve_model@tracing@parameter_covariance[starve_a,starve_a])       
    
    #calculate sdmTMB model values using Breheny and Burchett's contrast approach (2017)
    y_sdmTMB <- sdmTMB_model$sd_report$par.fixed[sdmTMB_a]*x
    se_sdmTMB <- x*sqrt(vcov(sdmTMB_model)[sdmTMB_a,sdmTMB_a])  
  }
  #if quadratic is allowed
  else
  {
    #calculate starve model values using Breheny and Burchett's contrast approach (2017)
    y_starve <- starve_model@TMB_out@sdr$par.fixed[starve_a]*x+starve_model@TMB_out@sdr$par.fixed[starve_b]*x^2
    se_starve <- sqrt(x^2*starve_model@tracing@parameter_covariance[starve_a,starve_a]+x^4*
                        starve_model@tracing@parameter_covariance[starve_b, starve_b]+2*x^3*
                        starve_model@tracing@parameter_covariance[starve_a,starve_b]) 
    
    #calculate sdmTMB model values using Breheny and Burchett's contrast approach (2017)
    y_sdmTMB <- sdmTMB_model$sd_report$par.fixed[sdmTMB_a]*x+sdmTMB_model$sd_report$par.fixed[sdmTMB_b]*x^2
    se_sdmTMB <- sqrt(x^2*vcov(sdmTMB_model)[sdmTMB_a,sdmTMB_a]+x^4*vcov(sdmTMB_model)[sdmTMB_b,sdmTMB_b]+2*x^3*
                        vcov(sdmTMB_model)[sdmTMB_a,sdmTMB_b]) 
  }
  
  #add starve model estimated effect and uncertainty to plot
  polygon(c(x, rev(x)), c(y_starve-1.96*se_starve, rev(y_starve+1.96*se_starve)),
          col = ci.colors[2], lty = 0)
  lines(x,y_starve, lwd = lwd, col = line_colors[2])
  
  #add sdmTMB model estimated effect and uncertainty to plot
  polygon(c(x, rev(x)), c(y_sdmTMB-1.96*se_sdmTMB, rev(y_sdmTMB+1.96*se_sdmTMB)),
          col = ci.colors[3], lty = 0)
  lines(x,y_sdmTMB, lwd = lwd, col = line_colors[3])
  
  #add significance stars
  if("mgcv"%in%sig)
  {
    mtext("*", col = sigstar_colors[1], adj = 0.15, cex = cex.main)
  }
  if("starve"%in%sig)
  {
    mtext("*", col = sigstar_colors[2], adj = 0.5, cex = cex.main)
  }
  if("sdmTMB"%in%sig)
  {
    mtext("*", col = sigstar_colors[3], adj = 0.85, cex = cex.main)
  }
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