#load necessary libraries
library(ggplot2); library(starve); library(sf)
library(sp); library(gtools); library(knitr)
library(viridis); library(dplyr); library(reshape)
library(parallel); library(DHARMa); library(raster)
library(ggnewscale); source("R Files/Helper Functions.R")

#increase memory limit
memory.limit(size = 400000)

#set number of decimals to round tables to
round_val <- 3

#set seed for residuals
set.seed(1)

#load in pre-processed survey data
survey_data_combined <- read.csv(paste(getwd(),"/Data/survey_data_combined.csv", sep = ""))
#transform some covariates
covariates <- c("BtmTempBNAM","BtmSalinityBNAM","BtmStressBNAMLog","RangeTemp","RangeStressLog",
                "DEM_log","sqrt_DEM_Slope","DEM_Easterness",
                "DEM_Northerness","DEM_RDMV")
survey_data_combined$RangeStressLog <- log(survey_data_combined$RangeStress)
survey_data_combined$BtmStressBNAMLog <- log(survey_data_combined$BtmStressBNAM)
survey_data_combined$sqrt_DEM_Slope <- sqrt(survey_data_combined$DEM_Slope)
survey_data_combined$DEM_log <- log(-survey_data_combined$DEM)
#scale covariates using all data
survey_data_combined[paste(covariates,"Scaled", sep="")] <- scale(survey_data_combined[covariates])

#grab NS, Quebec, New Brunswick, PEI, NFLD maps for later plots
map <- raster::getData(country = "CAN", level = 1)
map <- map[which(map$NAME_1 %in% c("Nova Scotia", "Prince Edward Island", 
                                   "New Brunswick", "Newfoundland and Labrador", "Québec")),]
map <- spTransform(map, crs("+datum=WGS84 +proj=utm +zone=20T +units=km"))
#grab Maine map for later plots
map1 <- raster::getData(country = "USA", level = 1)
map1 <- map1[which(map1$NAME_1 == "Maine"),]
map1 <- spTransform(map1, crs("+datum=WGS84 +proj=utm +zone=20T +units=km"))

#make duplicate UTMX, UTMY columns as originals will disappear
survey_data_combined$UTMX_Extra <- survey_data_combined$UTMX
survey_data_combined$UTMY_Extra <- survey_data_combined$UTMY
#get dataframe into an sf framework as starve requires
survey_data_combined <- st_as_sf(survey_data_combined, coords = c("UTMX", "UTMY"),
         crs = "+proj=utm +zone=20T +units=km")

#remove RV winter survey tows
#survey_data_combined <- survey_data_combined[-which(survey_data_combined$season_winter == 1),]

#set up node locations so that snowcrab stations aren't repeated
survey_data_combined_RV <- survey_data_combined[which(survey_data_combined$snowcrab==FALSE),]
survey_data_combined_snowcrab <- survey_data_combined[which(survey_data_combined$snowcrab==TRUE),]
survey_data_combined_snowcrab_no_dup <- survey_data_combined_snowcrab[match(unique(survey_data_combined_snowcrab$STATION), 
                                                                                  survey_data_combined_snowcrab$STATION),]
combined_locations <- rbind(survey_data_combined_RV, survey_data_combined_snowcrab_no_dup)
combined_locations <- combined_locations[order(combined_locations$year),]

##########################################ENCOUNTER SUB-MODEL###########################################

#set up response
survey_data_combined$presence <- as.numeric(survey_data_combined$presence)

#fit sub-model with all data... random walk so will converge... nu = 0.5
start_time <- Sys.time()
encounter_starve <- strv_prepare(presence ~ DEM_logScaled+I(DEM_logScaled^2)+
                                 BtmTempBNAMScaled+I(BtmTempBNAMScaled^2)+
                                 RangeTempScaled+I(RangeTempScaled^2)+BtmSalinityBNAMScaled+
                                 I(BtmSalinityBNAMScaled^2)+
                                 BtmStressBNAMLogScaled+RangeStressLogScaled+I(RangeStressLogScaled^2)+
                                 sqrt_DEM_SlopeScaled+I(sqrt_DEM_SlopeScaled^2)+DEM_NorthernessScaled+
                                 I(DEM_NorthernessScaled^2)+DEM_EasternessScaled+DEM_RDMVScaled+
                                 I(DEM_RDMVScaled^2)+as.numeric(snowcrab)+time(year, "rw"), 
                                  nodes = combined_locations, n_neighbours = 5, 
                               data = survey_data_combined[order(survey_data_combined$year),], 
                               distribution = "bernoulli")
encounter_starve <- strv_fit(encounter_starve, silent = F)
end_time <- Sys.time()
fit_time1 <- end_time - start_time #get time to fit

#create table to describe fixed effects, etc. - Table S3
fixed_effects_table <- round(fixed_effects(encounter_starve)$presence[1:2], round_val)
rownames(fixed_effects_table) <- c("Log(Depth)","Log(Depth) Squared","Bottom Temperature", 
                                   "Bottom Temperature Squared","Btm. Temperature Range",
                                   "Btm. Temperature Range Squared","Bottom Salinity", 
                                   "Bottom Salinity Squared", "Log(Bottom Stress)", "Log(Btm. Stress Range)", 
                                   "Log(Btm. Stress Range) Squared", "Sqrt(Slope)","Sqrt(Slope) Squared", 
                                   "Northerness","Northerness Squared","Easterness","RDMV","RDMV Squared",
                                   "Snow Crab Survey")
colnames(fixed_effects_table) <- c("Estimate", "Std. Error")
spatial_parameters_table <- round(space_parameters(encounter_starve)$presence[1:2], round_val)
colnames(spatial_parameters_table) <- c("Estimate", "Std. Error")
rownames(spatial_parameters_table) <- c("SD", "Range","Nu")
time_parameters_table <- round(time_parameters(encounter_starve)$presence[1:2], round_val)
rownames(time_parameters_table) <- c("Mu","AR(1)","SD")
colnames(time_parameters_table) <- c("Estimate", "Std. Error")
kable(list(fixed_effects_table, spatial_parameters_table, time_parameters_table), format = "latex")

#grab training predictions and set up data frame to hold them
prediction_df <- cbind(as.data.frame(data_predictions(encounter_starve)@locations), 
                       as.data.frame(data_predictions(encounter_starve)@values))
coordinates <- do.call(rbind, st_geometry(data_predictions(encounter_starve)@locations))
prediction_df <- cbind(prediction_df, coordinates)
colnames(prediction_df)[which(colnames(prediction_df)=="1"|colnames(prediction_df)=="2")] <- c("X","Y") 
#simulate 250 times to get DHARMa residuals
sim<-do.call(
  cbind,
  mclapply(
    seq(250),
    function(i) {
      sim<-strv_simulate(
        encounter_starve,
        conditional=TRUE
      )
      return(dat(sim)$presence)
    },
    mc.cores=1
  )
)
#get and plot the DHARMa residuals
encounter_residuals<-createDHARMa(
  simulated=sim,
  observed=as.numeric(dat(encounter_starve)$presence),
  fittedPredictedResponse = prediction_df$response,
  integer=TRUE
)
plot(encounter_residuals)

#spatio-temporally plot residuals
ggplot() + geom_point(size = 0.9, aes(x = X, y = Y, col = encounter_residuals$scaledResiduals),
                           data = prediction_df) + 
  scale_colour_gradientn(colours = viridis(100)) +
  labs(col = "Quantile\nResidual") + theme(text=element_text(size=16)) + 
  geom_polygon(data = fortify(map), aes(group = group, x = long, y = lat)) +
  xlab("Eastings") + ylab("Northings") + facet_wrap(~year) + 
  coord_cartesian(ylim = c(4650, 5250), xlim = c(125, 1000)) + 
  geom_polygon(data = fortify(map1), aes(group = group, x = long, y = lat))
#temporally plot residuals
ggplot(aes(x = year, y = encounter_residuals$scaledResiduals), data = prediction_df) + 
  geom_point(size = 0.15) + geom_smooth(method = "loess") + ylab("Quantile Residual") + xlab("Year") + 
  theme(text=element_text(size=16))

###################################CPUE SUB-MODEL################################

#grab encounters and set up response
survey_data_combined_CPUE <- survey_data_combined[which(is.na(survey_data_combined$std.WGT)==FALSE&
                                                          survey_data_combined$std.WGT!=0),]
survey_data_combined_CPUE$std.WGT_log <- log(survey_data_combined_CPUE$std.WGT)

#fit sub-model with all data... random walk so will converge... nu = 0.5
start_time <- Sys.time()
CPUE_starve <- strv_prepare(std.WGT_log ~ DEM_logScaled+I(DEM_logScaled^2)+BtmTempBNAMScaled+
                             I(BtmTempBNAMScaled^2)+RangeTempScaled+
                             I(RangeTempScaled^2)+BtmSalinityBNAMScaled+
                             I(BtmSalinityBNAMScaled^2)+BtmStressBNAMLogScaled+RangeStressLogScaled+
                             I(RangeStressLogScaled^2)+sqrt_DEM_SlopeScaled+I(sqrt_DEM_SlopeScaled^2)+
                             DEM_NorthernessScaled+I(DEM_NorthernessScaled^2)+DEM_EasternessScaled+DEM_RDMVScaled+
                             I(DEM_RDMVScaled^2)+as.numeric(snowcrab)+time(year, "rw"), 
                           nodes = combined_locations, 
                           data = survey_data_combined_CPUE[order(survey_data_combined_CPUE$year),], 
                           n_neighbours = 5, distribution = "gaussian")
CPUE_starve <- strv_fit(CPUE_starve, silent = F)
end_time <- Sys.time()
fit_time2 <- end_time - start_time #get time to fit

#create table to describe fixed effects, etc. - Table S4
fixed_effects_table <- round(fixed_effects(CPUE_starve)$std.WGT_log[1:2], round_val)
rownames(fixed_effects_table) <- c("Log(Depth)","Log(Depth) Squared","Bottom Temperature", 
                                   "Bottom Temperature Squared","Btm. Temperature Range",
                                   "Btm. Temperature Range Squared","Bottom Salinity", 
                                   "Bottom Salinity Squared", "Log(Bottom Stress)", "Log(Btm. Stress Range)", 
                                   "Log(Btm. Stress Range) Squared",
                                   "Sqrt(Slope)","Sqrt(Slope) Squared", "Northerness","Northerness Squared",
                                   "Easterness","RDMV","RDMV Squared", "Snow Crab Survey")
colnames(fixed_effects_table) <- c("Estimate", "Std. Error")
spatial_parameters_table <- round(space_parameters(CPUE_starve)$std.WGT_log[1:2], round_val)
colnames(spatial_parameters_table) <- c("Estimate", "Std. Error")
rownames(spatial_parameters_table) <- c("SD", "Range","Nu")
time_parameters_table <- round(time_parameters(CPUE_starve)$std.WGT_log[1:2], round_val)
rownames(time_parameters_table) <- c("Mu","AR(1)","SD")
colnames(time_parameters_table) <- c("Estimate", "Std. Error")
kable(list(fixed_effects_table, spatial_parameters_table, time_parameters_table), format = "latex")

#grab training predictions and set up data frame to hold them
prediction_df <- cbind(as.data.frame(data_predictions(CPUE_starve)@locations), 
                       as.data.frame(data_predictions(CPUE_starve)@values))
coordinates <- do.call(rbind, st_geometry(data_predictions(CPUE_starve)@locations))
prediction_df <- cbind(prediction_df, coordinates)
colnames(prediction_df)[which(colnames(prediction_df)=="1"|colnames(prediction_df)=="2")] <- c("X","Y") 
#simulate 250 times to get DHARMa residuals
sim<-do.call(
  cbind,
  mclapply(
    seq(250),
    function(i) {
      sim<-strv_simulate(
        CPUE_starve,
        conditional=TRUE
      )
      return(dat(sim)$std.WGT_log)
    },
    mc.cores=1
  )
)
#get and plot the DHARMa residuals
CPUE_residuals<-createDHARMa(
  simulated=sim,
  observed=dat(CPUE_starve)$std.WGT_log,
  fittedPredictedResponse = prediction_df$response,
  integer=FALSE
)
plot(CPUE_residuals) #Figure S4

#geographically plot the residuals for the CPUE sub-model
ggplot() + geom_point(size = 0.9, aes(x = X, y = Y, col = CPUE_residuals$scaledResiduals),
                           data = prediction_df) + 
  scale_colour_gradientn(colours = viridis(100)) +
  labs(col = "Quantile\nResidual\n") + theme(text=element_text(size=16)) + 
  geom_polygon(data = fortify(map), aes(group = group, x = long, y = lat)) +
  xlab("Eastings") + ylab("Northings") + facet_wrap(~year) + 
  coord_cartesian(ylim = c(4650, 5250), xlim = c(125, 1000)) + 
  geom_polygon(data = fortify(map1), aes(group = group, x = long, y = lat))
#temporally plot the residuals for the CPUE sub-model
ggplot(aes(x = year, y = CPUE_residuals$scaledResiduals), data = prediction_df) + 
  geom_point(size = 0.15) + geom_smooth(method = "loess") + ylab("Quantile Residual") + xlab("Year") + 
  theme(text=element_text(size=16))

################################################PREDICTIONS########################################

#load in raster data...
raster_data_df <- read.csv(paste(getwd(),"/Data/raster_data_df_timeseries.csv", sep = ""))

#transform some covariates
raster_data_df$RangeStressLog <- log(raster_data_df$RangeStress)
raster_data_df$BtmStressBNAMLog <- log(raster_data_df$BtmStressBNAM)
raster_data_df$sqrt_DEM_Slope <- sqrt(raster_data_df$DEM_Slope)
raster_data_df$DEM_log <- log(-raster_data_df$DEM)
#scale covariates using training data
raster_data_df[paste(covariates,"Scaled", sep="")] <- mapply(x = raster_data_df[covariates], 
                                                             y = st_set_geometry(survey_data_combined[covariates], 
                                                                                 NULL),
                                                       function(x, y) scale(x, center = mean(y), 
                                                                            scale = sd(y)))

#convert the data frame into starve-friendly format
raster_data_df <- st_as_sf(raster_data_df, coords = c("x", "y"), crs = "+proj=utm +zone=20T +units=km")
#create splits to divide up raster_data... splitting up may help with memory issues in extreme cases
split_num <- 1 #sets number of splits to use... if this is set to 1, no splitting occurs
a <- 1:nrow(raster_data_df)
splits <- split(a, ceiling(seq_along(a)/(length(a)/split_num)))

#make the encounter predictions (by going through splits and combining each result)
encounter_prediction_df <- data.frame()
for (i in 1:length(splits))
{
  print(i)
  predictions_encounter <- strv_predict(encounter_starve, raster_data_df[splits[[i]],])
  new_prediction_df<- cbind(do.call(cbind,(predictions_encounter@values[,,1])),
                                 locations(predictions_encounter))
  encounter_prediction_df <- rbind(encounter_prediction_df, new_prediction_df)
}
colnames(encounter_prediction_df)[1:6]<- c("w","w_se","linear","linear_se","response","response_se")

#make the CPUE predictions (by going through splits and combining each result)
CPUE_prediction_df <- data.frame()
for (i in 1:length(splits))
{
  print(i)
  predictions_CPUE <- strv_predict(CPUE_starve, raster_data_df[splits[[i]],])
  new_prediction_df<- cbind(do.call(cbind,(predictions_CPUE@values[,,1])),
                            locations(predictions_CPUE))
  CPUE_prediction_df <- rbind(CPUE_prediction_df, new_prediction_df)
}
colnames(CPUE_prediction_df)[1:6]<- c("w","w_se","linear","linear_se","response","response_se")

#create final combined predictions and add to encounter_prediction_df
encounter_prediction_df$combined_predict <- log(exp(CPUE_prediction_df$linear+
                                                     0.5*response_parameters(CPUE_starve)$std.WGT_log$par^2)*
                                                 encounter_prediction_df$response)
#get combined standard errors
encounter_prediction_df$se_combined <- sqrt(encounter_prediction_df$linear_se^2*(1-encounter_prediction_df$response)^2+
                                        CPUE_prediction_df$linear_se^2)
#create cutoff for low predictions
encounter_prediction_df$NAs <- ifelse(encounter_prediction_df$combined_predict<(-8), "< -8", NA)

#plot random fields for encounter model (Figure S20)
ggplot() + geom_raster(data = encounter_prediction_df, aes(x = UTMX, y = UTMY, fill = w)) +
  coord_equal() + scale_fill_gradientn(name = "", colours = viridis(100)) + 
  geom_polygon(data = fortify(map), aes(group = group, x = long, y = lat)) +
  theme(text=element_text(size=16), axis.text = element_text(color="black")) + 
  coord_cartesian(xlim = c(75, 1050), ylim = c(4600, 5250)) +
  geom_polygon(data = fortify(map1), aes(group = group, x = long, y = lat)) + 
  ylab("Northings") + xlab("Eastings") +
  ggtitle("Encounter Model Random Effects") + facet_wrap(~year, ncol = 4)

#plot random fields for conditional CPUE (Figure S21)
ggplot() + geom_raster(data = CPUE_prediction_df, aes(x = UTMX, y = UTMY, fill = w)) +
  coord_equal() + scale_fill_gradientn(name = "", colours = viridis(100)) + 
  geom_polygon(data = fortify(map), aes(group = group, x = long, y = lat)) +
  theme(text=element_text(size=16), axis.text = element_text(color="black")) + 
  coord_cartesian(xlim = c(75, 1050), ylim = c(4600, 5250)) +
  geom_polygon(data = fortify(map1), aes(group = group, x = long, y = lat)) + 
  ylab("Northings") + xlab("Eastings") +
  ggtitle("CPUE Model Random Effects") + facet_wrap(~year, ncol = 4)

#plot the encounter predictions (Figure S12)
ggplot() + geom_raster(data = encounter_prediction_df, aes(x = UTMX, 
                      y = UTMY, fill = response)) +
  coord_equal() + scale_fill_gradientn(name = "", colours = viridis(100), limits = c(0,1)) + 
  geom_polygon(data = fortify(map), aes(group = group, x = long, y = lat)) +
  theme(text=element_text(size=16), axis.text = element_text(color="black")) + 
  coord_cartesian(xlim = c(75, 1050), ylim = c(4600, 5250)) +
  geom_polygon(data = fortify(map1), aes(group = group, x = long, y = lat)) + 
  ylab("Northings") + xlab("Eastings") +
  ggtitle("Predicted Probability of an Encounter") + facet_wrap(~year, ncol = 4)

#plot the conditional CPUE predictions (Figure S13)
ggplot() + geom_raster(data = CPUE_prediction_df, 
                       aes(x = UTMX, y = UTMY, fill = linear)) +
  coord_equal() + scale_fill_gradientn(name = "", colours = viridis(100), limits = c(-6,11.5)) + 
  geom_polygon(data = fortify(map), aes(group = group, x = long, y = lat)) +
  theme(text=element_text(size=16), axis.text = element_text(color="black")) + 
  coord_cartesian(xlim = c(75, 1050), ylim = c(4600, 5250)) +
  geom_polygon(data = fortify(map1), aes(group = group, x = long, y = lat)) + 
  ylab("Northings") + xlab("Eastings") +
  ggtitle(expression(paste("Predicted Log(CPUE [kg/",km^{2},"]) [Given an Encounter]"))) + 
  facet_wrap(~year, ncol = 4)

#plot the unconditional CPUE predictions (Figure S10)
ggplot() + geom_raster(data = encounter_prediction_df, 
                       aes(x = UTMX, y = UTMY, fill = combined_predict)) +
  coord_equal() + scale_fill_gradientn(name = "", colours = viridis(100), limits = c(-8,12), 
                                       breaks = c(-8,-4,0,4,8,12)) + facet_wrap(~year)  + 
  geom_polygon(data = fortify(map), aes(group = group, x = long, y = lat)) +
  theme(text=element_text(size=16), axis.text = element_text(color="black")) + 
  coord_cartesian(xlim = c(75, 1050), ylim = c(4600, 5250)) +
  geom_polygon(data = fortify(map1), aes(group = group, x = long, y = lat)) + 
  ylab("Northings") + xlab("Eastings") +
  ggtitle(expression(paste("Log(Predicted CPUE [kg/",km^{2},"])"))) + facet_wrap(~year, ncol = 4) + 
  new_scale_fill() + scale_fill_manual(name=NULL, values="black") + 
  geom_raster(data = na.omit(encounter_prediction_df[c("UTMX","UTMY","NAs")]), 
              aes(x = UTMX, y = UTMY, fill = as.factor(NAs)))

#get uncertainty in final predictions (Figure S11)
ggplot() + geom_raster(data = encounter_prediction_df, 
                       aes(x = UTMX, y = UTMY, fill = se_combined)) +
  coord_equal() + scale_fill_gradientn(name = "", colours = viridis(100), limits = c(0,6))  + 
  geom_polygon(data = fortify(map), aes(group = group, x = long, y = lat)) +
  theme(text=element_text(size=16), axis.text = element_text(color="black")) + 
  coord_cartesian(xlim = c(75, 1050), ylim = c(4600, 5250)) +
  geom_polygon(data = fortify(map1), aes(group = group, x = long, y = lat)) + 
  ylab("Northings") + xlab("Eastings") +
  ggtitle(expression(paste("Log(Predicted CPUE [kg/",km^{2},"]) - Std. Errors"))) + facet_wrap(~year, ncol = 4)

#remove starve geometry
encounter_prediction_df <- st_set_geometry(encounter_prediction_df, NULL)
#export prediction and standard error rasters for external use (i.e., in ArcGIS... Figures 4 and 6)
starve_prediction_raster <- raster_create(encounter_prediction_df, "combined_predict", year_start = 2000)
crs(starve_prediction_raster) <- "+datum=WGS84 +proj=utm +zone=20T +units=km" 
writeRaster(starve_prediction_raster, paste(getwd(),"/Model Output/starve/starve_predictions.tif", sep = ""),
            format = "GTiff", bylayer = T, suffix = seq(2000,2019, by = 1), overwrite = T)
se_raster <- raster_create(encounter_prediction_df, "se_combined", year_start = 2000)
crs(se_raster) <- "+datum=WGS84 +proj=utm +zone=20T +units=km" 
writeRaster(se_raster, paste(getwd(),"/Model Output/starve/starve_standard_errors.tif", sep = ""),
            format = "GTiff", bylayer = T, suffix = seq(2000,2019, by = 1), overwrite = T)

#save model and raster predictions as .RData
save(encounter_starve, CPUE_starve, starve_prediction_raster, file = "starve.RData")