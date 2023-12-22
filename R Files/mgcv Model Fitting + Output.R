#load necessary libraries
library(mgcv); library(ggplot2); library(knitr)
library(dplyr); library(sp); library(reshape)
library(gtools); library(raster); library(DHARMa)
library(gratia); library(viridis); library(ggnewscale)
source("R Files/Helper Functions.R")

#increase memory limit
memory.limit(size = 400000)

#set number of decimals to round tables to
round_val <- 3

#load in pre-processed survey data
survey_data_combined <- read.csv(paste(getwd(),"/Data/survey_data_combined.csv", sep = ""))
#transform some covariates
covariates <- c("BtmTempBNAM","BtmSalinityBNAM","BtmStressBNAMLog","RangeTemp",
                "RangeStressLog","DEM_log","sqrt_DEM_Slope","DEM_Easterness",
                "DEM_Northerness","DEM_RDMV")
survey_data_combined$RangeStressLog <- log(survey_data_combined$RangeStress)
survey_data_combined$BtmStressBNAMLog <- log(survey_data_combined$BtmStressBNAM)
survey_data_combined$sqrt_DEM_Slope <- sqrt(survey_data_combined$DEM_Slope)
survey_data_combined$DEM_log <- log(-survey_data_combined$DEM)
#scale covariates using all data
survey_data_combined[paste(covariates,"Scaled", sep="")] <- scale(survey_data_combined[covariates])

#grab NS, New Brunswick, Quebec, PEI, NFLD maps for later plots
map <- raster::getData(country = "CAN", level = 1)
map <- map[which(map$NAME_1 %in% c("Nova Scotia", "Prince Edward Island", 
                                   "New Brunswick", "Newfoundland and Labrador", "Québec")),]
map <- spTransform(map, crs("+datum=WGS84 +proj=utm +zone=20T +units=km"))
#grab Maine map for later plots
map1 <- raster::getData(country = "USA", level = 1)
map1 <- map1[which(map1$NAME_1 == "Maine"),]
map1 <- spTransform(map1, crs("+datum=WGS84 +proj=utm +zone=20T +units=km"))

#remove RV winter survey tows
#survey_data_combined <- survey_data_combined[-which(survey_data_combined$season_winter == 1),]

#load in spatial domain shape file
shp <- shapefile("Data/SpatialDomain/MaritimesRegionEcosystemAssessmentStrata_SSsubset_ForSeaCuke.shp")
shp <- spTransform(shp, crs("+datum=WGS84 +proj=utm +zone=20T +units=km"))

##########################################ENCOUNTER SUB-MODEL###########################################

#fit the sub-model on all data
start_time <- Sys.time()
encounter_mgcv <- gam(presence~te(UTMX, UTMY, year, k = c(400, 5), 
                                  bs = c("ds","tp"), m = list(c(1,0.5),NA), d = c(2,1)) + 
                        s(DEM_logScaled, k = 3)+s(BtmTempBNAMScaled, k = 3) + s(RangeTempScaled, k = 3) + 
                        s(BtmSalinityBNAMScaled, k = 3) + s(BtmStressBNAMLogScaled, k = 3) + 
                        s(RangeStressLogScaled, k = 3) + s(sqrt_DEM_SlopeScaled, k = 3) + 
                        s(DEM_NorthernessScaled, k = 3) + s(DEM_EasternessScaled, k = 3) + 
                        s(DEM_RDMVScaled, k = 3) + as.numeric(snowcrab), 
                      family = "binomial", 
                      data = survey_data_combined, 
                      method = "REML") 
end_time <- Sys.time()
fit_time <- end_time - start_time #get time to fit

#plot the spatio-temporal smooth (Figure S18)
draw(encounter_mgcv, select = "te(UTMX,UTMY,year)", n_3d = 20, contour = F)+
  geom_polygon(data = fortify(shp), aes(group = group, x = long, y = lat), fill = NA, 
               col = "white") +
  geom_polygon(data = fortify(map), aes(group = group, x = long, y = lat)) +
  theme(text=element_text(size=16), axis.text=element_text(color="black")) +
  geom_polygon(data = fortify(map1), aes(group = group, x = long, y = lat)) + 
  coord_cartesian(ylim = c(4700, 5200), xlim = c(175, 925)) + 
  scale_fill_gradientn(colours = viridis(100)) + ggtitle("f(Eastings, Northings, Year)") + 
  xlab("Eastings") + ylab("Northings") + labs(fill = "Partial\nEffect")

#ready information to put into gam summary table for encounter sub-model
a <- data.frame(format(round(summary(encounter_mgcv)$p.table, round_val), digits = round_val))
colnames(a) <- c("Estimate", "Std. Error", "Z","P-value")
rownames(a) <- c("Intercept", "Snow Crab Survey")
b <- data.frame(format(round(summary(encounter_mgcv)$s.table, round_val), digits = round_val))
colnames(b) <- c("EDF","Ref. DF","Chi-sq.","P-value")
rownames(b) <- c("f(Eastings, Northings, Year)", "f(Log(Depth))", "f(Bottom Temperature)", 
                 "f(Btm. Temperature Range)", "f(Bottom Salinity)", "f(Log(Bottom Stress))", 
                 "f(Log(Btm. Stress Range))", "f(Sqrt(Slope))", "f(Northerness)","f(Easterness)", 
                 "f(RDMV)")
#make Table S1
c <- kable(list(a,b), format = "latex", caption = "Encounter Model", align = "rc")

#compute and plot the DHARMa residuals, use default of 250 simulations
survey_data_combined["residuals_encounter"] <- simulateResiduals(fittedModel = encounter_mgcv, 
                                                                plot = T)$scaledResiduals

#geographically plot the residuals for the encounter sub-model
ggplot() + geom_point(size = 0.9, aes(x = UTMX, y = UTMY, col = residuals_encounter),
                           data = survey_data_combined) + 
  scale_colour_gradientn(colours = viridis(100)) +
  labs(col = "Quantile\nResidual\n") + theme(text=element_text(size=16)) + 
  geom_polygon(data = fortify(map), aes(group = group, x = long, y = lat)) +
  xlab("Eastings") + ylab("Northings") + coord_cartesian(ylim = c(4650, 5250), xlim = c(125, 1000)) + 
  geom_polygon(data = fortify(map1), aes(group = group, x = long, y = lat))
#temporally plot the residuals for the encounter sub-model
ggplot(aes(x = year, y = residuals_encounter), data = survey_data_combined) + 
  geom_point(size = 0.15) + geom_smooth(method = "loess") + ylab("Quantile Residual") + xlab("Year") + 
  theme(text=element_text(size=16))

#################################################CPUE SUB-MODEL###########################################

#grab encounter tows and set up response
survey_data_combined_CPUE <- survey_data_combined[which(is.na(survey_data_combined$std.WGT)==FALSE&
                                                          survey_data_combined$std.WGT!=0),]
survey_data_combined_CPUE$std.WGT_log <- log(survey_data_combined_CPUE$std.WGT)

#fit the sub-model with all the data
start_time <- Sys.time()
CPUE_mgcv <- gam(log(std.WGT) ~ te(UTMX, UTMY, year, k = c(300, 5), 
                                      bs = c("ds","tp"), m = list(c(1,0.5),NA), d = c(2,1)) + 
                      s(DEM_logScaled, k = 3)+s(BtmTempBNAMScaled, k = 3) + s(RangeTempScaled, k = 3) + 
                      s(BtmSalinityBNAMScaled, k = 3) + s(BtmStressBNAMLogScaled, k = 3) + 
                      s(RangeStressLogScaled, k = 3) + s(sqrt_DEM_SlopeScaled, k = 3) + 
                      s(DEM_NorthernessScaled, k = 3) + s(DEM_EasternessScaled, k = 3) + 
                      s(DEM_RDMVScaled, k = 3) + as.numeric(snowcrab), 
                    family = "gaussian", 
                    data = survey_data_combined_CPUE, 
                    method = "REML") 
end_time <- Sys.time()
fit_time3 <- end_time - start_time #get time to fit

#plot the spatio-temporal smooth (Figure S19)
draw(CPUE_mgcv, select = "te(UTMX,UTMY,year)", n_3d = 20, contour = F)+
  geom_polygon(data = fortify(shp), aes(group = group, x = long, y = lat), fill = NA, col = "white") +
  geom_polygon(data = fortify(map), aes(group = group, x = long, y = lat)) +
  theme(text=element_text(size=16), axis.text=element_text(color="black")) +
  geom_polygon(data = fortify(map1), aes(group = group, x = long, y = lat)) + 
  coord_cartesian(ylim = c(4750, 5150), xlim = c(250, 900)) + 
  scale_fill_gradientn(colours = viridis(100)) + ggtitle("f(Eastings, Northings, Year)") +
  xlab("Eastings") + ylab("Northings") + labs(fill = "Partial\nEffect")

#ready information to put into gam summary table for CPUE sub-model
a <- data.frame(format(round(summary(CPUE_mgcv)$p.table, round_val), digits = round_val))
colnames(a) <- c("Estimate", "Std. Error", "Z","P-value")
rownames(a) <- c("Intercept", "Snow Crab Survey")
b <- data.frame(format(round(summary(CPUE_mgcv)$s.table, round_val), digits = round_val))
colnames(b) <- c("EDF","Ref. DF","F","P-value")
rownames(b) <- c("f(Eastings, Northings, Year)", "f(Log(Depth))", "f(Bottom Temperature)", 
                 "f(Btm. Temperature Range)", "f(Bottom Salinity)", "f(Log(Bottom Stress))", 
                 "f(Log(Btm. Stress Range))", "f(Sqrt(Slope))", "f(Northerness)","f(Easterness)", "f(RDMV)")
#make Table S2
c <- kable(list(a,b), format = "latex", caption = "CPUE Model", align = "rc")

#compute and plot the DHARMa residuals, use default of 250 simulations... Figure S3
survey_data_combined_CPUE["residuals_CPUE"] <- simulateResiduals(fittedModel = CPUE_mgcv, 
                                                                 plot = T)$scaledResiduals

#geographically plot the residuals for the CPUE sub-model
ggplot() + geom_point(size = 0.9, aes(x = UTMX, y = UTMY, col = residuals_CPUE),
                           data = survey_data_combined_CPUE) + 
  scale_colour_gradientn(colours = viridis(100)) +
  labs(col = "Quantile\nResidual") + theme(text=element_text(size=20)) + 
  geom_polygon(data = fortify(map), aes(group = group, x = long, y = lat)) +
  xlab("Eastings") + ylab("Northings") + 
  coord_cartesian(ylim = c(4650, 5250), xlim = c(75, 1050)) + 
  geom_polygon(data = fortify(map1), aes(group = group, x = long, y = lat))
#temporally plot the residuals for the CPUE sub-model
ggplot(aes(x = year, y = residuals_CPUE), data = survey_data_combined_CPUE) + 
  geom_point(size = 0.15) + geom_smooth(method = "loess") + ylab("Quantile Residual") + xlab("Year") + 
  theme(text=element_text(size=16))

############################################PREDICTIONS##################################################

#load in raster data...
raster_data_df <- read.csv(paste(getwd(),"/Data/raster_data_df_timeseries.csv", sep = ""))

#transform some covariates
raster_data_df$RangeStressLog <- log(raster_data_df$RangeStress)
raster_data_df$BtmStressBNAMLog <- log(raster_data_df$BtmStressBNAM)
raster_data_df$sqrt_DEM_Slope <- sqrt(raster_data_df$DEM_Slope)
raster_data_df$DEM_log <- log(-raster_data_df$DEM)
#scale covariates using all training data
raster_data_df[paste(covariates,"Scaled", sep="")] <- mapply(x = raster_data_df[covariates], 
                                                             y = survey_data_combined[covariates], 
                                                       function(x, y) scale(x, center = mean(y), 
                                                                            scale = sd(y)))

#get predictions...
stack_predict_CPUE <- predict(CPUE_mgcv, newdata = raster_data_df, se.fit = TRUE)
stack_predict_zero <- predict(encounter_mgcv, newdata = raster_data_df, se.fit = TRUE)
stack_predict_zero$fit <- 1/(1+exp(-stack_predict_zero$fit)) #convert log odds to probability
#combine to get final abundance predictions
stack_predict_combined <- stack_predict_zero
stack_predict_combined$fit <- log(stack_predict_zero$fit*exp(stack_predict_CPUE$fit+0.5*CPUE_mgcv$scale))
#get combined standard errors
stack_predict_combined$se.fit <- sqrt(stack_predict_zero$se.fit^2*(1-stack_predict_zero$fit)^2+
                                        stack_predict_CPUE$se.fit^2)
#create cutoff for low predictions
raster_data_df$NAs <- ifelse(stack_predict_combined$fit<(-8), "< -8", NA)

#plot the encounter probability predictions (Figure S8)
ggplot() + geom_raster(data = raster_data_df, aes(x = UTMX, y = UTMY, fill = stack_predict_zero$fit)) +
  coord_equal() + scale_fill_gradientn(name = "", colours = viridis(100), limits = c(0,1)) + 
  geom_polygon(data = fortify(map), aes(group = group, x = long, y = lat)) +
  theme(text=element_text(size=16), axis.text=element_text(color="black")) + 
  coord_cartesian(xlim = c(75, 1050), ylim = c(4600, 5250)) +
  geom_polygon(data = fortify(map1), aes(group = group, x = long, y = lat)) + 
  ylab("Northings") + xlab("Eastings") +
  ggtitle("Predicted Probability of an Encounter") + facet_wrap(~year, ncol = 4)

#plot the conditional CPUE predictions (Figure S9)
ggplot() + geom_raster(data = raster_data_df, 
                       aes(x = UTMX, y = UTMY, fill = stack_predict_CPUE$fit)) +
  coord_equal() + scale_fill_gradientn(name = "", colours = viridis(100), limits = c(-6,11.5)) + 
  geom_polygon(data = fortify(map), aes(group = group, x = long, y = lat)) +
  theme(text=element_text(size=16), axis.text=element_text(color="black")) + 
  coord_cartesian(xlim = c(75, 1050), ylim = c(4600, 5250)) +
  geom_polygon(data = fortify(map1), aes(group = group, x = long, y = lat)) + 
  ylab("Northings") + xlab("Eastings") +
  ggtitle(expression(paste("Predicted Log(CPUE [kg/",km^{2},"]) [Given an Encounter]"))) + 
  facet_wrap(~year, ncol = 4)

#plot the unconditional CPUE predictions (Figure S10)
ggplot() + geom_raster(data = raster_data_df, 
                       aes(x = UTMX, y = UTMY, fill = stack_predict_combined$fit)) +
  coord_equal() + scale_fill_gradientn(name = "", colours = viridis(100), limits = c(-8,12), 
                                       breaks = c(-8,-4,0,4,8,12)) + 
  geom_polygon(data = fortify(map), aes(group = group, x = long, y = lat)) +
  theme(text=element_text(size=16), axis.text=element_text(color="black")) + 
  coord_cartesian(xlim = c(75, 1050), ylim = c(4600, 5250)) +
  geom_polygon(data = fortify(map1), aes(group = group, x = long, y = lat)) + 
  ylab("Northings") + xlab("Eastings") +
  ggtitle(expression(paste("Log(Predicted CPUE [kg/",km^{2},"])"))) + facet_wrap(~year, ncol = 4) +
  new_scale_fill() + scale_fill_manual(name=NULL, values="black") + 
  geom_raster(data = na.omit(raster_data_df[c("UTMX","UTMY","NAs")]), 
              aes(x = UTMX, y = UTMY, fill = as.factor(NAs)))

#get uncertainty in final predictions (Figure S11)
ggplot() + geom_raster(data = raster_data_df, 
                       aes(x = UTMX, y = UTMY, fill = stack_predict_combined$se.fit)) +
  coord_equal() + scale_fill_gradientn(name = "", colours = viridis(100), limits = c(0,6)) + 
  geom_polygon(data = fortify(map), aes(group = group, x = long, y = lat)) +
  theme(text=element_text(size=16), axis.text=element_text(color="black")) + 
  coord_cartesian(xlim = c(75, 1050), ylim = c(4600, 5250)) +
  geom_polygon(data = fortify(map1), aes(group = group, x = long, y = lat)) + 
  ylab("Northings") + xlab("Eastings") +
  ggtitle(expression(paste("Log(Predicted CPUE [kg/",km^{2},"]) - Std. Errors"))) + 
  facet_wrap(~year, ncol = 4)

#export prediction and standard error rasters for external use (i.e., in ArcGIS... Figures 4 and 6)
mgcv_prediction_raster <- raster_create_gam(stack_predict_combined$fit, raster_data_df, year_start = 2000)
crs(mgcv_prediction_raster) <- "+datum=WGS84 +proj=utm +zone=20T +units=km" 
writeRaster(mgcv_prediction_raster, paste(getwd(),"/Model Output/mgcv/mgcv_predictions.tif", sep = ""),
            format = "GTiff", bylayer = T, suffix = seq(2000,2019, by = 1), overwrite = T)
se_raster <- raster_create_gam(stack_predict_combined$se.fit, raster_data_df, year_start = 2000)
crs(se_raster) <- "+datum=WGS84 +proj=utm +zone=20T +units=km" 
writeRaster(se_raster, paste(getwd(),"/Model Output/mgcv/mgcv_standard_errors.tif", sep = ""),
            format = "GTiff", bylayer = T, suffix = seq(2000,2019, by = 1), overwrite = T)

#save model and raster predictions as .RData
save(encounter_mgcv, CPUE_mgcv, mgcv_prediction_raster, file = "mgcv.RData")