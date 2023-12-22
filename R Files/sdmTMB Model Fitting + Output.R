#load necessary libraries
library(sdmTMB); library(INLA); library(raster)
library(knitr); library(sp); library(ggplot2)
library(inlabru); library(viridis); library(ggnewscale)
library(DHARMa); source("R Files/Helper Functions.R")

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

#grab NS, Quebec, New Brunswick, PEI, NFLD maps for later plots
map <- raster::getData(country = "CAN", level = 1)
map <- map[which(map$NAME_1 %in% c("Nova Scotia", "Prince Edward Island", 
                                   "New Brunswick", "Newfoundland and Labrador", "Québec")),]
map <- spTransform(map, crs("+datum=WGS84 +proj=utm +zone=20T +units=km"))
#grab Maine map for later plots
map1 <- raster::getData(country = "USA", level = 1)
map1 <- map1[which(map1$NAME_1 == "Maine"),]
map1 <- spTransform(map1, crs("+datum=WGS84 +proj=utm +zone=20T +units=km"))

#reorder the data by year
survey_data_combined <- survey_data_combined[order(survey_data_combined$year, decreasing = FALSE),]

#remove RV winter survey tows
#survey_data_combined <- survey_data_combined[-which(survey_data_combined$season_winter == 1),]

#set up the mesh
non_convex_bdry <- inla.nonconvex.hull(cbind(survey_data_combined$UTMX, survey_data_combined$UTMY), 
                                       -0.03, resolution = c(200, 50))
mesh <- inla.mesh.2d(boundary = non_convex_bdry, max.edge=c(25, 200), 
                     offset = c(100), cutoff = 11.5, 
                     loc = cbind(survey_data_combined$UTMX, survey_data_combined$UTMY))
#plot the mesh (Figure S1)
ggplot() + 
  geom_point(size = 0.9, aes(x = UTMX, y = UTMY), col = "orange", data = survey_data_combined) + 
  theme(text=element_text(size=24), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.text=element_text(color="black")) + 
  geom_polygon(data = fortify(map), aes(group = group, x = long, y = lat)) +
  xlab("Eastings") + ylab("Northings") + 
  coord_cartesian(ylim = c(4500, 5400), xlim = c(0, 1105)) + 
  geom_polygon(data = fortify(map1), aes(group = group, x = long, y = lat)) + 
  gg(mesh, int.linewidth = 0.8, edge.color = "grey44")

#########################################ENCOUNTER SUB-MODEL#########################################

#fit the encounter sub-model in sdmTMB using the mesh... use random walk
mesh_model <- make_mesh(survey_data_combined, xy_cols = c("UTMX", "UTMY"), mesh = mesh)
start_time <- Sys.time()
encounter_sdmTMB <- sdmTMB(presence ~ DEM_logScaled+I(DEM_logScaled^2)+BtmTempBNAMScaled+
                           I(BtmTempBNAMScaled^2)+RangeTempScaled+I(RangeTempScaled^2)+
                           BtmSalinityBNAMScaled+I(BtmSalinityBNAMScaled^2)+BtmStressBNAMLogScaled+
                           RangeStressLogScaled+I(RangeStressLogScaled^2)+sqrt_DEM_SlopeScaled+
                           I(sqrt_DEM_SlopeScaled^2)+DEM_NorthernessScaled+
                           I(DEM_NorthernessScaled^2)+DEM_EasternessScaled+DEM_RDMVScaled+
                           I(DEM_RDMVScaled^2)+as.numeric(snowcrab),
  data = survey_data_combined,
  mesh = mesh_model,
  family = binomial(),
  spatiotemporal = "rw", spatial = "on", share_range = TRUE, time = "year", 
  control = sdmTMBcontrol(newton_loops = 1, nlminb_loops = 2))
end_time <- Sys.time()
fit_time1 <- end_time - start_time #get time to fit
#sanity check
sanity(encounter_sdmTMB)

#create table to describe fixed effects, etc. - Table S5
fixed_effects_table <- round(tidy(encounter_sdmTMB)[2:3], round_val)
rownames(fixed_effects_table) <- c("Intercept","Log(Depth)","Log(Depth) Squared",
                                   "Bottom Temperature", "Bottom Temperature Squared",
                                   "Btm. Temperature Range","Btm. Temperature Range Squared",
                                   "Bottom Salinity", "Bottom Salinity Squared", 
                                   "Log(Bottom Stress)", "Log(Btm. Stress Range)", 
                                   "Log(Btm. Stress Range) Squared",
                                   "Sqrt(Slope)","Sqrt(Slope) Squared", "Northerness",
                                   "Northerness Squared","Easterness","RDMV","RDMV Squared",
                                   "Snow Crab Survey")
colnames(fixed_effects_table) <- c("Estimate", "Std. Error")
spatiotemporal_parameters_table <- round(tidy(encounter_sdmTMB, effects = "ran_pars", 
                                              conf.int = FALSE)[2], round_val)
se <- as.list(encounter_sdmTMB$sd_report, "Std. Error", report = TRUE) #grab SEs in natural space
spatiotemporal_parameters_table <- cbind(spatiotemporal_parameters_table, 
                                         round(c(se$range[1], se$sigma_O, se$sigma_E[1]), round_val))
colnames(spatiotemporal_parameters_table) <- c("Estimate", "Std. Error")
rownames(spatiotemporal_parameters_table) <- c("Range", "Marginal Spatial SD","Marginal Spatio-Temporal SD")
kable(list(fixed_effects_table, spatiotemporal_parameters_table), format = "latex")

#get DHARMa residuals, using 250 simulations
sim <- simulate(encounter_sdmTMB, nsim = 250, seed = 2)
residuals_encounter <- DHARMa::createDHARMa(
  simulatedResponse = sim,
  observedResponse = as.numeric(survey_data_combined$presence),
  integer = T,
  fittedPredictedResponse = predict(encounter_sdmTMB, type = "response")$est
)
plot(residuals_encounter)

#spatio-temporally plot residuals
ggplot() + geom_point(size = 0.9, aes(x = UTMX, y = UTMY, col = residuals_encounter$scaledResiduals),
                           data = survey_data_combined) + 
  scale_colour_gradientn(colours = viridis(100)) +
  labs(col = "Quantile\nResidual\n") + theme(text=element_text(size=16)) + 
  geom_polygon(data = fortify(map), aes(group = group, x = long, y = lat)) +
  xlab("Eastings") + ylab("Northings") + facet_wrap(~year) + 
  coord_cartesian(ylim = c(4650, 5250), xlim = c(125, 1000)) + 
  geom_polygon(data = fortify(map1), aes(group = group, x = long, y = lat))
#temporally plot residuals
ggplot(aes(x = year, y = residuals_encounter$scaledResiduals), data = survey_data_combined) + 
  geom_point(size = 0.15) + geom_smooth(method = "loess") + ylab("Quantile Residual") + xlab("Year") + 
  theme(text=element_text(size=16))

#########################################CPUE SUB-MODEL###########################################

#grab encounter tows and set up the response
survey_data_combined_CPUE <- survey_data_combined[which(is.na(survey_data_combined$std.WGT)==FALSE&
                                                          survey_data_combined$std.WGT!=0),]
survey_data_combined_CPUE$std.WGT_log <- log(survey_data_combined_CPUE$std.WGT)

#fit the CPUE sub-model in sdmTMB using the mesh... use random walk
mesh_model <- make_mesh(survey_data_combined_CPUE, xy_cols = c("UTMX", "UTMY"), mesh = mesh)
start_time <- Sys.time()
CPUE_sdmTMB <- sdmTMB(std.WGT_log ~ DEM_logScaled+I(DEM_logScaled^2)+BtmTempBNAMScaled+
                       I(BtmTempBNAMScaled^2)+RangeTempScaled+
                       I(RangeTempScaled^2)+BtmSalinityBNAMScaled+I(BtmSalinityBNAMScaled^2)+
                       BtmStressBNAMLogScaled+RangeStressLogScaled+I(RangeStressLogScaled^2)+
                       sqrt_DEM_SlopeScaled+I(sqrt_DEM_SlopeScaled^2)+DEM_NorthernessScaled+
                       I(DEM_NorthernessScaled^2)+DEM_EasternessScaled+DEM_RDMVScaled+
                       I(DEM_RDMVScaled^2)+as.numeric(snowcrab),
                     data = survey_data_combined_CPUE,
                     mesh = mesh_model,
                     family = gaussian(),
                     spatiotemporal = "rw", spatial = "on", share_range = TRUE, time = "year", 
                     control = sdmTMBcontrol(newton_loops = 1))
end_time <- Sys.time()
fit_time3 <- end_time - start_time #get time to fit
#sanity check
sanity(CPUE_sdmTMB)

#create table to describe fixed effects, etc. - Table S6
fixed_effects_table <- round(tidy(CPUE_sdmTMB)[2:3], round_val)
rownames(fixed_effects_table) <- c("Intercept","Log(Depth)","Log(Depth) Squared",
                                   "Bottom Temperature", "Bottom Temperature Squared",
                                   "Btm. Temperature Range","Btm. Temperature Range Squared",
                                   "Bottom Salinity", "Bottom Salinity Squared", 
                                   "Log(Bottom Stress)", "Log(Btm. Stress Range)", 
                                   "Log(Btm. Stress Range) Squared",
                                   "Sqrt(Slope)","Sqrt(Slope) Squared", "Northerness",
                                   "Northerness Squared","Easterness","RDMV","RDMV Squared",
                                   "Snow Crab Survey")
colnames(fixed_effects_table) <- c("Estimate", "Std. Error")
spatiotemporal_parameters_table <- round(tidy(CPUE_sdmTMB, effects = "ran_pars", 
                                              conf.int = FALSE)[-2,2], round_val)
se <- as.list(CPUE_sdmTMB$sd_report, "Std. Error", report = TRUE) #grab SEs in natural space
spatiotemporal_parameters_table <- cbind(spatiotemporal_parameters_table, 
                                         round(c(se$range[1], se$sigma_O, se$sigma_E[1]), round_val))
colnames(spatiotemporal_parameters_table) <- c("Estimate", "Std. Error")
rownames(spatiotemporal_parameters_table) <- c("Range", "Marginal Spatial SD","Marginal Spatio-Temporal SD")
kable(list(fixed_effects_table, spatiotemporal_parameters_table), format = "latex")

#get DHARMa residuals, using 250 simulations
sim <- simulate(CPUE_sdmTMB, nsim = 250, seed = 2)
residuals_CPUE <- DHARMa::createDHARMa(
  simulatedResponse = sim,
  observedResponse = survey_data_combined_CPUE$std.WGT_log,
  integer = F,
  fittedPredictedResponse = predict(CPUE_sdmTMB)$est
)
plot(residuals_CPUE) #Figure S5

#spatio-temporally plot residuals
ggplot() + geom_point(size = 0.9, aes(x = UTMX, y = UTMY, col = residuals_CPUE$scaledResiduals),
                           data = survey_data_combined_CPUE) + 
  scale_colour_gradientn(colours = viridis(100), limits = c(-5,5)) +
  labs(col = "Quantile\nResidual\n") + theme(text=element_text(size=16)) + 
  geom_polygon(data = fortify(map), aes(group = group, x = long, y = lat)) +
  xlab("Eastings") + ylab("Northings") + facet_wrap(~year) + 
  coord_cartesian(ylim = c(4650, 5250), xlim = c(125, 1000)) + 
  geom_polygon(data = fortify(map1), aes(group = group, x = long, y = lat))
#temporally plot residuals
ggplot(aes(x = year, y = residuals_CPUE$scaledResiduals), data = survey_data_combined_CPUE) + 
  geom_point(size = 0.15) + geom_smooth(method = "loess") + ylab("Quantile Residual") + xlab("Year") + 
  theme(text=element_text(size=16))

###########################################PREDICTIONS##########################################

#load in raster data...
raster_data_df <- read.csv(paste(getwd(),"/Data/raster_data_df_timeseries.csv", sep = ""))

#transform some covariates
raster_data_df$RangeStressLog <- log(raster_data_df$RangeStress)
raster_data_df$BtmStressBNAMLog <- log(raster_data_df$BtmStressBNAM)
raster_data_df$sqrt_DEM_Slope <- sqrt(raster_data_df$DEM_Slope)
raster_data_df$DEM_log <- log(-raster_data_df$DEM)
#scale covariates using training data
raster_data_df[paste(covariates,"Scaled", sep="")] <- mapply(x = raster_data_df[covariates], 
                                                             y = survey_data_combined[covariates], 
                                                       function(x, y) scale(x, center = mean(y), 
                                                                            scale = sd(y)))

#make predictions
encounter_prediction_df <- predict(encounter_sdmTMB, newdata = raster_data_df, type = "response")
CPUE_prediction_df <- predict(CPUE_sdmTMB, newdata = raster_data_df)
#get standard errors
encounter_se <- predict(encounter_sdmTMB, newdata = raster_data_df, nsim = 500)
CPUE_se <- predict(CPUE_sdmTMB, newdata = raster_data_df, nsim = 500)
if(nrow(encounter_se)>nrow(raster_data_df))#remove possible extra rows... sdmTMB quirk
{
  encounter_se <- encounter_se[-(nrow(raster_data_df)+1:nrow(encounter_se)),]
  CPUE_se <- CPUE_se[-(nrow(raster_data_df)+1:nrow(CPUE_se)),]
}
encounter_prediction_df$se <- apply(encounter_se, 1, sd)
CPUE_prediction_df$se <- apply(CPUE_se, 1, sd)
#create final combined predictions, standard errors, and add to encounter_prediction_df
encounter_prediction_df$combined_predict <- log(exp(CPUE_prediction_df$est+
                                                     0.5*CPUE_sdmTMB$sd_report$value[2]^2)*
                                                 encounter_prediction_df$est)
encounter_prediction_df$combined_se <- sqrt((encounter_prediction_df$se)^2*
                                             (1-encounter_prediction_df$est)^2+(CPUE_prediction_df$se)^2)
#create cutoff for low predictions
encounter_prediction_df$NAs <- ifelse(encounter_prediction_df$combined_predict<(-8), "< -8", NA)

#plot random fields (spatial + spatio-temporal) for encounter model (Figure S22)
ggplot() + geom_raster(data = encounter_prediction_df, aes(x = UTMX, y = UTMY, fill = est_rf)) +
  coord_equal() + scale_fill_gradientn(name = "", colours = viridis(100)) + 
  geom_polygon(data = fortify(map), aes(group = group, x = long, y = lat)) +
  theme(text=element_text(size=16), axis.text=element_text(color="black")) + 
  coord_cartesian(xlim = c(75, 1050), ylim = c(4600, 5250)) +
  geom_polygon(data = fortify(map1), aes(group = group, x = long, y = lat)) + 
  ylab("Northings") + xlab("Eastings") +
  ggtitle("Encounter Model Random Effects (Spatial + Spatio-Temporal)") + facet_wrap(~year, ncol = 4)

#plot random fields (spatial + spatio-temporal) for conditional CPUE (Figure S23)
ggplot() + geom_raster(data = CPUE_prediction_df, aes(x = UTMX, y = UTMY, fill = est_rf)) +
  coord_equal() + scale_fill_gradientn(name = "", colours = viridis(100)) + 
  geom_polygon(data = fortify(map), aes(group = group, x = long, y = lat)) +
  theme(text=element_text(size=16), axis.text=element_text(color="black")) + 
  coord_cartesian(xlim = c(75, 1050), ylim = c(4600, 5250)) +
  geom_polygon(data = fortify(map1), aes(group = group, x = long, y = lat)) + 
  ylab("Northings") + xlab("Eastings") +
  ggtitle("CPUE Model Random Effects (Spatial + Spatio-Temporal)") + facet_wrap(~year, ncol = 4)

#plot the encounter probability predictions (Figure S16)
ggplot() + geom_raster(data = encounter_prediction_df, aes(x = UTMX, y = UTMY, fill = est)) +
  coord_equal() + scale_fill_gradientn(name = "", colours = viridis(100), limits = c(0,1)) + 
  geom_polygon(data = fortify(map), aes(group = group, x = long, y = lat)) +
  theme(text=element_text(size=16), axis.text=element_text(color="black")) + 
  coord_cartesian(xlim = c(75, 1050), ylim = c(4600, 5250)) +
  geom_polygon(data = fortify(map1), aes(group = group, x = long, y = lat)) + 
  ylab("Northings") + xlab("Eastings") +
  ggtitle("Predicted Probability of an Encounter") + facet_wrap(~year, ncol = 4)

#plot the conditional CPUE predictions (Figure S17)
ggplot() + geom_raster(data = CPUE_prediction_df, 
                       aes(x = UTMX, y = UTMY, fill = est)) +
  coord_equal() + scale_fill_gradientn(name = "", colours = viridis(100), limits = c(-6,11.5)) + 
  geom_polygon(data = fortify(map), aes(group = group, x = long, y = lat)) +
  theme(text=element_text(size=16), axis.text=element_text(color="black")) + 
  coord_cartesian(xlim = c(75, 1050), ylim = c(4600, 5250)) +
  geom_polygon(data = fortify(map1), aes(group = group, x = long, y = lat)) + 
  ylab("Northings") + xlab("Eastings") + 
  ggtitle(expression(paste("Log(Predicted CPUE [kg/",km^{2},"]) [Given an Encounter]"))) + 
  facet_wrap(~year, ncol = 4)

#plot the unconditional CPUE predictions (Figure S14)
ggplot() + geom_raster(data = encounter_prediction_df, 
                       aes(x = UTMX, y = UTMY, fill = combined_predict)) +
  coord_equal() + scale_fill_gradientn(name = "", colours = viridis(100), 
                                       limits = c(-8,12), breaks = c(-8,-4,0,4,8,12))  + 
  geom_polygon(data = fortify(map), aes(group = group, x = long, y = lat)) +
  theme(text=element_text(size=16), axis.text=element_text(color="black")) + 
  coord_cartesian(xlim = c(75, 1050), ylim = c(4600, 5250)) +
  geom_polygon(data = fortify(map1), aes(group = group, x = long, y = lat)) + 
  ylab("Northings") + xlab("Eastings") +
  ggtitle(expression(paste("Log(Predicted CPUE [kg/",km^{2},"])"))) + facet_wrap(~year, ncol = 4) +
  new_scale_fill() + scale_fill_manual(name=NULL, values="black") + 
  geom_raster(data = na.omit(encounter_prediction_df[c("UTMX","UTMY","NAs")]), 
              aes(x = UTMX, y = UTMY, fill = as.factor(NAs)))

#get uncertainty in final predictions (Figure S15)
ggplot() + geom_raster(data = encounter_prediction_df, 
                       aes(x = UTMX, y = UTMY, fill = combined_se)) +
  coord_equal() + scale_fill_gradientn(name = "", colours = viridis(100), limits = c(0,6))  + 
  geom_polygon(data = fortify(map), aes(group = group, x = long, y = lat)) +
  theme(text=element_text(size=16), axis.text=element_text(color="black")) + 
  coord_cartesian(xlim = c(75, 1050), ylim = c(4600, 5250)) +
  geom_polygon(data = fortify(map1), aes(group = group, x = long, y = lat)) + 
  ylab("Northings") + xlab("Eastings") +
  ggtitle(expression(paste("Log(Predicted CPUE [kg/",km^{2},"]) - Std. Errors"))) + 
  facet_wrap(~year, ncol = 4)

#export prediction and standard error rasters for external use (i.e., in ArcGIS... Figures 4 and 6)
sdmTMB_prediction_raster <- raster_create(encounter_prediction_df, "combined_predict", year_start = 2000)
crs(sdmTMB_prediction_raster) <- "+datum=WGS84 +proj=utm +zone=20T +units=km" 
writeRaster(sdmTMB_prediction_raster, paste(getwd(),"/Model Output/sdmTMB/sdmTMB_predictions.tif", sep = ""),
            format = "GTiff", bylayer = T, suffix = seq(2000,2019, by = 1), overwrite = T)
se_raster <- raster_create(encounter_prediction_df, "combined_se", year_start = 2000)
crs(se_raster) <- "+datum=WGS84 +proj=utm +zone=20T +units=km" 
writeRaster(se_raster, paste(getwd(),"/Model Output/sdmTMB/sdmTMB_standard_errors.tif", sep = ""),
            format = "GTiff", bylayer = T, suffix = seq(2000,2019, by = 1), overwrite = T)

#save model and raster predictions as .RData
save(encounter_sdmTMB, CPUE_sdmTMB, sdmTMB_prediction_raster, file = "sdmTMB.RData")