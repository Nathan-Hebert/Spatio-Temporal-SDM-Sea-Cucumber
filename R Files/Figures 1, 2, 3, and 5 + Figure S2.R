#load necessary libraries
library(starve); library(sdmTMB); library(mgcv)
library(raster); library(reshape); library(ggplot2)
library(dplyr); library(scales); library(viridis)
library(ggridges); library(ggthemes); library(rgdal)
source("R Files/Helper Functions.R")

#load in pre-processed survey data
survey_data_combined <- read.csv(paste(getwd(),"/Data/survey_data_combined.csv", sep = ""))

#load in the models and prediction rasters
load("mgcv.RData")
load("starve.RData")
load("sdmTMB.RData")

#load in spatial domain shape file
shp <- shapefile("Data/SpatialDomain/MaritimesRegionEcosystemAssessmentStrata_SSsubset_ForSeaCuke.shp")

#get world data
world <- map_data('world') 
n_am <- world[world$region %in% c('Canada','USA'), ] #subset North America

######################FIGURE 1 - PLOTS OF SPATIAL DOMAIN, TOWS, AND DEM############################

#load in DEM data
DEM = raster("Data/TASSE_attributes_MaritimeDEM/BathyCHS_GEBCO_SEAM_mixedData_MartitimeExtentClip_100m_LatLong.asc")
#crop to spatial domain and cut out part around Sable that we don't want
DEM <- mask(crop(DEM, shp), shp)
DEM <- reclassify(DEM, cbind(-23, 0, NA), left = TRUE, right = TRUE)

#plot wider study area (Figure 1a)
ggplot() + scale_x_continuous(breaks = c(-70,-63,-56)) + scale_y_continuous(breaks = c(41, 46, 51)) +
  theme(text=element_text(size=33), axis.text=element_text(color="black"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.key.size = unit(1.5,"line"), legend.text=element_text(size=33), 
        legend.position = c(0.184,0.86),
        legend.background = element_rect(color = "black"), legend.key = element_rect(fill = "light grey")) + 
  xlab("Longitude") + ylab("Latitude") + geom_polygon(data = shp, mapping = 
                                                        aes(x = long, y = lat, group = group), 
                                                      col = "orange", fill = NA, size = 1.25) + 
  coord_map(xlim = c(-71.5,-53), ylim = c(39,52)) +
  geom_polygon(data = n_am, mapping = aes(x = long, y = lat, group = group), fill = "grey44", col = "black") +
  ggsn::scalebar(x.min = -71.5, x.max = -53, y.min = 39, y.max = 52, dist = 350, transform = TRUE, 
                 model = "WGS84", height = 0.02, 
                 st.dist = 0.05, st.bottom = F, dist_unit = "km" , location = "bottomleft", st.size = 11) + 
  ggspatial::annotation_north_arrow(location = "tr", which_north = "true",
                                    pad_x = unit(-0.01, "in"), pad_y = unit(0.16, "in"), 
                                    height = unit(3.2, "cm"), width = unit(3.2, "cm"), 
                                    style = ggspatial::north_arrow_nautical(fill = c("grey40", "white"), 
                                                                            line_col = "grey20", 
                                                                            text_size = 33)) +
  theme(plot.margin = margin(0.25, 0.25, 0.25, 0.25, "cm")) + 
  annotate("rect", xmin = -68, xmax = -56, ymin = 41.1, ymax = 47.8, fill = NA, col = "black", size = 1.5)

#plot spatial domain and tows (delineated by survey) (Figure 1b)
survey_data_combined$Survey <- ifelse(survey_data_combined$snowcrab==FALSE, "RV","Snow Crab")
ggplot() + scale_x_continuous(breaks = c(-66,-62,-58)) +
  theme(text=element_text(size=33), axis.text=element_text(color="black"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.key.size = unit(1.5,"line"), legend.text=element_text(size=33), legend.position = c(0.184,0.86),
        legend.background = element_rect(color = "black"), legend.key = element_rect(fill = "light grey")) + 
  xlab("Longitude") + ylab("Latitude") + 
  geom_polygon(data = shp, mapping = aes(x = long, y = lat, group = group), fill = "black") + 
  geom_point(data = survey_data_combined, aes(x = mid.lon, y = mid.lat, col = Survey)) + 
  coord_map(xlim = c(-67.5,-56.5), ylim = c(41.2,47.5)) +
  geom_polygon(data = n_am, mapping = aes(x = long, y = lat, group = group), fill = "grey44", col = "black") +
  ggsn::scalebar(x.min = -67.5, x.max = -56.5, y.min = 41.2, y.max = 47.5, dist = 200, 
                 transform = TRUE, model = "WGS84", height = 0.02, 
                 st.dist = 0.04, st.bottom = F, dist_unit = "km" , location = "bottomleft", st.size = 11) + 
  ggspatial::annotation_north_arrow(location = "tr", which_north = "true",
                                    pad_x = unit(0.1, "in"), pad_y = unit(0.42, "in"), 
                                    height = unit(3.2, "cm"), width = unit(3.2, "cm"), 
                                    style = ggspatial::north_arrow_nautical(fill = c("grey40", "white"), 
                                                                            line_col = "grey20", 
                                                                            text_size = 33)) +
  theme(plot.margin = margin(0.25, 0.25, 0.25, 0.25, "cm"))+ 
  scale_color_manual(values = c(rgb(0,154/255,240/255),"orange")) + 
  guides(colour = guide_legend(override.aes = list(size = 5)))

#plot DEM (Figure 1c)
DEM <- aggregate(-DEM, fact = 20)
DEM_plot <- as.data.frame(DEM, xy = TRUE) %>% melt(id.vars = c('x','y'))
ggplot() + geom_tile(data = na.omit(DEM_plot), aes(x = x, y = y, fill = value)) + 
  scale_x_continuous(breaks = c(-66,-62,-58)) +
  scale_fill_gradientn(name = "Depth (m)", trans = "log", 
                       colours = paletteer::paletteer_c("ggthemes::Classic Blue", 14), na.value = NA, 
                       breaks = c(30,150,1100))+
  theme(text=element_text(size=33), axis.text=element_text(color="black"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.key.size = unit(1.5,"line"), legend.text=element_text(size=33), legend.position = c(0.137,0.82),
        legend.background = element_rect(color = "black")) + xlab("Longitude") + ylab("Latitude") + 
  geom_polygon(data = n_am, mapping = aes(x = long, y = lat, group = group), fill = "grey44", col = "black") + 
  coord_map(xlim = c(-67.5,-56.5), ylim = c(41.2,47.5)) + 
  ggsn::scalebar(x.min = -67.5, x.max = -56.5, y.min = 41.2, y.max = 47.5, dist = 200, 
                 transform = TRUE, model = "WGS84", height = 0.02, 
                 st.dist = 0.04, st.bottom = F, dist_unit = "km" , location = "bottomleft", st.size = 11) + 
  ggspatial::annotation_north_arrow(location = "tr", which_north = "true",
                                    pad_x = unit(0.1, "in"), pad_y = unit(0.42, "in"), 
                                    height = unit(3.2, "cm"), width = unit(3.2, "cm"), 
                                    style = ggspatial::north_arrow_nautical(fill = c("grey40", "white"), 
                                                                            line_col = "grey20", 
                                                                            text_size = 33)) + 
  theme(plot.margin = margin(0.25, 0.25, 0.25, 0.25, "cm"))

######################FIGURES 2 AND 3 - PARTIAL EFFECT PLOTS###########################

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

#grab encounter tows and set up response
survey_data_combined_CPUE <- survey_data_combined[which(is.na(survey_data_combined$std.WGT)==FALSE&
                                                          survey_data_combined$std.WGT!=0),]

#plot the encounter sub-model terms (Figure 2)
par(mfrow = c(4,3), mar = c(5,5,4,2))
effects_plot_combined(encounter_mgcv, encounter_starve, encounter_sdmTMB, mgcv_a = 2, sdmTMB_a = 2, 
                      sdmTMB_b = 3,  starve_a = 5, starve_b = 6, 
                      original_vector = survey_data_combined$DEM_logScaled, ylim = c(-10,1.3),
                      xlab = "Log(Depth)")
effects_plot_combined(encounter_mgcv, encounter_starve, encounter_sdmTMB, mgcv_a = 3, sdmTMB_a = 4, 
                      sdmTMB_b = 5, starve_a = 7, starve_b = 8, 
                      original_vector = survey_data_combined$BtmTempBNAMScaled, ylim = c(-2,3),
                      xlab = "Bottom Temperature")
effects_plot_combined(encounter_mgcv, encounter_starve, encounter_sdmTMB, mgcv_a = 4, sdmTMB_a = 6, 
                      sdmTMB_b = 7, starve_a = 9, starve_b = 10, 
                      original_vector = survey_data_combined$RangeTempScaled, ylim = c(-2,3),
                      xlab = "Btm. Temperature Range", sig = c("mgcv", "sdmTMB"))
effects_plot_combined(encounter_mgcv, encounter_starve, encounter_sdmTMB, mgcv_a = 5, sdmTMB_a = 8, 
                      sdmTMB_b = 9, starve_a = 11, starve_b = 12, 
                      original_vector = survey_data_combined$BtmSalinityBNAMScaled, ylim = c(-2,3),
                      xlab = "Bottom Salinity")
effects_plot_combined(encounter_mgcv, encounter_starve, encounter_sdmTMB, mgcv_a = 6, sdmTMB_a = 10, 
                      starve_a = 13, original_vector = survey_data_combined$BtmStressBNAMLogScaled, 
                      ylim = c(-2,3),xlab = "Log(Bottom Stress)")
effects_plot_combined(encounter_mgcv, encounter_starve, encounter_sdmTMB, mgcv_a = 7, sdmTMB_a = 11, 
                      sdmTMB_b = 12, starve_a = 14, starve_b = 15, 
                      original_vector = survey_data_combined$RangeStressLogScaled, ylim = c(-2,3),
                      xlab = "Log(Btm. Stress Range)", sig = c("mgcv", "sdmTMB"))
effects_plot_combined(encounter_mgcv, encounter_starve, encounter_sdmTMB, mgcv_a = 8, sdmTMB_a = 13, 
                      sdmTMB_b = 14, starve_a = 16, starve_b = 17, 
                      original_vector = survey_data_combined$sqrt_DEM_SlopeScaled, ylim = c(-2,3),
                      xlab = "Square Root of Slope", sig = c())
effects_plot_combined(encounter_mgcv, encounter_starve, encounter_sdmTMB, mgcv_a = 9, sdmTMB_a = 15, 
                      sdmTMB_b = 16, starve_a = 18, starve_b = 19, 
                      original_vector = survey_data_combined$DEM_NorthernessScaled, ylim = c(-2,3),
                      xlab = "Northerness", sig  = c("mgcv"))
effects_plot_combined(encounter_mgcv, encounter_starve, encounter_sdmTMB, mgcv_a = 10, sdmTMB_a = 17,
                      starve_a = 20, original_vector = survey_data_combined$DEM_EasternessScaled, 
                      ylim = c(-2,3), xlab = "Easterness", sig = c())
effects_plot_combined(encounter_mgcv, encounter_starve, encounter_sdmTMB, mgcv_a = 11, sdmTMB_a = 18, 
                      sdmTMB_b = 19, starve_a = 21, starve_b = 22, 
                      original_vector = survey_data_combined$DEM_RDMVScaled, ylim = c(-2,3),
                      xlab = "RDMV", sig = c())

#plot the CPUE sub-model terms (Figure 3)
par(mfrow = c(4,3), mar = c(5,5,4,2))
effects_plot_combined(CPUE_mgcv, CPUE_starve, CPUE_sdmTMB, mgcv_a = 2, sdmTMB_a = 2, sdmTMB_b = 3, 
                      starve_a = 6, starve_b = 7, 
                      original_vector = survey_data_combined_CPUE$DEM_logScaled, ylim = c(-4,1.3),
                      xlab = "Log(Depth)")
effects_plot_combined(CPUE_mgcv, CPUE_starve, CPUE_sdmTMB, mgcv_a = 3, sdmTMB_a = 4, sdmTMB_b = 5, 
                      starve_a = 8, starve_b = 9, 
                      original_vector = survey_data_combined_CPUE$BtmTempBNAMScaled, ylim = c(-2,3),
                      xlab = "Bottom Temperature", sig = c("starve","sdmTMB"))
effects_plot_combined(CPUE_mgcv, CPUE_starve, CPUE_sdmTMB, mgcv_a = 4, sdmTMB_a = 6, sdmTMB_b = 7, 
                      starve_a = 10, starve_b = 11, 
                      original_vector = survey_data_combined_CPUE$RangeTempScaled, ylim = c(-2,3),
                      xlab = "Btm. Temperature Range")
effects_plot_combined(CPUE_mgcv, CPUE_starve, CPUE_sdmTMB, mgcv_a = 5, sdmTMB_a = 8, sdmTMB_b = 9, 
                      starve_a = 12, starve_b = 13, 
                      original_vector = survey_data_combined_CPUE$BtmSalinityBNAMScaled, ylim = c(-2,3),
                      xlab = "Bottom Salinity", sig = c("starve","sdmTMB"))
effects_plot_combined(CPUE_mgcv, CPUE_starve, CPUE_sdmTMB, mgcv_a = 6, sdmTMB_a = 10, 
                      starve_a = 14, 
                      original_vector = survey_data_combined_CPUE$BtmStressBNAMLogScaled, ylim = c(-2,3),
                      xlab = "Log(Bottom Stress)", sig = c())
effects_plot_combined(CPUE_mgcv, CPUE_starve, CPUE_sdmTMB, mgcv_a = 7, sdmTMB_a = 11, sdmTMB_b = 12, 
                      starve_a = 15, starve_b = 16, 
                      original_vector = survey_data_combined_CPUE$RangeStressLogScaled, ylim = c(-2,3),
                      xlab = "Log(Btm. Stress Range)", sig = c())
effects_plot_combined(CPUE_mgcv, CPUE_starve, CPUE_sdmTMB, mgcv_a = 8, sdmTMB_a = 13, sdmTMB_b = 14, 
                      starve_a = 17, starve_b = 18, 
                      original_vector = survey_data_combined_CPUE$sqrt_DEM_SlopeScaled, ylim = c(-2,3),
                      xlab = "Square Root of Slope", sig = c("sdmTMB"))
effects_plot_combined(CPUE_mgcv, CPUE_starve, CPUE_sdmTMB, mgcv_a = 9, sdmTMB_a = 15, sdmTMB_b = 16, 
                      starve_a = 19, starve_b = 20, 
                      original_vector = survey_data_combined_CPUE$DEM_NorthernessScaled, ylim = c(-2,3),
                      xlab = "Northerness", sig = c())
effects_plot_combined(CPUE_mgcv, CPUE_starve, CPUE_sdmTMB, mgcv_a = 10, sdmTMB_a = 17,
                      starve_a = 21, 
                      original_vector = survey_data_combined_CPUE$DEM_EasternessScaled, ylim = c(-2,3),
                      xlab = "Easterness", sig = c())
effects_plot_combined(CPUE_mgcv, CPUE_starve, CPUE_sdmTMB, mgcv_a = 11, sdmTMB_a = 18, sdmTMB_b = 19, 
                      starve_a = 22, starve_b = 23, 
                      original_vector = survey_data_combined_CPUE$DEM_RDMVScaled, ylim = c(-2,3),
                      xlab = "RDMV", sig = c())

######################FIGURES 6 AND APPENDIX S1: FIGURE S2 - RESERVE PLOTS###########################

#load polygons for reserves
load("Data/reserves_shape_file.RData")

#setup vector to store the rasters
model <- c("mgcv","starve","sdmTMB")
prediction_rasters <- c(mgcv_prediction_raster, starve_prediction_raster, 
                        sdmTMB_prediction_raster)

#subset all models' predictions to the reserve areas, and combine these predictions
#into one dataframe
prediction_df_list <- list()
for (i in 1:length(prediction_rasters)) {
  
  #subset to reserve area
  raster <- mask(prediction_rasters[[i]], reserve)
  
  #turn into dataframe
  raster_df <- as.data.frame(raster, xy = TRUE) %>%
    melt(id.vars = c('x', 'y')) %>%
    filter(!is.na(value))
  levels(raster_df$variable) <- seq(2000, 2019, by = 1)
  raster_df$model <- model[i]
  
  #add dataframe to the list
  prediction_df_list[[length(prediction_df_list) + 1]] <- raster_df
}
combined_prediction_df <- do.call(rbind, prediction_df_list)

#reorder the 'model' variable factor levels
combined_prediction_df$model <- factor(
  combined_prediction_df$model,
  levels = c("mgcv", "starve", "sdmTMB")
)

#density plots of the predictions inside reserves for all years and models (Figure 6)
ggplot(combined_prediction_df, 
       aes(x = value, fill = model, color = model)) +
  geom_density(alpha = 0.17, size = 1.2) +
  theme_base() +
  theme(
    text = element_text(size = 20),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "top",
    legend.key = element_blank(),
    plot.background = element_rect(colour = NA, fill = NA, size = 5)
  ) +
  scale_fill_manual(
    name = "",
    values = c("mgcv" = "blue", "starve" = "black", "sdmTMB" = "red"),
    labels = c(
      expression(italic("mgcv") ~ "model"),
      expression(italic("starve") ~ "model"),
      expression(italic("sdmTMB") ~ "model")
    )
  ) +
  scale_color_manual(values = c("blue", "black", "red")) +
  ylab("Density") +
  xlab("Log(Predicted CPUE [kg/sq km]) Inside Reserves") +
  xlim(c(2, 11.5)) +
  facet_wrap(~as.factor(variable), scales = "free_y", ncol = 4) +
  guides(color = FALSE, fill = guide_legend(override.aes = list(alpha = 1, linetype = 0)))

#plot illustrating the location of the reserves (Appendix S1: Figure S2)
reserve <- spTransform(reserve, 
                       "+proj=longlat +datum=WGS84 +no_defs") #convert to latitude/longitude from UTM
ggplot() + scale_x_continuous(breaks = c(-66,-62,-58)) +
  theme(text=element_text(size=24), axis.text=element_text(color="black"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.key.size = unit(1.5,"line"), legend.text=element_text(size=24), legend.position = c(0.184,0.86),
        legend.background = element_rect(color = "black"), legend.key = element_rect(fill = "light grey")) + 
  xlab("Longitude") + ylab("Latitude") + 
  geom_polygon(data = fortify(reserve), aes(group = group, x = long, y = lat), 
               col = "black", fill = NA, size = 1) +
  geom_polygon(data = shp, mapping = aes(x = long, y = lat, group = group), fill = NA, 
               col = "orange", size = 1) +
  coord_map(xlim = c(-67.5,-56.5), ylim = c(41.2,47.5)) +
  geom_polygon(data = n_am, mapping = aes(x = long, y = lat, group = group), fill = "grey44", col = "black") +
  ggsn::scalebar(x.min = -67.5, x.max = -56.5, y.min = 41.2, y.max = 47.5, dist = 200, 
                 transform = TRUE, model = "WGS84", height = 0.02, 
                 st.dist = 0.04, st.bottom = F, dist_unit = "km" , location = "bottomleft", st.size = 7) + 
  ggspatial::annotation_north_arrow(location = "tr", which_north = "true",
                                    pad_x = unit(0.1, "in"), pad_y = unit(0.42, "in"), 
                                    height = unit(3.2, "cm"), width = unit(3.2, "cm"), 
                                    style = ggspatial::north_arrow_nautical(fill = c("grey40", "white"), 
                                                                            line_col = "grey20", 
                                                                            text_size = 24)) +
  theme(plot.margin = margin(0.25, 0.25, 0.25, 0.25, "cm"))+ 
  scale_color_manual(values = c(rgb(0,154/255,240/255),"orange")) + 
  guides(colour = guide_legend(override.aes = list(size = 5)))