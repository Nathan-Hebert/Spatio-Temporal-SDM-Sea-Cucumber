#load necessary libraries
library(starve)
library(sdmTMB)
library(raster)
library(reshape)
library(ggplot2)
library(dplyr)
library(scales)
library(viridis)

######################FIGURE 1 - PLOTS OF SPATIAL DOMAIN, TOWS, AND DEM###############################

#load in pre-processed survey data
survey_data_combined <- read.csv(paste(getwd(),"/Data/survey_data_combined.csv", sep = ""))

#load in DEM data
DEM = raster("Data/TASSE_attributes_MaritimeDEM/BathyCHS_GEBCO_SEAM_mixedData_MartitimeExtentClip_100m_LatLong.asc")
#crop to spatial domain and cut out part around Sable that we don't want
DEM <- mask(crop(DEM, shp), shp)
DEM <- reclassify(DEM, cbind(-23, 0, NA), left = TRUE, right = TRUE)

#load in spatial domain shape file
shp <- shapefile("Data/SpatialDomain/MaritimesRegionEcosystemAssessmentStrata_SSsubset_ForSeaCuke.shp")

#get world data
world <- map_data('world') 
n_am <- world[world$region %in% c('Canada','USA'), ] #subset North America

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

################FIGURE 7 - PLOTS COMPARING ESTIMATES FROM sdmTMB and starve###########################

#CPUE SUB-MODELS - DEPTH (Figure 7a)
load(paste(getwd(),"/Workspace Files/sdmTMB.RData", sep = ""))
source("R Files/Helper Functions.R")
#plot effect estimate from sdmTMB model with 95% confidence bounds via 
#Breheny and Burchett's contrast approach (2017)
par(mar = c(5,5,4,2), mfrow = c(1,1))
effects_plot(CPUE_model, 2, 3, original_vector = survey_data_combined_CPUE$DEM_logScaled,
             xlab = "Log(Depth)", ylab = "Partial Effect", cex.main = 2,
             cex.axis = 2.45, cex.lab = 3, ylim = c(-5,2), line_color = "blue", 
             ci.color = "light blue", lwd = 3)
#add starve model line and 95% confidence bounds via Breheny and Burchett's contrast approach (2017)
load(paste(getwd(),"/Workspace Files/starve.RData", sep = ""))
x <- seq(min(survey_data_combined_CPUE$DEM_logScaled), max(survey_data_combined_CPUE$DEM_logScaled), 
         length.out = 1000)
y <- CPUE_model@TMB_out@sdr$par.fixed[6]*x+CPUE_model@TMB_out@sdr$par.fixed[7]*x^2
se <- sqrt(x^2*CPUE_model@tracing@parameter_covariance[6,6]+
             x^4*CPUE_model@tracing@parameter_covariance[7,7]+
             2*x^3*presence_model@tracing@parameter_covariance[6,7]) 
polygon(c(x, rev(x)), c(y-1.96*se, rev(y+1.96*se)),
        col = adjustcolor("orange", alpha.f=0.3), lty = 0)
lines(x,y, lwd = 3, col = "darkorange2")
#add legend differentiating the models
legend("topright",legend=c(expression(paste(italic('sdmTMB')," Model")), 
                        expression(paste(italic('starve')," Model"))), 
       fill = c("blue","darkorange2"), cex = 3)

#PRESENCE SUB-MODELS - BTM TEMP RANGE (Figure 7b)
load(paste(getwd(),"/Workspace Files/sdmTMB.RData", sep = ""))
source("R Files/Helper Functions.R")
#plot effect estimate from sdmTMB model with 95% confidence bounds via 
#Breheny and Burchett's contrast approach (2017)
par(mar = c(5,5,4,2), mfrow = c(1,1))
effects_plot(presence_model, 6,7, original_vector = survey_data_combined$RangeTempScaled,
             xlab = "Bottom Temperature Range", ylab = "Partial Effect", cex.main = 2, 
             cex.axis = 2.45, cex.lab = 3, ylim = c(-4,2.5), line_color = "blue", 
             ci.color = "light blue", lwd = 3)
#add starve model line and 95% confidence bounds via Breheny and Burchett's contrast approach (2017)
load(paste(getwd(),"/Workspace Files/starve.RData", sep = ""))
x <- seq(min(survey_data_combined$RangeTempScaled), max(survey_data_combined$RangeTempScaled), 
         length.out = 1000)
y <- presence_model@TMB_out@sdr$par.fixed[9]*x+presence_model@TMB_out@sdr$par.fixed[10]*x^2
se <- sqrt(x^2*presence_model@tracing@parameter_covariance[9,9]+
             x^4*presence_model@tracing@parameter_covariance[10,10]+
             2*x^3*presence_model@tracing@parameter_covariance[9,10]) 
polygon(c(x, rev(x)), c(y-1.96*se, rev(y+1.96*se)),
        col = adjustcolor("orange", alpha.f=0.3), lty = 0)
lines(x,y, lwd = 3, col = "darkorange2")
#add legend differentiating the models
legend("topright",legend=c(expression(paste(italic('sdmTMB')," Model")), 
                           expression(paste(italic('starve')," Model"))), 
       fill = c("blue","darkorange2"), cex = 3)