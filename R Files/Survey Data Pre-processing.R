#import necessary libraries
library(raster)
library(sp)
source("R Files/Helper Functions.R")

###LOAD IN BNAM DATA###

#grab list of desired BNAM files
full_file_list <- list.files("Data/BNAM/BNAM_2000-2019_monthly", full.names = TRUE, 
                             recursive = TRUE)
BtmTempFiles <- full_file_list[grepl("BtmTemp",full_file_list)==TRUE&
                                 grepl(".asc",full_file_list)==TRUE]
BtmSalinityFiles <- full_file_list[grepl("BtmSalinity",full_file_list)==TRUE&
                                     grepl(".asc",full_file_list)==TRUE]
BtmStressFiles <- full_file_list[grepl("BtmStress",full_file_list)==TRUE&
                                   grepl(".asc",full_file_list)==TRUE]

#load in spatial domain shape file
shp <- shapefile("Data/SpatialDomain/MaritimesRegionEcosystemAssessmentStrata_SSsubset_ForSeaCuke.shp")
#import the raster stacks of interest
BtmTempStack <- crop(stack(BtmTempFiles[121:360]), shp)
BtmSalinityStack <- crop(stack(BtmSalinityFiles[121:360]), shp)
BtmStressStack <- crop(stack(BtmStressFiles[121:360]), shp)

#calculate the mean, min, max, and range for bottom temperature
AverageTemp <- calc(BtmTempStack, fun = mean)
MinTemp <- calc(BtmTempStack, fun = min)
MaxTemp <- calc(BtmTempStack, fun = max)
RangeTemp <- MaxTemp - MinTemp
#calculate the mean, min, max, and range for bottom salinity
AverageSalinity <- calc(BtmSalinityStack, fun = mean)
MinSalinity <- calc(BtmSalinityStack, fun = min)
MaxSalinity <- calc(BtmSalinityStack, fun = max)
RangeSalinity <- MaxSalinity - MinSalinity
#calculate the mean, min, max, and range for bottom stress
AverageStress <- calc(BtmStressStack, fun = mean)
MinStress <- calc(BtmStressStack, fun = min)
MaxStress <- calc(BtmStressStack, fun = max)
RangeStress <- MaxStress - MinStress

#grab the average min, max, and range for bottom temperature
years <- seq(2000, 2019, by = 1)
average_min_max_range_temp <- average_min_max_range_calc(BtmTempFiles[121:360], years)
AverageMinTemp <- average_min_max_range_temp[[1]]
AverageMaxTemp <- average_min_max_range_temp[[2]]
AverageRangeTemp <- average_min_max_range_temp[[3]]
#grab the average min, max, and range for bottom salinity
years <- seq(2000, 2019, by = 1)
average_min_max_range_salinity <- average_min_max_range_calc(BtmSalinityFiles[121:360], years)
AverageMinSalinity <- average_min_max_range_salinity[[1]]
AverageMaxSalinity <- average_min_max_range_salinity[[2]]
AverageRangeSalinity <- average_min_max_range_salinity[[3]]
#grab the average min, max, and range for bottom stress
average_min_max_range_stress <- average_min_max_range_calc(BtmStressFiles[121:360], years)
AverageMinStress <- average_min_max_range_stress[[1]]
AverageMaxStress <- average_min_max_range_stress[[2]]
AverageRangeStress <- average_min_max_range_stress[[3]]

###LOAD IN DEM DATA (AND DERIVED LAYERS)###
DEM = raster("Data/TASSE_attributes_MaritimeDEM/BathyCHS_GEBCO_SEAM_mixedData_MartitimeExtentClip_100m_LatLong.asc")
Easterness = raster("Data/TASSE_attributes_MaritimeDEM/TASSE_Attributes/east_bathych/w001001.adf")
LocalMean = raster("Data/TASSE_attributes_MaritimeDEM/TASSE_Attributes/mean_bathych/w001001.adf")
Northerness = raster("Data/TASSE_attributes_MaritimeDEM/TASSE_Attributes/nort_bathych/w001001.adf")
Slope = raster("Data/TASSE_attributes_MaritimeDEM/TASSE_Attributes/slop_bathych/w001001.adf")
StandardDeviation = raster("Data/TASSE_attributes_MaritimeDEM/TASSE_Attributes/stde_bathych/w001001.adf")
RDMV = raster("Data/TASSE_attributes_MaritimeDEM/TASSE_Attributes/RDMV_bathych/w001001.adf")

###COMBINE THE TWO SURVEY DATA FRAMES###

#load in data from the two surveys...
RVseacukes <- read.csv(paste(getwd(),"/Data/RVsurvey.2000to2020.seacukes.Nov232021.final.csv", 
                             sep = ""))
snowcrab <- read.csv(paste(getwd(),"/Data/snowcrab_data_clean_Jan242022.csv", sep = ""))

#add important columns to RV survey data
RVseacukes$presence <- ifelse(is.na(RVseacukes$std.WGT) == FALSE | 
                                is.na(RVseacukes$std.No) == FALSE, TRUE, FALSE)
RVseacukes$season_winter <- ifelse(RVseacukes$month == 2 | RVseacukes$month == 3 | 
                                     RVseacukes$month == 4, 1, 0)
#preprocess snow crab survey data
snowcrab$presence <- ifelse(((is.na(snowcrab$Code6600.no) == FALSE) & (snowcrab$Code6600.no !=0)) 
                            | ((is.na(snowcrab$Code6600.wgt) == FALSE) & (snowcrab$Code6600.wgt != 0)), 
                            TRUE, FALSE)

#convert weight units for RV survey
RVseacukes$std.WGT_km_per_tow <- RVseacukes$std.WGT/(3.241)
RVseacukes$std.WGT <- RVseacukes$std.WGT_km_per_tow/(0.0124968)
#convert counts for RV survey
RVseacukes$std.No_km_per_tow <- RVseacukes$std.No/(3.241)
RVseacukes$std.No <- RVseacukes$std.No_km_per_tow/(0.0124968)

#combine the two sets of data
snowcrab$snowcrab <- TRUE
colnames(snowcrab)[c(6,7,10,11)] <- c("mid.lat","mid.lon","std.WGT","std.No") 
survey_data_combined <- smartbind(RVseacukes, snowcrab)
survey_data_combined$snowcrab <- ifelse(is.na(survey_data_combined$snowcrab) == TRUE, FALSE, TRUE)
#do other preprocessing
survey_data_combined$DEM_RDMV <- ifelse(is.na(survey_data_combined$DEM_RDMV), 0, 
                                        survey_data_combined$DEM_RDMV)
survey_data_combined$month_name <- month.abb[survey_data_combined$month]

#restrict spatial extent
coordinates <- SpatialPoints(cbind(survey_data_combined$mid.lon, survey_data_combined$mid.lat))
proj4string(coordinates) <- "+proj=longlat +datum=WGS84 +no_defs "
survey_data_combined <- survey_data_combined[which(is.na(over(coordinates, 
                                                              as(shp,"SpatialPolygons")))==FALSE),]
#drop 2020 tows as there is no BNAM data available for them currently
survey_data_combined <- survey_data_combined[-which(survey_data_combined$year == 2020),]

#extract the locations and convert to UTM-coordinates
XY = data.frame(survey_data_combined$mid.lon, survey_data_combined$mid.lat)
names(XY) = c("X","Y")
latLon = SpatialPoints(XY, proj4string = CRS(as.character(NA)))
proj4string(latLon) <- CRS("+proj=longlat +datum=WGS84") 
survey_data_combined[c("UTMX","UTMY")] <- as.data.frame(spTransform(latLon,
                                                                    paste("+proj=utm +zone=20T +units=km",
                                                                          sep = "")))

#remove tow in middle of domain that likely recorded sea cucumber weights in error
survey_data_combined <- survey_data_combined[-which(survey_data_combined$MISSIONSET=="TEM2008830_94"),]
#remove tow on slopes that likely recorded sea cucumber weights in error
survey_data_combined <- survey_data_combined[-which(survey_data_combined$MISSIONSET=="NED2010027_221"),]

###ADD BNAM/DEM TO COMBINED SURVEY DATA###

#add temporally/geographically matched BNAM values to the data frame
survey_data_combined$BtmTempBNAM <- match_BNAM(survey_data_combined, BtmTempStack, 
                                               "mid.lon", "mid.lat")
survey_data_combined$BtmSalinityBNAM <- match_BNAM(survey_data_combined, BtmSalinityStack, 
                                                   "mid.lon", "mid.lat")
survey_data_combined$BtmStressBNAM <- match_BNAM(survey_data_combined, BtmStressStack, 
                                                 "mid.lon", "mid.lat")

#add other derived temp layers
survey_data_combined$AverageTemp <- extract(AverageTemp, cbind(survey_data_combined$mid.lon, 
                                                               survey_data_combined$mid.lat))
survey_data_combined$MinTemp <- extract(MinTemp, cbind(survey_data_combined$mid.lon, 
                                                       survey_data_combined$mid.lat))
survey_data_combined$MaxTemp <- extract(MaxTemp, cbind(survey_data_combined$mid.lon, 
                                                       survey_data_combined$mid.lat))
survey_data_combined$RangeTemp <- extract(RangeTemp, cbind(survey_data_combined$mid.lon, 
                                                           survey_data_combined$mid.lat))
survey_data_combined$AverageMinTemp <- extract(AverageMinTemp, cbind(survey_data_combined$mid.lon, 
                                                                     survey_data_combined$mid.lat))
survey_data_combined$AverageMaxTemp <- extract(AverageMaxTemp, cbind(survey_data_combined$mid.lon, 
                                                                     survey_data_combined$mid.lat))
survey_data_combined$AverageRangeTemp <- extract(AverageRangeTemp, cbind(survey_data_combined$mid.lon, 
                                                                         survey_data_combined$mid.lat))

#add other derived salinity layers
survey_data_combined$AverageSalinity <- extract(AverageSalinity, cbind(survey_data_combined$mid.lon, 
                                                                       survey_data_combined$mid.lat))
survey_data_combined$MinSalinity <- extract(MinSalinity, cbind(survey_data_combined$mid.lon, 
                                                               survey_data_combined$mid.lat))
survey_data_combined$MaxSalinity <- extract(MaxSalinity, cbind(survey_data_combined$mid.lon, 
                                                               survey_data_combined$mid.lat))
survey_data_combined$RangeSalinity <- extract(RangeSalinity, cbind(survey_data_combined$mid.lon, 
                                                                   survey_data_combined$mid.lat))
survey_data_combined$AverageMinSalinity <- extract(AverageMinSalinity, cbind(survey_data_combined$mid.lon, 
                                                                             survey_data_combined$mid.lat))
survey_data_combined$AverageMaxSalinity <- extract(AverageMaxSalinity, cbind(survey_data_combined$mid.lon, 
                                                                             survey_data_combined$mid.lat))
survey_data_combined$AverageRangeSalinity <- extract(AverageRangeSalinity, cbind(survey_data_combined$mid.lon, 
                                                                                 survey_data_combined$mid.lat))

#add other derived stress layers
survey_data_combined$AverageStress <- extract(AverageStress, cbind(survey_data_combined$mid.lon, 
                                                                   survey_data_combined$mid.lat))
survey_data_combined$MinStress <- extract(MinStress, cbind(survey_data_combined$mid.lon, 
                                                           survey_data_combined$mid.lat))
survey_data_combined$MaxStress <- extract(MaxStress, cbind(survey_data_combined$mid.lon, 
                                                           survey_data_combined$mid.lat))
survey_data_combined$RangeStress <- extract(RangeStress, cbind(survey_data_combined$mid.lon, 
                                                               survey_data_combined$mid.lat))
survey_data_combined$AverageMinStress <- extract(AverageMinStress, cbind(survey_data_combined$mid.lon, 
                                                                         survey_data_combined$mid.lat))
survey_data_combined$AverageMaxStress <- extract(AverageMaxStress, cbind(survey_data_combined$mid.lon, 
                                                                         survey_data_combined$mid.lat))
survey_data_combined$AverageRangeStress <- extract(AverageRangeStress, cbind(survey_data_combined$mid.lon, 
                                                                             survey_data_combined$mid.lat))

#attach DEM to the data frame
survey_data_combined$DEM <- extract(DEM, cbind(survey_data_combined$mid.lon, 
                                               survey_data_combined$mid.lat))
survey_data_combined$DEM_Easterness <- extract(Easterness, cbind(survey_data_combined$mid.lon, 
                                                                 survey_data_combined$mid.lat))
survey_data_combined$DEM_LocalMean <- extract(LocalMean, cbind(survey_data_combined$mid.lon, 
                                                               survey_data_combined$mid.lat))
survey_data_combined$DEM_Northerness <- extract(Northerness, cbind(survey_data_combined$mid.lon, 
                                                                   survey_data_combined$mid.lat))
survey_data_combined$DEM_Slope <- extract(Slope, cbind(survey_data_combined$mid.lon, 
                                                       survey_data_combined$mid.lat))
survey_data_combined$DEM_StandardDeviation <- extract(StandardDeviation, cbind(survey_data_combined$mid.lon, 
                                                                               survey_data_combined$mid.lat))
survey_data_combined$DEM_RDMV <- extract(RDMV, cbind(survey_data_combined$mid.lon, 
                                                     survey_data_combined$mid.lat))

#if RDMV is NA, set to zero
survey_data_combined$DEM_RDMV <- ifelse(is.na(survey_data_combined$DEM_RDMV), 0, 
                                        survey_data_combined$DEM_RDMV)

###WRITE TO CSV###

write.csv(survey_data_combined, paste(getwd(),"/Data/survey_data_combined.csv", sep = ""))