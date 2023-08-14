#load necessary libraries
library(raster)
source("R Files/Helper Functions.R")

#increase memory limit
memory.limit(size = 400000)

#load in spatial domain shape file
shp <- shapefile("Data/SpatialDomain/MaritimesRegionEcosystemAssessmentStrata_SSsubset_ForSeaCuke.shp")

#grab list of desired BNAM files
full_file_list <- list.files("Data/BNAM/BNAM_1990-2019_monthly", full.names = TRUE, 
                             recursive = TRUE)
BtmTempFiles <- full_file_list[grepl("BtmTemp",full_file_list)==TRUE&
                                 grepl(".asc",full_file_list)==TRUE]
BtmSalinityFiles <- full_file_list[grepl("BtmSalinity",full_file_list)==TRUE&
                                     grepl(".asc",full_file_list)==TRUE]
BtmStressFiles <- full_file_list[grepl("BtmStress",full_file_list)==TRUE&
                                   grepl(".asc",full_file_list)==TRUE]

#load in DEM files
DEM = raster("Data/TASSE_attributes_MaritimeDEM/BathyCHS_GEBCO_SEAM_mixedData_MartitimeExtentClip_100m_LatLong.asc")
Easterness = raster("Data/TASSE_attributes_MaritimeDEM/TASSE_Attributes/east_bathych/w001001.adf")
Northerness = raster("Data/TASSE_attributes_MaritimeDEM/TASSE_Attributes/nort_bathych/w001001.adf")
Slope = raster("Data/TASSE_attributes_MaritimeDEM/TASSE_Attributes/slop_bathych/w001001.adf")
RDMV = raster("Data/TASSE_attributes_MaritimeDEM/TASSE_Attributes/RDMV_bathych/w001001.adf")

#crop and mask DEM rasters to our spatial domain
DEM <- mask(crop(DEM, shp), shp)
Easterness <- mask(crop(Easterness, shp), shp)
Northerness <- mask(crop(Northerness, shp), shp)
Slope <- mask(crop(Slope, shp), shp)
RDMV <- mask(crop(RDMV, shp), shp)

#get an average temp. raster for each year and then crop 
raster_data_BtmTemp <- yearly_average_raster_data(BtmTempFiles[121:360])
raster_data_BtmTemp <- mask(crop(raster_data_BtmTemp, shp), shp)
#get an average salinity raster for each year and then crop
raster_data_BtmSalinity <- yearly_average_raster_data(BtmSalinityFiles[121:360])
raster_data_BtmSalinity <- mask(crop(raster_data_BtmSalinity, shp), shp)
#get an average stress raster for each year and then crop 
raster_data_BtmStress <- yearly_average_raster_data(BtmStressFiles[121:360])
raster_data_BtmStress <- mask(crop(raster_data_BtmStress, shp), shp)

#create raster_datas for 2000-2019 bottom stress and crop to spatial domain
BtmStressraster_data <- stack(BtmStressFiles[121:360])
BtmStressraster_data <- mask(crop(BtmStressraster_data, shp), shp)
#create raster_datas for 2000-2019 bottom temp and crop to spatial domain
BtmTempraster_data <- stack(BtmTempFiles[121:360])
BtmTempraster_data <- mask(crop(BtmTempraster_data, shp), shp)
#create raster_datas for 2000-2019 bottom salinity and crop to spatial domain
BtmSalinityraster_data <- stack(BtmSalinityFiles[121:360])
BtmSalinityraster_data <- mask(crop(BtmSalinityraster_data, shp), shp)
#calculate the range for bottom stress, salinity, and temp
RangeStress <- calc(BtmStressraster_data, fun = max) - calc(BtmStressraster_data, fun = min)
RangeTemp <- calc(BtmTempraster_data, fun = max) - calc(BtmTempraster_data, fun = min)
RangeSalinity <- calc(BtmSalinityraster_data, fun = max) - calc(BtmSalinityraster_data, fun = min)

#project the rasters
raster_data_BtmTemp <- projectRaster(raster_data_BtmTemp, 
                                     crs = "+proj=utm +zone=20T +datum=WGS84 +units=km")
raster_data_BtmSalinity <- projectRaster(raster_data_BtmSalinity, 
                                         crs = "+proj=utm +zone=20T +datum=WGS84 +units=km")
raster_data_BtmStress <- projectRaster(raster_data_BtmStress, 
                                       crs = "+proj=utm +zone=20T +datum=WGS84 +units=km")
RangeStress <- projectRaster(RangeStress, 
                             crs = "+proj=utm +zone=20T +datum=WGS84 +units=km")
RangeTemp <- projectRaster(RangeTemp, 
                           crs = "+proj=utm +zone=20T +datum=WGS84 +units=km")
RangeSalinity <- projectRaster(RangeSalinity, 
                               crs = "+proj=utm +zone=20T +datum=WGS84 +units=km")
DEM <- projectRaster(DEM, 
                     crs = "+proj=utm +zone=20T +datum=WGS84 +units=km")
Easterness <- projectRaster(Easterness, 
                            crs = "+proj=utm +zone=20T +datum=WGS84 +units=km")
Northerness <- projectRaster(Northerness, 
                             crs = "+proj=utm +zone=20T +datum=WGS84 +units=km")
Slope <- projectRaster(Slope, 
                       crs = "+proj=utm +zone=20T +datum=WGS84 +units=km")
RDMV <- projectRaster(RDMV, 
                      crs = "+proj=utm +zone=20T +datum=WGS84 +units=km")

#cut out part of domain around Sable that we don't want
DEM <- reclassify(DEM, cbind(-23, 0, NA), left = TRUE, right = TRUE)

#aggregate DEM files
aggregate_factor = 60
DEM_timeseries <- aggregate(DEM, fact = aggregate_factor, na.rm = FALSE)
Easterness_timeseries <- aggregate(Easterness, fact = aggregate_factor)
Northerness_timeseries <- aggregate(Northerness, fact = aggregate_factor)
Slope_timeseries <- aggregate(Slope, fact = aggregate_factor)
RDMV_timeseries <- aggregate(RDMV, fact = aggregate_factor)

#make BNAM rasters comparable to DEM
raster_data_BtmTemp_timeseries <- resample(raster_data_BtmTemp, DEM_timeseries)
raster_data_BtmSalinity_timeseries <- resample(raster_data_BtmSalinity, DEM_timeseries)
raster_data_BtmStress_timeseries <- resample(raster_data_BtmStress, DEM_timeseries)
RangeStress_timeseries <- resample(RangeStress, DEM_timeseries)
RangeTemp_timeseries <- resample(RangeTemp, DEM_timeseries)
RangeSalinity_timeseries <- resample(RangeSalinity, DEM_timeseries)

#create a raster for snowcrab indicator
snowcrab_ind_timeseries <- setValues(DEM_timeseries, TRUE)

#create a raster_data of predictor values for each year and convert the raster_data into a data frame...
#for 2019 only, change 1:20 to 20:20
raster_data_df_timeseries <- data.frame()
for (i in 1:20)
{
  pred_raster_data <- stack(raster_data_BtmTemp_timeseries[[i]], raster_data_BtmSalinity_timeseries[[i]], 
                            DEM_timeseries, Easterness_timeseries, Northerness_timeseries, Slope_timeseries, 
                            RDMV_timeseries, snowcrab_ind_timeseries, raster_data_BtmStress_timeseries[[i]], 
                            RangeTemp_timeseries, RangeStress_timeseries, RangeSalinity_timeseries)
  names(pred_raster_data) <- c("BtmTempBNAM","BtmSalinityBNAM","DEM", "DEM_Easterness",
                         "DEM_Northerness","DEM_Slope","DEM_RDMV", "snowcrab", "BtmStressBNAM",
                         "RangeTemp", "RangeStress","RangeSalinity")
  raster_data_df_new <- as.data.frame(na.omit(rasterToPoints(pred_raster_data)))
  raster_data_df_new$year <- 1999 + i #increment year
  raster_data_df_timeseries <- rbind(raster_data_df_timeseries, raster_data_df_new)
}
raster_data_df_timeseries$UTMX <- raster_data_df_timeseries$x
raster_data_df_timeseries$UTMY <- raster_data_df_timeseries$y
#trick to remove NA's caused by log
raster_data_df_timeseries$log_depth <- log(-raster_data_df_timeseries$DEM)
raster_data_df_timeseries <- na.omit(raster_data_df_timeseries)

#write to csv
write.csv(raster_data_df_timeseries, paste(getwd(),"/Data/raster_data_df_timeseries.csv", sep = ""))