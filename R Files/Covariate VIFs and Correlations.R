#load necessary libraries
library(corrplot); library(car)

#load in pre-processed survey data
survey_data_combined <- read.csv(paste(getwd(),"/Data/survey_data_combined.csv", sep = ""))

#correlation plot for temperature layers alone
cor <- cor(survey_data_combined[grep("Temp",names(survey_data_combined))])
rownames(cor) <- c("Temporally-matched Temp.", "Temp. Overall Average","Temp. Minimum", 
                   "Temp. Maximum", "Temp. Range", "Temp. Average Minimum", "Temp. Average Maximum"
                   , "Temp. Average Range")
colnames(cor) <- rownames(cor)
corrplot(cor, method = "number", type = "upper")
#correlation plot for salinity layers alone
cor <- cor(survey_data_combined[grep("Salinity",names(survey_data_combined))])
rownames(cor) <- c("Temporally-matched Salinity", "Salinity Overall Average","Salinity Minimum", 
                   "Salinity Maximum", "Salinity Range", "Salinity Average Minimum", "Salinity Average Maximum"
                   , "Salinity Average Range")
colnames(cor) <- rownames(cor)
corrplot(cor, method = "number", type = "upper")
#correlation plot for stress layers alone
cor <- cor(survey_data_combined[grep("Stress",names(survey_data_combined))])
rownames(cor) <- c("Temporally-matched Stress", "Stress Overall Average","Stress Minimum", 
                   "Stress Maximum", "Stress Range", "Stress Average Minimum", "Stress Average Maximum"
                   , "Stress Average Range")
colnames(cor) <- rownames(cor)
corrplot(cor, method = "number", type = "upper")

#correlation plot with predictors used in final models
survey_data_combined$RangeStressLog <- log(survey_data_combined$RangeStress)
survey_data_combined$BtmStressBNAMLog <- log(survey_data_combined$BtmStressBNAM)
survey_data_combined$sqrt_DEM_Slope <- sqrt(survey_data_combined$DEM_Slope)
survey_data_combined$DEM_log <- log(-survey_data_combined$DEM)
cor <- cor(survey_data_combined[c("BtmTempBNAM","BtmSalinityBNAM","BtmStressBNAMLog","RangeTemp",
                                "RangeStressLog","DEM_log","sqrt_DEM_Slope","DEM_RDMV", 
                                "DEM_Easterness","DEM_Northerness")])
rownames(cor) <- c("Btm. Temperature", "Btm. Salinity","Log(Btm. Stress)", 
                   "Btm. Temperature Range", "Log(Btm. Stress Range)", 
                   "Log(Depth)", "Square Root of Slope", "RDMV","Easterness", "Northerness")
colnames(cor) <- rownames(cor)
corrplot(cor, method = "number", type = "upper", tl.col = "black")

#check VIFs for final encounter models
b <- glm(presence~BtmTempBNAM+
           BtmSalinityBNAM+
           log(BtmStressBNAM)+
           log(RangeStress)+
           RangeTemp+
           log(-DEM)+
           sqrt(DEM_Slope)+
           DEM_RDMV+
           DEM_Easterness+
           DEM_Northerness+
           snowcrab, family = "binomial", data = survey_data_combined)
vif(b)

#check VIFs for final CPUE models
survey_data_combined_CPUE <- survey_data_combined[which(is.na
                                                        (survey_data_combined$std.WGT)==FALSE&
                                                          survey_data_combined$std.WGT!=0),]
b <- glm(log(std.WGT)~BtmTempBNAM+
           BtmSalinityBNAM+
           log(BtmStressBNAM)+
           log(RangeStress)+
           RangeTemp+
           log(-DEM)+
           sqrt(DEM_Slope)+
           DEM_RDMV+
           DEM_Easterness+
           DEM_Northerness+
           snowcrab, family = "gaussian", data = 
           survey_data_combined_CPUE)
vif(b)