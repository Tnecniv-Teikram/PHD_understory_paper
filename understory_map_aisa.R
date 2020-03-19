# ------------------This Script calculates understory BRF as a function of fertility class based on an upsampled 10m Aisa image ----------------# 

library(raster)
library(sp)
library(sf)
library(rgdal)
library(ggplot2)
library(dplyr)

# ---------------------------------------------------------------------------------------------------------------#
# Import data
#---------------------------------------------------------------------------------------------------------------#

# set file directory
setwd("C:/Users/VincentMarkiet/OneDrive - Advian Oy/Tiedostot/Personal/Education/phd")

# load vector data
fieldplots <- shapefile(x = "data/hyytiala_field_data/shapefiles/field_plots_titta_lauri_10m_buffer.shp")
#fieldplots <- st_read("data/hyytiala_field_data/shapefiles/field_plots_titta_lauri_10m_buffer.shp")
fielddata <- data.frame(fieldplots)
plot_U3_area <- st_read("data/hyytiala_field_data/shapefiles/plot_U3.shp")

# subset fieldplot
plot_test <- fieldplots[fieldplots$ID == "U3",]

# load Aisa raster data
aisa <- brick("C:/Users/VincentMarkiet/OneDrive - Advian Oy/Tiedostot/Personal/Education/phd/data/raster_data/AISA_20150703mosaic_10m")

#plotRGB(aisa, r = 90, g = 65, b = 35, stretch = 'hist') # inspect Aisa brick

# read Envi file header
header <- paste(readLines("C:/Users/VincentMarkiet/OneDrive - Advian Oy/Tiedostot/Personal/Education/phd/data/raster_data/AISA_20150703mosaic_10m.hdr"), collapse=" ")
header = strsplit(header, "}")
header = strsplit(header[[1]], "=")

#---------------------------------------------------------------------------------------------------------------#
# Pre-process Aisa & wavelength header files
#---------------------------------------------------------------------------------------------------------------#
WL=header[[6]][2]
WL=gsub("[^0-9., ]","",WL) #remove non-numeric, but leaves ".", comma and space
WL=strsplit(WL, ",")
WL=as.numeric(WL[[1]])
WL = WL*1000

band_name= header[[5]][2]
band_name=gsub("[^a-zA-Z0-9., ]","",band_name) #remove non-numeric, but leaves ".", comma and space
band_name=strsplit(band_name, ",")
band_name=band_name[[1]]
band_name <- gsub('Mosaic', "", band_name) # remove "mosaic" from band name
band_name <- gsub(" ","", band_name) 

# remove water vapor bands (753 - 763 nm) from wavelength variable
WL_noWV_1 <- WL[1:77]; # subset first set of bands
WL_noWV_2 <- WL[81:100]; # subset second set of bands
WL_noWV <- c(WL_noWV_1,WL_noWV_2); # combine the two sets
rm(WL_noWV_1, WL_noWV_2) # remove temp variable

# remove decimals from wavelength list
WL_noWV <- round(WL_noWV,0)
aisa_1_100 <- subset(aisa,1:100) # create subset of first 100 spectral bands (858nm - 991.43nm)
aisa_no_wv <- dropLayer(aisa_1_100,c(77,78,79)) # remove water vapor bands (749nm,753nm,758nm)
aisa <- aisa_no_wv
rm(aisa_1_100) # remove temp variable 

# raster_plot <- crop(aisa, plot_U3_area) # crop raster

#---------------------------------------------------------------------------------------------------------------#
# Create NDVI mask
#---------------------------------------------------------------------------------------------------------------#

# 
# # create forest mask using ndvi 0.7    # source for NDVI bands: https://www.researchgate.net/post/How_to_calculate_NDVI_using_all_possible_hyperspectral_band_combination
# ndvi <- (raster_plot$Mosaic..Band.100...0.858530. - raster_plot$Mosaic..Band.64...0.687550.) / (raster_plot$Mosaic..Band.100...0.858530. + raster_plot$Mosaic..Band.64...0.687550.) # calculate NDVI
# m <- ndvi
# m[m < 0.7] <- NA # mask values 

# mask all non-forest pixels using the NDVI mask
#raster_plot_mskd <- mask(raster_plot, m, maskvalue = NA)

#---------------------------------------------------------------------------------------------------------------#
# Extract plot BRF values for each pixel
#---------------------------------------------------------------------------------------------------------------#
#brf <- extract(aisa,fieldplots %>% filter(ID == "U18"),df = T, snap = 'near')
#brf <- extract(aisa, plot_U3_area, df = T, snap = 'near') # plot_test or 
#brf <- brf[,2:98] # remove first plot ID column
#brf <- brf / 10000
#save(brf, file = "C:/Users/VincentMarkiet/OneDrive - Advian Oy/Tiedostot/Personal/Education/phd/data/Rdata/aisa_spectra_plot.Rdata")
load(file = "C:/Users/VincentMarkiet/OneDrive - Advian Oy/Tiedostot/Personal/Education/phd/data/Rdata/aisa_spectra_plot.Rdata")

#------------------- test for extracting pixel BRF for CGF using all plots ------------------------------------ #
#brf_all <- extract(aisa, fieldplots, na.rm = T, df = T)
#save(brf_all, file = "C:/Users/VincentMarkiet/OneDrive - Advian Oy/Tiedostot/Personal/Education/phd/data/Rdata/aisa_spectra_10m_all_plots.Rdata")

# load all  aisa BRF values for all plots 
load("C:/Users/VincentMarkiet/OneDrive - Advian Oy/Tiedostot/Personal/Education/phd/data/Rdata/aisa_spectra_10m_all_plots.Rdata")

brf_all[,2:(ncol(brf_all))] <- brf_all[,2:(ncol(brf_all))] / 10000 # resample BRF values to a range between 0-1 and skip the plot ID column.
# replace all ID numbers with plot ID names
brf_all$ID <- fielddata$ID[match(brf_all$ID, unique(brf_all$ID))]

#---------------------------------------------------------------------------------------------------------------#
# clean up fieldplot database 
#---------------------------------------------------------------------------------------------------------------#
# First make sure  both brf dataframe and fieldplot dataframe have similar dimensions
# copy all fieldplot data to each pixel that is part of that same plot. 

brf_all_IDs <- data.frame(brf_all$ID )
colnames(brf_all_IDs) <- "ID"
plots <- merge(fielddata,brf_all_IDs ,by.x = "ID", by.y = "ID",sort = F) #, incomparables = NA

# Find all cells in fieldplot and mean plot BRF dataframe that have 0 values and NA values and return a dataframe with a TRUE or FALSE value.  
gooddata <- complete.cases(plots) & complete.cases(brf_all) 

# remove incomplete cells (empty or NA) 
plots <- plots[gooddata,] 
brf <- brf_all[gooddata,]

# create a list of erroneous plots. See excel: "list_of_erroneous_plots_titta_Lauri_hyytiala.xlsx" for reasons of removal.
error_list <- c("A2","B2","B3", "D5", "D6","D7", "E7","G1","U12","U13","U21","U24","U25",
                "61","111","114","125",
                "312","313","321","332","333",
                "411","412","413","421","425","431","441",
                "614","622","632","641",
                "712","721","723","725","734",
                "811","813","821","822","823","825","832","833","834","845",
                "912","925","931","932","933","934","941","942","943",
                "1014","1023","1032","1041","1042","1043",
                "1115","1141","1142","1145",
                "1211","1223","1232","1233","1234","1242","1243","1244",
                "1311","1322","1323","1331",
                "1411","1412","1413","1414","1422","1424","1431","1433","1434","1444",
                "1511","1512","1513","1521","1522","1532","1534","1543","1544","1545",
                "1622","1633","1641","1642")

# remove all erroneous plots from plot dataframe
for (i in (1:nrow(plots))){
  plots <- plots[! plots$ID %in% error_list,] # remove bad plots from field plot dataframe
  brf <- brf[! brf$ID %in% error_list,] # remove bad plots from mean plot BRF dataframe
  # brf_all <- brf_all[! brf$ID %in% error_list,] # remove bad plots from all pixel BRF database # this is commented out due to large processing time. The cleaned dataframe is saved and loaded in the next lines.
}


#---------------------------------------------------------------------------------------------------------------#
# Estimate sunlit fraction 
#---------------------------------------------------------------------------------------------------------------#

# import reference albedo data for estimating sunlit fraction 
#spec_ref <- read.csv("C:/Users/VincentMarkiet/OneDrive - Advian Oy/Tiedostot/Personal/Education/phd/IDL_scripts/leafspectra.csv", sep = "\t", header = T) # this prospect leaf albedo was discarded because it showed absorption features in the NIR region which shouldn't be there. 
#colnames(spec_ref)[1] <- "wavelength"
spec_ref <- read.csv("C:/Users/VincentMarkiet/OneDrive - Advian Oy/Tiedostot/Personal/Education/phd/IDL_scripts/prospect_Matti_november_2018/refalbedo_PROSPECTCp.txt", sep = " ", header = T, col.names = c("wavelength", "PROSPECT") ) 
spec_ref_prospect <- spec_ref$PROSPECT
wl_ref <- spec_ref$wavelength# wavelength array of reference spectrum
wl_image <- WL_noWV[69:82] # wavelengths of image between 710 and 790

# subset aisa sunlit understory brf to 709-788nm
BRF_subset <- brf[,69:82]  # create subset for wavelengths > 709nm and < 788 nm

# interpolate reference spectra to image wavelength
specref_int<- approx(wl_ref,spec_ref_prospect, wl_image)
specref <- specref_int$y

# convert solar zenith angle from degrees to radians
sol_zen_rad <- 48.80 / 180 * pi  

slf_df <- data.frame() # create empty dataframe to store sunlit fraction, slope, intercept and R2

# estimate slope, intercept of mean TOC BRF and mean sunlit fraction for each plot in fieldplot dataframe.

# for (i in 1:nrow(plots)) { # loop over plots
#   
#   for (j in 1:nrow(BRF_subset)){ # loop over plot pixels
#   
#   # test for NAs  
#   all(is.na(j))    
#   # estimate canopy sunlit fraction using prospect leaf albedo 
#   yy <- unlist(BRF_subset[j,]/ specref)    # estimate ratio of brf from leaf albedo 
#   bk <- lm(yy ~ unlist(BRF_subset[j,] ) ) # fit ordinary linear model between ratio of TOC BRF and leaf albedo. 
#   
#   # results
#   coefs_slf <- coef(bk) # linear model coeffecients
#   slope_slf <- coefs_slf[2] # slope
#   intercept_slf <- coefs_slf[1] # intercept
#   r2_slf <- round(summary(bk)$r.squared,2) # r2 or gooddness of fit
#   
#   # calculate sunlit fraction
#   slf <- 4 * cos(sol_zen_rad) * intercept_slf #  reference: 6, equation A7, page 5114  
#   slf_df <- rbind(slf_df, c(slope_slf,intercept_slf,r2_slf,slf)) # add row for each result 
# 
#   }
# }

for (i in 1:nrow(plots)) { 
  
  # estimate canopy sunlit fraction using prospect leaf albedo 
  yy <- unlist(BRF_subset[i,]/ specref)    # estimate ratio of brf from leaf albedo 
  bk <- lm(yy ~ unlist(BRF_subset[i,] ) ) # fit ordinary linear model between ratio of TOC BRF and leaf albedo. 
  
  # results
  coefs_slf <- coef(bk) # linear model coeffecients
  slope_slf <- coefs_slf[2] # slope
  intercept_slf <- coefs_slf[1] # intercept
  r2_slf <- round(summary(bk)$r.squared,2) # r2 or gooddness of fit
  
  # calculate sunlit fraction
  slf <- 4 * cos(sol_zen_rad) * intercept_slf #  reference: 6, equation A7, page 5114  
  slf_df <- rbind(slf_df, c(slope_slf,intercept_slf,r2_slf,slf)) # add row for each result 
  
}


# for (i in 1:nrow(plots)) { 
#   
#   # estimate canopy sunlit fraction using prospect leaf albedo 
#   
#   results_slf <- data.frame()
#   
#   s = BRF_subset[BRF_subset$ID == plots$ID[i],]
#   
#   for (j in 1:nrow(s)){
# 
#     yy <- unlist(s[j,]) / specref    # estimate ratio of brf from leaf albedo 
#     bk <- lm(yy ~ unlist(s[j,])) # fit ordinary linear model between ratio of TOC BRF and leaf albedo. 
#     # results
#     coefs_slf <- coef(bk) # linear model coeffecients
#     slope_slf <- coefs_slf[2] # slope
#     intercept_slf <- coefs_slf[1] # intercept
#     r2_slf <- round(summary(bk)$r.squared,2) # r2 or gooddness of fit
#     
#     # calculate sunlit fraction
#     slf <- 4 * cos(sol_zen_rad) * intercept_slf #  reference: 6, equation A7, page 5114      
#     
#     results_slf <- rbind(results_slf, c(slope_slf,intercept_slf,r2_slf,slf)) # add row for each result  )
#     colnames(results_slf) <- c("slope_slf", "intercept_slf","r2_slf", "sunlit_fraction")
#     
#   }
#   
#slf_df <- rbind(slf_df, results_slf) # add row for each result 

#}

colnames(slf_df) <- c("slope_slf", "intercept_slf","r2_slf", "sunlit_fraction") # add column name for each result

slf_df <- data.frame(slf_df,plots$ID) # add plot IDs to sunlit fraction dataframe

#----------------------------------------------------------------------------------------------------------------------------------------#
# calculate canopy gap fraction for plot 
#----------------------------------------------------------------------------------------------------------------------------------------#
#load libraries
library(betareg)
#library(gamlss)
#library(reshape)

source(file = "C:/Users/VincentMarkiet/OneDrive - Advian Oy/Tiedostot/Personal/Education/phd/backup_rscript/understory_reflectance/functions/error_bars_function.R", echo = FALSE)
# loads error bar function to plot error bars for the following wavelengths: 552nm, 645nm, 682nm, 739nm, 820
 
# references 

# 1) Hernandez et al.2016 IEEE Trans. Geosci. Remote Sens. 2016, 54, 5105–5116.
# 2) Knyazikhin, Y. et al., 2013. Proceedings of the National Academy of Sciences of the United States of America, 110(3), pp.E185-92.


results_r2 <- data.frame() # create empty dataframe to store results

# loop over each wavelength to estimate the goodness of fit (R2) between canopy gap fraction from field measurements
# and TOC BRF values from Aisa data. 

for (j in 1:(ncol(brf))){ # loop over each wavelength
  print(colnames(results_r2))
  brf_j <- brf[,2:ncol(brf)]
  brf_i <- c(brf_j[,j]) # get the BRF values for all pixels for a specific wavelength (j)
  r2_brf_can_model <- lm( brf_i ~ plots$gaps1) # fit a linear model between BRF and canopy gaps field measurements
  r2_can <- round( summary(r2_brf_can_model)$r.squared,2 ) # estimate R2 for each pixel
  
  r2_brf_can_model_gaps2 <- lm( brf_i ~ plots$gaps2) # fit a linear model between BRF and canopy gaps field measurements
  r2_can_gaps2 <- round( summary(r2_brf_can_model_gaps2)$r.squared,2 ) # estimate R2 for each pixel
  
  r2_brf_can_model_gaps3 <- lm( brf_i ~ plots$gaps3) # fit a linear model between BRF and canopy gaps field measurements
  r2_can_gaps3 <- round( summary(r2_brf_can_model_gaps3)$r.squared,2 ) # estimate R2 for each pixel
  
  r2_brf_can_model_gaps4 <- lm( brf_i ~ plots$gaps4) # fit a linear model between BRF and canopy gaps field measurements
  r2_can_gaps4 <- round( summary(r2_brf_can_model_gaps4)$r.squared,2 ) # estimate R2 for each pixel
  
  r2_brf_can_model_gaps5 <- lm( brf_i ~ plots$gaps5) # fit a linear model between BRF and canopy gaps field measurements
  r2_can_gaps5 <- round( summary(r2_brf_can_model_gaps5)$r.squared,2 ) # estimate R2 for each pixel
  
  results_r2 <- rbind(results_r2 ,c(r2_can,r2_can_gaps2,r2_can_gaps3,r2_can_gaps4,r2_can_gaps5) ) # add R2 for each pixel to dataframe
  
}

colnames(results_r2) <- c("R2_gaps1","R2_gaps2","R2_gaps3","R2_gaps4","R2_gaps5") # add column names

results_r2$WV <- WL_noWV

# print highest R2 values for each CGF view angle
R2_bands_gaps1 <- results_r2[results_r2$R2_gaps1 == max(results_r2$R2_gaps1),c(6,1)]
R2_bands_gaps2 <- results_r2[results_r2$R2_gaps2 == max(results_r2$R2_gaps2),c(6,2)]
R2_bands_gaps3 <- results_r2[results_r2$R2_gaps3 == max(results_r2$R2_gaps3),c(6,3)]
R2_bands_gaps4 <- results_r2[results_r2$R2_gaps4 == max(results_r2$R2_gaps4),c(6,4)]
R2_bands_gaps5 <- results_r2[results_r2$R2_gaps5 == max(results_r2$R2_gaps5),c(6,5)]


###################################################################################################################################################
# perform a beta regression between field estimated canopy cover (dependent variable) and BRF (621nm) + sunlit fraction 
###################################################################################################################################################
# select spectral bands with highest correlation with canopy gap fraction.
BRF_621 <- brf[,51]
BRF_612 <- brf[,49]


beta_can_brf_model_gaps1 <- betareg(plots$gaps1 ~ BRF_621 + slf_df$slope_slf + slf_df$intercept_slf) # fit a linear model
beta_can_brf_model_gaps2 <- betareg(plots$gaps2 ~ BRF_621 + slf_df$slope_slf + slf_df$intercept_slf) # fit a linear model
beta_can_brf_model_gaps3 <- betareg(plots$gaps3 ~ BRF_621 + slf_df$slope_slf + slf_df$intercept_slf) # fit a linear model
beta_can_brf_model_gaps4 <- betareg(plots$gaps4 ~ BRF_612 + slf_df$slope_slf + slf_df$intercept_slf) # fit a linear model


index_gaps5 <- plots$gaps5 == 0
plots$gaps5[index_gaps5] <- 0.000000001
beta_can_brf_model_gaps5 <- betareg(plots$gaps5 ~ BRF_612 + slf_df$slope_slf + slf_df$intercept_slf) # fit a linear model

# fit a model
fitted_model_gaps1 <- lm(fitted(beta_can_brf_model_gaps1) ~ plots$gaps1) # get the predicted values for the model
fitted_model_gaps2 <- lm(fitted(beta_can_brf_model_gaps2) ~ plots$gaps2) # get the predicted values for the model
fitted_model_gaps3 <- lm(fitted(beta_can_brf_model_gaps3) ~ plots$gaps3) # get the predicted values for the model
fitted_model_gaps4 <- lm(fitted(beta_can_brf_model_gaps4) ~ plots$gaps4) # get the predicted values for the model
fitted_model_gaps5 <- lm(fitted(beta_can_brf_model_gaps5) ~ plots$gaps5) # get the predicted values for the model

# plot measured canopy gap fraction against explanatory variables. 
plot(plots$gaps1 , fitted(beta_can_brf_model_gaps1) , main = " Measured canopy gap fraction vs TOC BRF 621nm + SLF slope + SLF intercept",
     xlab = "Field measured CGF", ylab = "Predicted CGF",
     xlim = c(0,1), ylim = c(0,1))
abline(fitted_model_gaps1, col = "red") # fit a line using the estimated values.
legend("topleft", bty = "n", paste("R^2",round(summary(fitted_model_gaps1)$adj.r.squared),3))
R2 <- summary(fitted_model_gaps1)$adj.r.squared
legend=bquote(atop(italic(R)^2 == .(format(R2, digits = 2)),"n" == .(nrow(plots$gaps1))))


CGF_DF <- data.frame(plots$gaps1,fitted(beta_can_brf_model_gaps1),plots$gaps2,fitted(beta_can_brf_model_gaps2),
                     plots$gaps3,fitted(beta_can_brf_model_gaps3),plots$gaps4,fitted(beta_can_brf_model_gaps4),
                     plots$gaps5, fitted(beta_can_brf_model_gaps5))

# change dataframe column names 
colnames(CGF_DF) <- c("CGF_field_gaps1","CGF_pred_BR_gaps1","CGF_field_gaps2","CGF_pred_BR_gaps2","CGF_field_gaps3","CGF_pred_BR_gaps3",
                      "CGF_field_gaps4","CGF_pred_BR_gaps4","CGF_field_gaps5","CGF_pred_BR_gaps5")

# add plot ID names for future reference
CGF_DF$ID <- plots$ID

# load MASS library 

# create function to plot scatterplot and fit a linear line
BR_CGF_plot <- function(CGF_field_gaps,CGF_pred_BR){
  
  m <- lm(CGF_pred_BR ~ CGF_field_gaps)  
  r2 <- as.character( format(summary(m)$r.squared, digits = 3) )
  
  # now plot linear model and R2
  ggplot(CGF_DF, aes(x = CGF_field_gaps, y = CGF_pred_BR )) +
    geom_point(col = "red") + geom_smooth(method = lm, col = "blue") + # plot data points
    xlab("Canopy gap fraction measured") + ylab("Canopy gap fraction predicted")  + xlim(0,(max(CGF_field_gaps)+0.1)) + ylim(0,(max(CGF_pred_BR))+ 0.1) +
    geom_text(x = 0.0, y = max(CGF_field_gaps), label = paste("r2", r2)) + 
    theme(plot.title = element_text(size = 18) , plot.background = element_rect(fill = "white", colour = "white"))  # + ggtitle("Measured CGF vs predicted CGF) center title above plot
  
}

BR_CGF_plot(CGF_DF$CGF_field_gaps1,CGF_DF$CGF_pred_BR_gaps1) # Canopy gap fraction observed vs predicted 1st view angle. View angle rings correspond to the LAI-2000 device to measure canopy gap fraction.
BR_CGF_plot(CGF_DF$CGF_field_gaps2,CGF_DF$CGF_pred_BR_gaps2) # ... 2nd view angle
BR_CGF_plot(CGF_DF$CGF_field_gaps3,CGF_DF$CGF_pred_BR_gaps3) # ... 3rd view angle
BR_CGF_plot(CGF_DF$CGF_field_gaps4,CGF_DF$CGF_pred_BR_gaps4) # ... 4th view angle
BR_CGF_plot(CGF_DF$CGF_field_gaps5,CGF_DF$CGF_pred_BR_gaps5) # ... 5th view angle


#----------------------------------------------------------------------------------------------------------------------------------------#
# estimate forest understory 
#----------------------------------------------------------------------------------------------------------------------------------------#

# load library needed to calculate all possible unique combinations of classes
library("gtools")

# references 

# 1) Hernandez et al.2016 IEEE Trans. Geosci. Remote Sens. 2016, 54, 51055116.
# 2) Knyazikhin, Y. et al., 2013. Proceedings of the National Academy of Sciences of the United States of America, 110(3), pp.E185-92.

###### 1) create subset

# create index of subset plots based on fertility class and tree species fraction.
# create subset indices 
plots_fert_1_index <- which( plots$Kasvupaikk == 1 ) #&  c(plots$pine_fraction > 0.9 | plots$spruce_fraction > 0.9 | plots$birch_fraction > 0.9)) 
plots_fert_2_index <- which( plots$Kasvupaikk == 2 ) #&  c(plots$pine_fraction > 0.9 | plots$spruce_fraction > 0.9 | plots$birch_fraction > 0.9)) 
plots_fert_3_index <- which( plots$Kasvupaikk == 3 ) #&  c(plots$pine_fraction > 0.9 | plots$spruce_fraction > 0.9 | plots$birch_fraction > 0.9)) 
plots_fert_4_index <- which( plots$Kasvupaikk == 4 ) #&  c(plots$pine_fraction > 0.9 | plots$spruce_fraction > 0.9 | plots$birch_fraction > 0.9)) 
plots_fert_5_index <- which( plots$Kasvupaikk == 5 ) #&  c(plots$pine_fraction > 0.9 | plots$spruce_fraction > 0.9 | plots$birch_fraction > 0.9)) 
# plots_fert_6_index <- which( plots$Kasvupaikka == 6 ) #&  c(plots$pine_fraction > 0.9 | plots$spruce_fraction > 0.9 | plots$birch_fraction > 0.9)) 

# create list with subsetplot indices 
indexlist <- list( plots_fert_1_index,  plots_fert_2_index,  plots_fert_3_index,  plots_fert_4_index, plots_fert_5_index ) # ,  plots_fert_6_index

##### 2) calculate BRF understory for subset plots

# create an empty list for results
plot_BRFu_list <- list()

# create empty list for bright pixel sample size tables
canopy_brf_list <- list()

# create list for albedo results
plot_bright_alb_est <- list()

# create list of fertility class names
fert_class_list <- c("herb-rich","moderately rich","moist upland","dryish upland","dry upland") # removed class ,"nutrient poor upland"  due to lack of field plot data for this fertility class

# set canopy gap fraction input. Either from field measurements (LAI-2000) or predicted canopy gap fraction using Hadi's Betaregression between BRF in red and measured canopy gap fraction.

# predicted
can_gaps1 <- CGF_DF$CGF_pred_BR_gaps1 #
can_gaps2 <- CGF_DF$CGF_pred_BR_gaps2 #
can_gaps3 <- CGF_DF$CGF_pred_BR_gaps3 #
can_gaps4 <- CGF_DF$CGF_pred_BR_gaps4 #
can_gaps5 <- CGF_DF$CGF_pred_BR_gaps5 #

#  interpolate leaf prospect data to aisa image wavelengths
specref_int_aisa<- approx(wl_ref,spec_ref_prospect, WL_noWV) 
rc <- specref_int_aisa$y

####### temp remove after ####
# p.ref = mean(slf_df$slope_slf) # use the mean p value using the mean p value for all plots
# p.eff = 0.65
# p.rel = (p.ref - p.eff) / (1 - p.eff)
# 
# rc <- (1 - p.rel) * rc / 1 - p.rel * rc 
##############################

# set brightest pixel bright_treshold. e.g. 0.8 = Use 20% brightest pixels
bright_treshold <- 0.9

##### input parameters for BRF understory model

# set up parameters used for calculating canopy gap fractions in the solar and sensor direction. 
x <- c(7,23,38,53,68) # LAI-2000 viewing angles
theta <- 48.8 # solar zenith angle (degrees)

# create dataframe to store p, LAI, tree species fraction values for each plot 
p_alb_results <- data.frame()

for ( i_index in 1:length(indexlist)){ 
  
  i_fert <- indexlist[[i_index]] # loop over list of plot indices 
  
  # create empty dataframe for temporary results
  plot_result <- data.frame() # create empty dataframe
  
  bright_pix_table <- data.frame() # create empty dataframe
  
  albedo_est_plots <- data.frame() # create empty dataframe for calculated leaf albedos
  
  canopy_brf_plot <- data.frame() # ....... for calculated canopy BRF
  
  
  for ( i_plot in (i_fert)){ # loop over each plot
    
    # 1) Use theory of spectral invariants to calculate leaf albedo using the brightest plot pixels, PROSPECT leaf albedo, and pixel sunlit fraction 
    # suggested by Knyazikhin et al. 2013 to use PROSPECT derived reference leaf albedo to make the model insensitive to tree species. Reference 2)
    
    # select all plot pixels
    #plot_pixels_mean <- colMeans(brf[which(brf$ID == plots$ID[i_plot] ),2:ncol(brf)])
    
    # subset brightest 10% pixels to calculate leaf albedo spectrum
    
    # 1.1) select plot pixels using field plot ID
    plot_pixels <- brf_all[which(brf_all$ID == plots$ID[i_plot]),2:ncol(brf_all)] # skip column 1 (ID column)
    
    # 1.2) Calculate the mean of brightest 10% pixels
    bright_TOC_mean <- colMeans(subset(plot_pixels, plot_pixels$Mosaic..Band.85...0.787090. >= quantile(plot_pixels$Mosaic..Band.85...0.787090., bright_treshold)))
    
    # correct for actual optical properties (bark etc). For this, retrieve p and rho, and invert the spectral invariant equation for leaf albedo for all wavelengths
    # retrieve p and rho using prospect leaf albedo 
    yy <- unlist(bright_TOC_mean[69:82] / rc[69:82])    # estimate ratio of brf from leaf albedo using wavelength 710-790 nm
    bk <- lm(yy ~ unlist(bright_TOC_mean[69:82] ) ) # fit ordinary linear model between ratio of TOC BRF and leaf albedo. 
    coefs_slf <- coef(bk) # linear model coeffecients
    slope_slf <- coefs_slf[2] # slope / p
    intercept_slf <- coefs_slf[1] # intercept / rho
    aa <- slope_slf * bright_TOC_mean + intercept_slf # Eq. 5 in manuscript Markiet,Mottus .2019
    rc <- as.data.frame(t(bright_TOC_mean  / aa )) #store estimated leaf albedo in dataframe and tranpose rows to columns.
    
    
    
    # -------------------------------------------------------------------------------- #
    
    # set sunlit fraction
    slf <- slf_df$sunlit_fraction[i_plot]
    
    
    # 2) Calculate canopy gap fractions for sensor direction 
    
    Ti <- c(can_gaps1[i_plot], can_gaps2[i_plot], can_gaps3[i_plot],
            can_gaps4[i_plot], can_gaps5[i_plot]) # field measured canopy gap fraction 
    
    
    
    # interpolate viewing angles of LAI-2000  to the solar zenith angle to get gap fraction at solar zenith angle
    sz_gap_fr <- approx(x , Ti, theta)
    t0 <- sz_gap_fr$y #  gap fraction in solar zenith direction:  incoming radiation through canopy to forest floor (direct transmittance). Use transmittance measurements.
    tv <- can_gaps1[i_plot] # gap fraction in sensor direction
    
    # 2) set four reflectance conversion factors (k) and four proportions (f).
    
    ##### estimate fraction of canopy and understory visible in a pixel  ###########
    
    # 2.1) understory fractions      
    fu.sl <- t0 * tv
    fu.sh <- (1 - t0) * tv
    
    # 2.2) canopy fractions
    # alpha_slC =  i_slf[i_plot,4] # fraction of sunlit canopy foliage in all visible canopy foliage: assuming non-green understory
    #alpha_slC = (slf - t0 * tv ) / ( 1 - tv ) # fraction of sunlit canopy foliage in all visible canopy foliage: assuming green understory
    #fc.sl <- alpha_slC  # fraction of sunlit canopy in a pixel  ( we assume that fc.sl & fc.sh are approx. == sunlit fraction)
    #fc.sh <-  (1 - alpha_slC )  # fraction of shaded canopy in a pixel
    
    #fc.sl <- slf    # fraction of sunlit canopy in non-green understory scenario. Calculated as the ratio of total sunlit fraction to visible canopy in a pixel. 
    fc.sl <- slf - t0 * tv # we assume that the plot has green understory
    fc.sh <- 1 - fu.sh - fu.sl - fc.sl # fraction of shaded canopy where all fractions should add up to 1.
    
    # 2.3) calculate all the proportions of reflectance for each component   
    kc.sl <- slf_df$slope_slf[i_plot] * plot_pixels_mean + 1/( 4 * cos(sol_zen_rad) ) #           irradiance for sunlit canopy
    kc.sh <- slf_df$slope_slf[i_plot] * plot_pixels_mean   # FC / Fin                 ...... shaded canopy 
    
    
    # the effect of various p-values in ku.sh on understory BRF. This is done to test to correct for the difference in true leaf albedo and PROSPECT p.  
    
    # 2.4) calculate understory 
    p_test_slope_slf <- 1
    
    # calculate a p value that can be used independent of tree species. Eq. 22 in underestory BRF manuscript Markiet,Mottus. 2019 
    ku.sh <- p_test_slope_slf * plot_pixels_mean # irradiance for shaded understory # 
    ku.sl <- ku.sh + 1 # .... sunlit understory
    
    # 3) Calculate canopy BRF
    
    # 3.1) break down the canopy component
    BRF_can <-  rc * ( kc.sh * fc.sh   +  kc.sl * fc.sl ) # * (1 - can_gaps1[i_plot] ) #  field_plots_meas$gaps1[i_plot]
    und_component <- (ku.sl * fu.sl + ku.sh * fu.sh) # understory
    zz <- plot_pixels_mean - BRF_can
    
    # 4) calculate understory BRF 
    ru <-  as.data.frame(zz / und_component) # estimate understory reflectance
    
    # 5)  store results
    
    ######## temp #############
    
    
    p_alb_results <- rbind(p_alb_results,"p" = c(round(slope_slf,2),"LAI" = plots$LAI[i_plot], "pine_fraction" = plots$pine_fraction[i_plot],
                                                 "spruce_fraction" = plots$spruce_fraction[i_plot], "birch_fraction" = plots$birch_fraction[i_plot]))
    rownumber = dim(plot_result)[1]
    rownames(p_alb_results)[rownumber] <- paste(brf$fieldplot_buffers.ID[i_plot] )
    
    #########################
    
    # 5.1) store results, one row per plot
    if (dim(plot_result)[1] == 0 ) { # insert check to test if there are more than 0 rows.
      plot_result <- ru
      albedo_est_plots <- rc
      canopy_brf_plot <- BRF_can
      
      
      
    } else { # Add rows to existing result dataframe if number of rows is higher than 0.
      #temp_data <- data.frame(ru,Nr_bright_pix) # create temporary vector to store BRFus and Nr of bright pixels used for each plot. 
      # colnames(temp_data) <- colnames(plot_result) # assign result dataframe column names to temporary dataframe so they can be merged later.
      plot_result <- rbind( plot_result,ru) # add temporary dataframe to result dataframe
      albedo_est_plots <- rbind(albedo_est_plots,rc)
      canopy_brf_plot <- rbind(canopy_brf_plot, BRF_can)
      
      
    }
    
    # 5.3) store plot ID as row number
    rownumber = dim(plot_result)[1]
    rownames(plot_result)[rownumber] <-  paste(brf$fieldplot_buffers.ID[i_plot] )
    #
    
    rownumber = dim(albedo_est_plots)[1]
    rownames(albedo_est_plots)[rownumber] <-  paste(brf$fieldplot_buffers.ID[i_plot] )
    
    rowwnumber = dim(canopy_brf_plot)[1]
    rownames(canopy_brf_plot)[rownumber] <- paste(brf$fieldplot_buffers.ID[i_plot])
    
  }
  
  # 5.4) store dataframe results in a list 
  plot_BRFu_list[[i_index]] <- plot_result # store results in list.
  canopy_brf_list[[i_index]] <- bright_pix_table # store dataframe with number of bright pixels in a list.
  plot_bright_alb_est[[i_index]] <- albedo_est_plots # store leaf albedo dataframe in list
  canopy_brf_list[[i_index]] <- canopy_brf_plot # add each result dataframe to list
  
  # assign fertility class names to each column 
  #colnames(canopy_brf_list[[i_index]]) <- fert_class_list[i_index] 
  
  ######## temp remove after #####
  #p_list[[i_index]] <- p_alb_results
  ################################2
}



# Create map based on each pixel 



# calculate RMSE for plots