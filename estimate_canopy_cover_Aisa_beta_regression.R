# estimate canopy cover using AISA BRF and beta regression.  reference: Hadi et al. (2016)

# steps of script:
# 1) load input data: Mean Top of canopy BRF values (3-7-2015) for 254 plots (20m buffer pine, spruce, birch) for Hyytiälä, Finland.
# 2) Estimate which spectral band has the highest correlation with measured canopy gap fraction.
# 3) Estimate canopy cover with TOC BRF, sunlit fraction slope and intercept using a beta regression model.

#load library
library(betareg)
library(ggplot2)
library(gamlss)
library(reshape)

##################################################################################################################
# 1) load input data: Aisa top of BRF_can_subset (TOC) BRF for 20m buffers for pine, spruce and birch plots. 

source('C:/Users/VincentMarkiet/OneDrive - Advian Oy/Tiedostot/Personal/Education/phd/backup_rscript/understory_reflectance/import_aisa_data_for_understory_brf.R',
       encoding = 'UTF-8', echo=FALSE)

##################################################################################################################       


##################################################################################################################       
# 2)  Estimate which spectral band has the highest correlation with measured canopy gap fraction.

# create dataframe for results
results_r2 <- data.frame()

# loop over each wavelength to estimate the goodness of fit (R2) between canopy gap fraction from field measurements
# and TOC BRF values from Aisa data. 

for (j in 1:(ncol(brf)-1)){ # loop over each wavelength
  
  brf_i <- c(brf[,j]) # get the BRF values for all plots for a specific wavelength (j)
  r2_brf_can_model <- lm( brf_i ~ plots$gaps1) # fit a linear model between BRF and canopy gaps field measurements
  r2_can <- round( summary(r2_brf_can_model)$r.squared,2 ) # estimate R2 for each plot
  
  r2_brf_can_model_gaps2 <- lm( brf_i ~ plots$gaps2) # fit a linear model between BRF and canopy gaps field measurements
  r2_can_gaps2 <- round( summary(r2_brf_can_model_gaps2)$r.squared,2 ) # estimate R2 for each plot
  
  r2_brf_can_model_gaps3 <- lm( brf_i ~ plots$gaps3) # fit a linear model between BRF and canopy gaps field measurements
  r2_can_gaps3 <- round( summary(r2_brf_can_model_gaps3)$r.squared,2 ) # estimate R2 for each plot
  
  r2_brf_can_model_gaps4 <- lm( brf_i ~ plots$gaps4) # fit a linear model between BRF and canopy gaps field measurements
  r2_can_gaps4 <- round( summary(r2_brf_can_model_gaps4)$r.squared,2 ) # estimate R2 for each plot
  
  r2_brf_can_model_gaps5 <- lm( brf_i ~ plots$gaps5) # fit a linear model between BRF and canopy gaps field measurements
  r2_can_gaps5 <- round( summary(r2_brf_can_model_gaps5)$r.squared,2 ) # estimate R2 for each plot
  
  results_r2 <- rbind(results_r2 ,c(r2_can,r2_can_gaps2,r2_can_gaps3,r2_can_gaps4,r2_can_gaps5) ) # add R2 for each plot to dataframe
  
}

# change column names
colnames(results_r2) <- c("R2_gaps1","R2_gaps2","R2_gaps3","R2_gaps4","R2_gaps5")

# plot R2 values as a function of wavelength for gaps in the solar and sensor direction
# for (i in 1:(ncol(results_r2))){
#   plot(type = "p", x = WL_noWV , y = results_r2[,i],
#        col = i, xlab = "Wavelength (nm)",
#        ylab = "R2", main = paste("R2 as function of wavelength: BRF vs CGF", colnames(results_r2[i]), sep = " "))
#   lines(x = WL_noWV, y = rep(mean(results_r2[,i]),97), col = "black") # draw the mean R2 line
# }

# add wavelengths to R2 dataframe
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
BRF_621 <- brf[,50]
BRF_612 <- brf[,48]

# create the beta regression model: dependent variable field measured canopy gaps, independent variables: sunlit fraction slope and intercept.
beta_can_brf_model_gaps1 <- betareg(plots$gaps1 ~ BRF_621 + slf_df$slope_slf + slf_df$intercept_slf) # fit a linear model
beta_can_brf_model_gaps2 <- betareg(plots$gaps2 ~ BRF_621 + slf_df$slope_slf + slf_df$intercept_slf) # fit a linear model
beta_can_brf_model_gaps3 <- betareg(plots$gaps3 ~ BRF_621 + slf_df$slope_slf + slf_df$intercept_slf) # fit a linear model
beta_can_brf_model_gaps4 <- betareg(plots$gaps4 ~ BRF_612 + slf_df$slope_slf + slf_df$intercept_slf) # fit a linear model

# gaps5 has 0 values which are not supported by the betaregression. Thus change 0 values to tiny values
index_gaps5 <- plots$gaps5 == 0
plots$gaps5[index_gaps5] <- 0.000000001
beta_can_brf_model_gaps5 <- betareg(plots$gaps5 ~ BRF_612 + slf_df$slope_slf + slf_df$intercept_slf) # fit a linear model

# plot coefficients 
#coefficients(beta_can_brf_model_gaps1)

# fit a model
fitted_model_gaps1 <- lm(fitted(beta_can_brf_model_gaps1) ~ plots$gaps1) # get the predicted values for the model
fitted_model_gaps2 <- lm(fitted(beta_can_brf_model_gaps2) ~ plots$gaps2) # get the predicted values for the model
fitted_model_gaps3 <- lm(fitted(beta_can_brf_model_gaps3) ~ plots$gaps3) # get the predicted values for the model
fitted_model_gaps4 <- lm(fitted(beta_can_brf_model_gaps4) ~ plots$gaps4) # get the predicted values for the model
fitted_model_gaps5 <- lm(fitted(beta_can_brf_model_gaps5) ~ plots$gaps5) # get the predicted values for the model

# plot measured canopy gap fraction against explanatory variables. 
# plot(plots$gaps1 , fitted(beta_can_brf_model_gaps1) , main = " Measured canopy gap fraction vs TOC BRF 621nm + SLF slope + SLF intercept",
#       xlab = "Field measured CGF", ylab = "Predicted CGF",
#       xlim = c(0,1), ylim = c(0,1))
# abline(fitted_model_gaps1, col = "red") # fit a line using the estimated values.
legend("topleft", bty = "n", paste("R^2",round(summary(fitted_model_gaps1)$adj.r.squared),2))
R2 <- summary(fitted_model_gaps1)$adj.r.squared
legend=bquote(atop(italic(R)^2 == .(format(R2, digits = 2)),"n" == .(nrow(plots$gaps1))))

##### plot results using ggplot library #########
# combine results in a dataframe
CGF_DF <- data.frame(plots$gaps1,fitted(beta_can_brf_model_gaps1),plots$gaps2,fitted(beta_can_brf_model_gaps2),
                     plots$gaps3,fitted(beta_can_brf_model_gaps3),plots$gaps4,fitted(beta_can_brf_model_gaps4),
                     plots$gaps5, fitted(beta_can_brf_model_gaps5))

# change dataframe column names 
colnames(CGF_DF) <- c("CGF_field_gaps1","CGF_pred_BR_gaps1","CGF_field_gaps2","CGF_pred_BR_gaps2","CGF_field_gaps3","CGF_pred_BR_gaps3",
                      "CGF_field_gaps4","CGF_pred_BR_gaps4","CGF_field_gaps5","CGF_pred_BR_gaps5")

# add plot ID names for future reference
CGF_DF$ID <- plots$ID

# set plot border margin 
#par(mar = c(0, 0, 0, 0))

# create function to plot scatterplot and fit a linear line
# BR_CGF_plot <- function(CGF_field_gaps,CGF_pred_BR){
#   ggplot(CGF_DF, aes(x = CGF_field_gaps, y = CGF_pred_BR )) +
#    geom_point(col = "red") + geom_smooth(method = rlm, col = "blue") + # plot data points
#    xlab("Canopy gap fraction measured") + ylab("Canopy gap fraction predicted") +  # x and y axis labels
#    ggtitle("Measured CGF vs predicted CGF") + xlim(0,0.8) + ylim(0,1) + # plot title
#    theme(plot.title = element_text(size = 18) , plot.background = element_rect(fill = "white", colour = "white"))  # center title above plot
# }

# create function to plot scatterplot and fit a linear line
BR_CGF_plot <- function(CGF_field_gaps,CGF_pred_BR){

  m <- lm(CGF_pred_BR ~ CGF_field_gaps)  
  r2 <- as.character( format(summary(m)$r.squared, digits = 3) )

  # now plot linear model and R2
  ggplot(CGF_DF, aes(x = CGF_field_gaps, y = CGF_pred_BR )) +
    geom_point(col = "red") + geom_smooth(method = rlm, col = "blue") + # plot data points
    xlab("Canopy gap fraction measured") + ylab("Canopy gap fraction predicted")  + xlim(0,(max(CGF_field_gaps)+0.1)) + ylim(0,(max(CGF_pred_BR))+ 0.1) +
    geom_text(x = 25, y = 300, label = paste("r2", r2)) + 
    theme(plot.title = element_text(size = 18) , plot.background = element_rect(fill = "white", colour = "white"))  # + ggtitle("Measured CGF vs predicted CGF) center title above plot
  
}

# print all scatterplots field measured vs predicted canopy gap fraction
#BR_CGF_plot(CGF_DF$CGF_field_gaps1,CGF_DF$CGF_pred_BR_gaps1) # Canopy gap fraction observed vs predicted 1st view angle. View angle rings correspond to the LAI-2000 device to measure canopy gap fraction.
#BR_CGF_plot(CGF_DF$CGF_field_gaps2,CGF_DF$CGF_pred_BR_gaps2) # ... 2nd view angle
#BR_CGF_plot(CGF_DF$CGF_field_gaps3,CGF_DF$CGF_pred_BR_gaps3) # ... 3rd view angle
#BR_CGF_plot(CGF_DF$CGF_field_gaps4,CGF_DF$CGF_pred_BR_gaps4) # ... 4th view angle
#BR_CGF_plot(CGF_DF$CGF_field_gaps5,CGF_DF$CGF_pred_BR_gaps5) # ... 5th view angle


# end 
##############################################################################################################

