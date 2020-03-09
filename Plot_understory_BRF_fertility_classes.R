
#----------------------- model to calculate understory BRF as a function of fertility class. ---------------------------------#

#----------- import data -----------------#
# source('C:/Users/vmvincent/OneDrive/PhD/R_analysis/understory_reflectance/Extract_bright_canopy_pixels_for_can_albedo.R')  
#source("c:/Users/vmvincent/OneDrive/PhD/R_analysis/understory_reflectance/import_aisa_data_for_understory_brf.R", echo = FALSE)
#source("c:/Users/vmvincent/OneDrive/PhD/R_analysis/understory_reflectance/estimate_canopy_cover_Aisa_beta_regression.R",echo = FALSE)

source("C:/Users/VincentMarkiet/OneDrive - Advian Oy/Tiedostot/Personal/Education/phd/backup_rscript/understory_reflectance/estimate_canopy_cover_Aisa_beta_regression.R")
source("C:/Users/VincentMarkiet/Documents/Test_rscript/", echo = FALSE)
source('C:/Users/vmvincent/OneDrive/PhD/R_analysis/understory_reflectance/functions/error_bars_function.R', echo=FALSE)
#-----------------------------------------#
# load library needed to calculate all possible unique combinations of classes
library("gtools")

# This script makes available the following variables:

# 1) BRF for all pixels for each plot  (data frame name: BRF_all_plots)
# 2) sunlit fraction, slope, intercept for each plot (dataframe name: slf_df)
# 3) plot dataframe which contains: plot ID, fertility class, tree species fraction, LAI, gap fraction, etc. (dataframe name: plots)
# 4) PROSPECT leaf albedo dataframe 
# 5) loads error bar function to plot error bars for the following wavelengths: 552nm, 645nm, 682nm, 739nm, 820
# 6) imports predicted canopy gap fraction

# references 

# 1) Hernandez et al.2016 IEEE Trans. Geosci. Remote Sens. 2016, 54, 5105–5116.
# 2) Knyazikhin, Y. et al., 2013. Proceedings of the National Academy of Sciences of the United States of America, 110(3), pp.E185-92.


###### 1) create subset

# create index of subset plots based on fertility class and tree species fraction.
# create subset indices 
plots_fert_1_index <- which( plots$Kasvupaikka == 1 ) #&  c(plots$pine_fraction > 0.9 | plots$spruce_fraction > 0.9 | plots$birch_fraction > 0.9)) 
plots_fert_2_index <- which( plots$Kasvupaikka == 2 ) #&  c(plots$pine_fraction > 0.9 | plots$spruce_fraction > 0.9 | plots$birch_fraction > 0.9)) 
plots_fert_3_index <- which( plots$Kasvupaikka == 3 ) #&  c(plots$pine_fraction > 0.9 | plots$spruce_fraction > 0.9 | plots$birch_fraction > 0.9)) 
plots_fert_4_index <- which( plots$Kasvupaikka == 4 ) #&  c(plots$pine_fraction > 0.9 | plots$spruce_fraction > 0.9 | plots$birch_fraction > 0.9)) 
plots_fert_5_index <- which( plots$Kasvupaikka == 5 ) #&  c(plots$pine_fraction > 0.9 | plots$spruce_fraction > 0.9 | plots$birch_fraction > 0.9)) 
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

# field measured 
# can_gaps1 <- plots$gaps1
# can_gaps2 <- plots$gaps2
# can_gaps3 <- plots$gaps3
# can_gaps4 <- plots$gaps4
# can_gaps5 <- plots$gaps5


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
    plot_pixels_mean <- colMeans(brf_all_plots[which(brf_all_plots$ID == plots$ID[i_plot] ),2:ncol(brf_all_plots)])
    
    # subset brightest 10% pixels to calculate leaf albedo spectrum
    
    # 1.1) select plot using plot ID
    plot_pixels <- brf_all_plots[which(brf_all_plots$ID == plots$ID[i_plot]),2:ncol(brf_all_plots)]
    
    # 1.2) Calculate the mean of  brightest 10% pixels
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

# give column names to p results dataframe
colnames(p_alb_results) <- c("p", "LAI", "pine_fraction", "spruce_fraction","birch_fraction")
# create scatter plot of p against LAI 

#-------------------------- plot results ------------------------------------------------------------------------------------------------------#

# set plot settings
#par(mfrow = c(1,2))  

# plot linegraph with BRF understory as a function of wavelength.
# plot(type = "l", x = WL_noWV, y = ru,
#      main = paste("rU plot",plots$ID[i_plot]," ","Fert:",plots$Kasvupaikka[i_plot], " ","CGF1", can_gaps1[i_plot],
#                   "TOC >",brightness.lvl, "n = ",Nr_bright_pix ),
#      xlab = "Wavelength (nm)", col = "red")  

# par(mfrow=c(2,3)) # set row,column plot settings
# 
# # plot for each fertility class BRF understory graphs 
# for ( i in (1:length(plot_BRFu_list))){ # loop over list of dataframes with each containing BRF understory for a fertility class
#   
#   if (nrow(plot_BRFu_list[[i]] > 0)){
#     
#     # first, plot the fertility group mean rU
#     plot(type = "l", x = WL_noWV, y = lapply(plot_BRFu_list,colMeans)[[i]][1:length(bright_TOC_mean)] , main = c(paste("rU: fertility class",i,fert_class_list[i] )),
#          ylim = c(0,1), xlab = "Wavelength (nm)", ylab = "BRF", col = 1, lty = 2, lwd = 2) # plot for each fertility class a new plot
#     
#     for ( j in (1:10)){ # add rU reflectance lines for all plots  
#       
#       lines(x = WL_noWV, y = plot_BRFu_list[[i]][j,1:length(bright_TOC_mean)], col = j+1, lwd = 2) # add rU lines for each plot 
#       legend("topleft", bty = "n", c(paste("Mean fertility class rU",",","n of plots = ",nrow(plot_BRFu_list[[i]])) # add legend: Number of plots used to calculate mean rU + plot name, number of pixels for each plot 
#               ,paste("plot",rownames(plot_BRFu_list[[i]][1:5,]   ))),# add Plot ID. Only show the first 10 plot IDs
#               lty = c(2,rep(1,j)) ,lwd = 1, col = 1:(j+1)    )  # no legend border, line width, line color.
#       
#     } 
#     
#   } else { # stop the for loop if there is an empty dataframe
#     break
#   }
#   
# }


# 2) #-------------------- One plot with all Mean rU reflectance curves for each fertility class -----------------------------------------------------------------#

par(mfrow=c(1,1)) # set plot settings

# plot first Mean rU 
plot(type = "l", x = WL_noWV, y = lapply(plot_BRFu_list,colMeans)[[1]] ,  main = paste("rU", "slope slf p value test", p_test_slope_slf), #  "brightest %",1 - bright_treshold, "Mean rU AHS for each fertility class"
     xlab = "Wavelength (nm)", ylab = "BRF", lwd =2, ylim = c(0,1))

# plot other mean rU lines
for (i in (1:length(plot_BRFu_list))){
  lines(x = WL_noWV, y = colMeans(plot_BRFu_list[[i]]), col = i , lwd =2)
  legend("topleft",c(paste(1:5,fert_class_list[1:5])), lwd = 2, col = 1:i , bty = "n", lty = 1) #,fert_class_list[i] length(plot_BRFu_list)
}

# add error bars to plot 
error_bar.f(plot_BRFu_list,2) 


# -------------------- plot mean leaf spectra for each fertility class --------------------------------------------------------------------------#

plot(type = "l", x = WL_noWV, y = colMeans(plot_bright_alb_est[[1]]), xlab = "Wavelength (nm)", ylab = "BRF" , lty = 1, lwd = 2) #, 1 - bright_treshold # ,
#main = paste( "Mean leaf albedo for each fertility class", "brightest 10% pixels")

# add PROSPECT reference albedo to plot
lines(x = WL_noWV, y = specref_int_aisa$y, lwd = 2, col = "purple", lty = 6)

# create color index list for lines (skip the 5th color as light blue is not suitable for printing)
col_nr <- c(1,2,3,4,6)

# plot legend and leaf spectra lines for other fertility classes
for ( i in (2:length(plot_bright_alb_est))){
  lines(x = WL_noWV, y = colMeans(plot_bright_alb_est[[i]]), col = col_nr[i], lwd = 2,lty = i )  
  legend("topleft", bty = "n", c(paste( 1:5, fert_class_list[1:5]), "PROSPECT leaf albedo" ), col = c(1:4,6,"purple") , # we skip the 5th color as it is light blue which does not work well with printing.
         lty = c(1:6), x.intersp = 0.5, lwd = 2 ) # add legend, x.int is the space between symbol and legend text, seq_len is the the width of the line symbol.
}

#legend("topleft", bty = "n","  herb-rich",lty = 2, lwd = 2, seg.len = 0.5, x.intersp = 0.5) # first plot a line legend
#legend("topleft", bty = "n", c("","moderately rich","moist upland","dryish upland","dry upland","PROSPECT leaf albedo" ),
#       col = c("white",2:5,"purple"), pch = c(0,2,3,4,5,6)  ) # plot all other data as characters.
# add error bars to plot 
error_bar.f(plot_bright_alb_est,2) 

# -------------------- plot mean canopy BRF for each fertility class --------------------------------------------------------------------------#

plot(type = "l", x = WL_noWV, y = colMeans(canopy_brf_list[[1]]), xlab = "Wavelength (nm)", ylab = "BRF",
     main = "Overstory BRF", ylim = c(0,0.3))
# add canopy BRF for other fertility classes to plot
for ( i in (1:length(canopy_brf_list))){ 
  lines( x = WL_noWV, y = colMeans(canopy_brf_list[[i]]), col = i, lwd = 2 )
  legend("topleft", bty = "n",c(paste(1:5, fert_class_list[1:5])) , col = c(1:i), lwd = 2 ) # add legend
}
error_bar.f(canopy_brf_list,2) # add error bars



# -------------------- plot for each fertility class all plots and their BRF understory spectra --------------------------------------------------------------------------#

#par(mfrow=c(2,3)) # set row,column plot settings

# for ( i in (1:length(plot_bright_alb_est))){ # loop over list of dataframes with each containing BRF understory for a fertility class
# 
#    if (nrow(plot_bright_alb_est[[i]] > 0)){
# 
#      # first, plot the fertility group mean rU
#      plot(type = "l", x = WL_noWV, y = lapply(plot_bright_alb_est,colMeans)[[i]][1:length(bright_TOC_mean)] , main = c(paste("Leaf albedo: fertility class",i,fert_class_list[i], "brightest",1 - bright_treshold,"%" )),
#           ylim = c(0,1), xlab = "Wavelength (nm)", ylab = "Leaf albedo", col = 1, lty = 2, lwd = 2) # plot for each fertility class a new plot
# 
#      for ( j in (1:10)){ # add rU reflectance lines for all plots
# 
#        lines(x = WL_noWV, y = plot_bright_alb_est[[i]][j,1:length(bright_TOC_mean)], col = j+1, lwd = 2) # add rU lines for each plot
#        legend("topleft", bty = "n", c(paste("Mean fertility class Leaf albedo",",","n of plots = ",nrow(plot_bright_alb_est[[i]])) # add legend: Number of plots used to calculate mean rU + plot name, number of pixels for each plot
#                ,paste("plot",rownames(plot_bright_alb_est[[i]][1:5,]   ))),# add Plot ID. Only show the first 10 plot IDs
#                lty = c(2,rep(1,j)) ,lwd = 1, col = 1:(j+1)    )  # no legend border, line width, line color.
# 
#      }
# 
#    } else { # stop the for loop if there is an empty dataframe
#      break
#    }
# 
# }
# 


# -------------------- plot 4 understory measured plots with each: 1) rU field measurement, 2) rU from aisa 3) Mean rU for fert. class. ---------#

par(mfrow=c(2,2)) # set plot settings
lwd.thick = 3

line.wd = 3 # change line width for line charts

# plot fertility class 1
#plot(type = "l",x = WL_noWV, y = colMeans(plot_BRFu_list[[1]]) , lwd = line.wd, xlab = "Wavelength (nm)",
#     ylab = "Ru", main = "rU fertility class 1:Herb rich", col = "blue", lty = 3, ylim = c(0,1)) #plot mean rU for fertility class 1
#legend("topleft", c("Mean rU AHS"), lwd = line.wd, col = "blue", lty = 3, bty = "n")

##### plot fertility class 2: moderately rich upland
plot(type="l",x = field_spectra$wavelength, y = field_spectra$h3 , ylim = c(0,1) , # plot first one line
     main = c("Moderately rich upland"), lwd = line.wd,
     xlab = "Wavelength (nm)", ylab = "BRF")
# plot Aisa derived rU for plot H3
lines(x = WL_noWV, y = colMeans(plot_BRFu_list[[2]][which(rownames(plot_BRFu_list[[2]]) == "H3"),])  , col = "red" ,lwd = line.wd, lty = 2) # add Aisa derived understory BRF
# plot Aisa derived mean rU for this fertility class
lines(x = WL_noWV, y = colMeans(plot_BRFu_list[[2]]) , lwd = line.wd, lty =3, col = "blue")
# add legend
legend("topleft",lty = c(1,2,3),col = c("black","red","blue"), c("rU fieldplot","rU AHS fieldplot", "Mean AHS rU fertility class"),
       bty = "n", lwd = line.wd)

##### plot fert. class 3 Moist upland
plot(type="l",x = field_spectra$wavelength, y = field_spectra$u26 , ylim = c(0,1) , # plot first one line
     main = c("Moist upland"), lwd = line.wd,
     xlab = "Wavelength (nm)", ylab = "BRF")
# plot Aisa derived rU for plot U26
lines(x = WL_noWV, y = colMeans(plot_BRFu_list[[3]][which(rownames(plot_BRFu_list[[3]]) == "U26"),]) , col = "red" ,lwd = line.wd, lty = 2) # add Aisa derived understory BRF
# plot Aisa derived mean rU for this fertility class
lines(x = WL_noWV, y = colMeans(plot_BRFu_list[[3]]) , lwd = line.wd, lty =3, col = "blue")
# add legend
legend("topleft",lty = c(1,2,3),col = c("black","red","blue"), c("rU fieldplot","rU AHS fieldplot", "Mean AHS rU fertility class"),
       bty = "n", lwd = line.wd)


###### plot fert. class 4 Dryish upland
plot(type="l",x = field_spectra$wavelength, y = field_spectra$u18 , ylim = c(0,1) , # plot first one line
     main = c("Dryish upland"), lwd = line.wd,
     xlab = "Wavelength (nm)", ylab = "BRF")
# plot Aisa derived rU for plot U18
lines(x = WL_noWV, y = colMeans(plot_BRFu_list[[4]][which(rownames(plot_BRFu_list[[4]]) == "U18" ),]) , col = "red" ,lwd = line.wd, lty = 2) # add Aisa derived understory BRF
# plot Aisa derived mean rU for this fertility class
lines(x = WL_noWV, y = colMeans(plot_BRFu_list[[4]]) , lwd = line.wd, lty =3, col = "blue")
# add legend
legend("topleft",lty = c(1,2,3),col = c("black","red","blue"), c("rU fieldplot","rU AHS fieldplot", "Mean AHS rU fertility class"),
       bty = "n", lwd = line.wd)


###### plot fert. class type 5 Dry upland
# plot field measured underestory BRF
plot(type="l",x = field_spectra$wavelength, y = field_spectra$u10 , ylim = c(0,1) , # plot first one line
     main = c("Dry upland"), lwd = line.wd,
     xlab = "Wavelength (nm)", ylab = "BRF")
# plot Aisa derived rU for plot U17
lines(x = WL_noWV, y = colMeans(plot_BRFu_list[[which(fert_class_list == "dry upland")]][which(rownames(plot_BRFu_list[[which(fert_class_list == "dry upland")]]) == "U10"),])
      , col = "red" ,lwd = line.wd, lty = 2) # add Aisa derived understory BRF
# plot Aisa derived mean rU for this fertility class
lines(x = WL_noWV, y = colMeans(plot_BRFu_list[[which(fert_class_list == "dry upland")]]) , lwd = line.wd, lty =3, col = "blue")
# add legend
legend("topleft",lty = c(1,2,3),col = c("black","red","blue"), c("rU fieldplot","rU AHS fieldplot", "Mean AHS rU fertility class"),
       bty = "n", lwd = line.wd)


# --------------- plot 4 field measured understory BRF in one plot --------------------------------------------------------------------------#

# set plot parameters
par(mfrow=c(1,1)) # 1 column , 1 row

plot(type="l",x = field_spectra$wavelength, y = field_spectra$h3 , ylim = c(0,0.6), col = 2 , lwd = 2, xlab = "Wavelength (nm)", ylab = "BRF")
lines( x = field_spectra$wavelength, y = field_spectra$u26, col = 3, lwd = 2)
lines( x = field_spectra$wavelength, y = field_spectra$u18, col = 4, lwd = 2)
lines( x = field_spectra$wavelength, y = field_spectra$u10, col = 5, lwd = 2)
legend("topleft", bty = "n", lty = 1, c(fert_class_list[2],fert_class_list[3],fert_class_list[4],fert_class_list[5]), col = c(2,3,4,5), lwd = 2)     





# ----------------- plot sunlit fraction slope as a function of tree species ----------------------------------------------------------- #
# we do this to test if we can use the sunlit fraction slope to seperate tree species. 

# subset dataframes based on trees species purity and plot boxplots for sunlit fraction p slope for each tree species.

# combine dataframe with sunlit fraction slope (p slope) and all fieldplots
plots_df_p <- data.frame(plots,slf_df$slope_slf)

# create subsets for pure plots with tree species fraction above 80% and mixed plots with tree species fraction below 80%.

# all mixed plots
mix_birch.df <- data.frame(subset(plots_df_p, birch_fraction < 0.8))
mix_pine.df <- data.frame(subset(plots_df_p, pine_fraction < 0.8))
mix_spruce.df <- data.frame(subset(plots_df_p, spruce_fraction < 0.8))

# all pure plots
pure_birch.df <- data.frame(subset(plots_df_p,birch_fraction > 0.8 ) )
pure_pine.df <- data.frame(subset(plots_df_p,pine_fraction > 0.8 ) )
pure_spruce.df <- data.frame(subset(plots_df_p,spruce_fraction > 0.8 ) )

# create list of pure, and mix plots.
pure_p_list <- list(pure_pine = pure_pine.df$slf_df.slope_slf,pure_spruce = pure_spruce.df$slf_df.slope_slf,pure_birch = pure_birch.df$slf_df.slope_slf)
mix_p_list <- list(mix_pine = mix_pine.df$slf_df.slope_slf, mix_spruce = mix_spruce.df$slf_df.slope_slf, mix_birch = mix_birch.df$slf_df.slope_slf)

# create boxplots for all plots and for pure and mixed plots
boxplot(pure_p_list,col = rep(c("red","blue","green"), len = length(pure_p_list)), horizontal = T , names = c("pure pine", "pure spruce", "pure birch"),
        main = "p values for pure tree species plots (> 80%)", xlab = "p")
boxplot(mix_p_list, col = rep(c("red", "blue","green"), len = length(mix_p_list)), horizontal = T, names = c("mixed pine", "mixed spruce", "mixed birch"),
        main = "p values for mixed tree species plots (< 80%)", xlab = "p")

boxplot(slf_df$slope_slf, horizontal = T, col = "green", main = "p value distribution all 252 plots. New reference albedo")


# now plot p against LAI for all tree species and fit a line
par(mfrow = c(3,1))

# select font size for legend font
font.size = 2

#pine
plot(x = pure_pine.df$LAI, y = pure_pine.df$slf_df.slope_slf, col = "red", lwd = 2, xlab = "LAI", ylab = "p", main = "p as a function of LAI for pure tree trecies plots (> 80%)")
lm.purepine <- lm(pure_pine.df$slf_df.slope_slf ~ pure_pine.df$LAI)
abline(lm.purepine, col = "black")
r2 <- summary(lm.purepine)$adj.r.squared # store R2 for legend display
legend("topright", legend=bquote(atop(italic(R)^2 == .(format(r2, digits = 2)),"n" == .(nrow(pure_pine.df)))),bty="n", cex = font.size) # add legend with n and R2
legend("topleft", "pure pine", col = "red", pch = 1, bty = "n", cex = font.size)

# spruce
plot(x = pure_spruce.df$LAI, y = pure_spruce.df$slf_df.slope_slf, col = "blue", lwd = 2, xlab = "LAI", ylab = "p")
lm.purespruce <- lm(pure_spruce.df$slf_df.slope_slf ~ pure_spruce.df$LAI)
abline(lm.purespruce, col = "black") # calculate R2
r2 <- summary(lm.purespruce)$adj.r.squared
legend("bottomright", legend=bquote(atop(italic(R)^2 == .(format(r2, digits = 2)),"n" == .(nrow(pure_spruce.df)))),bty="n",cex = font.size) # add legend with n and R2
legend("topleft", "pure spruce", col = "blue", pch = 1, bty = "n",cex = font.size)

#birch
plot(x = pure_birch.df$LAI, y = pure_birch.df$slf_df.slope_slf, col = "green", lwd = 2, xlab = "LAI", ylab = "p")
lm.purebirch <- lm(pure_birch.df$slf_df.slope_slf ~ pure_birch.df$LAI)
abline(lm.purebirch, col = "black")
r2 <- summary(lm.purebirch)$adj.r.squared
legend("bottomright", legend=bquote(atop(italic(R)^2 == .(format(r2, digits = 2)),"n" == .(nrow(pure_birch.df)))),bty="n",cex = font.size) # add legend with n and R2
legend("topleft", "pure birch", col = "green", pch = 1, bty = "n",cex = font.size)

# now we do the same but for mixed plots
#pine
plot(x = mix_pine.df$LAI, y = mix_pine.df$slf_df.slope_slf, col = "red", lwd = 2, xlab = "LAI", ylab = "p", main = "p as a function of LAI for mixed tree trecies plots (< 80%)")
lm.mixpine <- lm(mix_pine.df$slf_df.slope_slf ~ mix_pine.df$LAI)
abline(lm.mixpine, col = "black")
r2 <- summary(lm.mixpine)$adj.r.squared # store R2 for legend display
legend("topright", legend=bquote(atop(italic(R)^2 == .(format(r2, digits = 2)),"n" == .(nrow(mix_pine.df)))),bty="n",cex = font.size) # add n and R2 to plot
legend("topleft", "mixed pine", col = "red", pch = 1, bty = "n",cex = font.size)

# spruce
plot(x = mix_spruce.df$LAI, y = mix_spruce.df$slf_df.slope_slf, col = "blue", lwd = 2, xlab = "LAI", ylab = "p")
lm.mixspruce <- lm(mix_spruce.df$slf_df.slope_slf ~ mix_spruce.df$LAI)
abline(lm.mixspruce, col = "black")
r2 <- summary(lm.mixspruce)$adj.r.squared # store R2 for legend display
legend("topright", legend=bquote(atop(italic(R)^2 == .(format(r2, digits = 2)),"n" == .(nrow(mix_spruce.df)))),bty="n",cex = font.size) # add legend with n and R2
legend("topleft", "mixed spruce", col = "blue", pch = 1, bty = "n",cex = font.size)

#birch
plot(x = mix_birch.df$LAI, y = mix_birch.df$slf_df.slope_slf, col = "green", lwd = 2, xlab = "LAI", ylab = "p")
lm.mixbirch <- lm(mix_birch.df$slf_df.slope_slf ~ mix_birch.df$LAI)
abline(lm.mixbirch, col = "black")
r2 <- summary(lm.mixbirch)$adj.r.squared # store R2 for legend display
legend("topright", legend=bquote(atop(italic(R)^2 == .(format(r2, digits = 2)),"n" == .(nrow(mix_birch.df)))),bty="n",cex = font.size) # add legend with n and R2
legend("topleft", "mixed birch", col = "green", pch = 1, bty = "n",cex = font.size)




# ------------------------------------------ Tets rU seperability among fertility classes ------------------------------------------------- #
 
# objective: To test the separability of the rU reflectance curves for the different fertility classes at selected wavelengths.

# The outcome should be a table with p-values and t-statistics for each unique fertility class combination.  
 
# make sure to first test whether the sample variance  is equal or unequal. This determines what t-test to use. 
# source: https://en.wikipedia.org/wiki/Student%27s_t-test#Independent_two-sample_t-test

# e.g.  Welch's student t.test comparing mean rU of all plots of the first wavelength of one fertility class with all rU for all plots in the second fertility class.

# Set selection of wavelengths for statistical test: 552nm,645nm,682nm,739nm,820)
wvl_sel <- c(35,55,63,75,89)

# find all possible unique combinations
fert_class_combn <- permutations(n = 5, r = 2,v = 1:length(fert_class_list) )# unique combinations for the 5 site fertility rU dataframes 

# create empty dataframe to store t-test results 
p_result.df <- data.frame(row.names = c("p.value","t-statistic","degrees of freedom"))

for ( i_row in 1:nrow(fert_class_combn)){ # loop over row in possible dataframe combination
  
    # select first dataframe in list by using first index position in matrix with all unique possible dataframe combinations
    subset1 <- plot_BRFu_list[[fert_class_combn[i_row,1]]]
    
    # select second dataframe using the second index position ......
    subset2 <- plot_BRFu_list[[fert_class_combn[i_row,2]]] 
  
    temp_p_results <- data.frame() # create empty dataframe for temporary p-values for a specific rU dataframe combination
    
    for (i_wav in wvl_sel){ # loop over selection of wavelengths and perform student's t-test
      
      # subset rU plots for specific wavelength for two dataframes
      subset1_wvl <- subset1[,i_wav] # select rU results from one dataframe for specific wavelength 
      subset2_wvl <- subset2[,i_wav] # ... for second dataframe
      
      # perform Welch's Student t.test on BRF understory plots for selected wavelength
      t.test_result <- t.test(subset1_wvl, subset2_wvl, alternative = "two.sided") 
      temp_p_results <- rbind(temp_p_results,c("p.value" = round(t.test_result$p.value,10), "t-statistic" = round(t.test_result$statistic,2),
                                               "parameter" = round(t.test_result$parameter)) )
      
      # add names of dataframe combination as rowname 
      rownumber = nrow(temp_p_results)
      rownames(temp_p_results)[rownumber] <- paste(fert_class_combn[i_row,1], "vs",fert_class_combn[i_row,2], WL_noWV[i_wav],"nm") 
      # colnames(temp_p_results) <- c("p.value", "t-statistic") # add column names
      
      }
  
    temp_p_results_transposed <- t(temp_p_results) # transpose p-value dataframe    
    
  p_result.df <- cbind(p_result.df, temp_p_results_transposed)
}

# print all significant p-values below 0.05 for all possible fertility class combinations and five selected wavelengths.
p_rU_subset <- p_result.df[1,which(p_result.df[1,] < 0.05)]

# store p values for each wavelength in seperate variable
p_ru_552nm <- p_rU_subset[,grepl("552", names(p_rU_subset) )]
p_ru_645nm <- p_rU_subset[,grepl("645", names(p_rU_subset) )]
p_ru_683nm <- p_rU_subset[,grepl("683", names(p_rU_subset) )]
p_ru_740nm <- p_rU_subset[,grepl("740", names(p_rU_subset) )]
p_ru_820nm <- p_rU_subset[,grepl("820", names(p_rU_subset) )]


# create table and print to html file
setwd("C:/Users/vmvincent/PhD/data/understory_reflectance_project/") # set working directory to store p-value table 
p_table_html <- xtable(p_result.df[1,which(p_result.df[1,] <= 0.05)])
print.xtable(p_table_html, type="html", file="p_values_table_BRFu.html")


#### testing 

# calculate t-test for each wavelength and each unique fertility combination
t.test.f <- function(dataframe_result){
  
  # select all rU plots in one wavelength of one dataframe
  
  # selecct all plots in one wavelength of second dataframe
  
  # perform t-test on both subsets and return p-value and t-statistic.
  
  t.test_temp <- apply(dataframe_result, 2 ,t.test , alternative = "two.sided" ) # apply t.test function on each plot rU as a function of wavelength  
  #names(t.test_temp) <- apply()
  results_t.test <- data.frame(sapply(t.test_temp, "[", c("statistic", "p.value") )) #store p-value and t-statistic in a results dataframe for each wavelength
  #rownames(results_t.test) <- c("t-statistic", "p-value") # add rownames
}

# store t-test results for each dataframe
results_t.test_BRFu_fert1 <- t.test.f(plot_BRFu_list[[1]])
results_t.test_BRFu_fert2 <- t.test.f(plot_BRFu_list[[2]])
results_t.test_BRFu_fert3 <- t.test.f(plot_BRFu_list[[3]])
results_t.test_BRFu_fert4 <- t.test.f(plot_BRFu_list[[4]])
results_t.test_BRFu_fert5 <- t.test.f(plot_BRFu_list[[5]])


# print p-values for each dataframe at selected wavelenghts
p_table <-xtable(results_t.test_BRFu_fert1[2,wvl_sel])
print.xtable(p_table, type="html", file="p_values_table.html")


# create table with number of plots for each fertility class
BRFu_plots_table <- as.table(c(length(plots_fert_1_index) , length(plots_fert_2_index) ,length(plots_fert_3_index), length(plots_fert_4_index), length(plots_fert_5_index ) ) )
names(BRFu_plots_table) <- fert_class_list[1:5]


# add understory pixel values to spatial dataframe