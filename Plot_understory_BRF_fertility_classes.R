
#----------------------- model to calculate understory BRF as a function of fertility class. ---------------------------------#

#----------- import data -----------------#
# source('C:/Users/vmvincent/OneDrive/PhD/R_analysis/understory_reflectance/Extract_bright_canopy_pixels_for_can_albedo.R')  
#source("c:/Users/vmvincent/OneDrive/PhD/R_analysis/understory_reflectance/import_aisa_data_for_understory_brf.R", echo = FALSE)
source("C:/git_projects/PHD_understory_paper/estimate_canopy_cover_Aisa_beta_regression.R")
#source("C:/Users/VincentMarkiet/Documents/Test_rscript/", echo = FALSE)
source('C:/Users/VincentMarkiet/OneDrive - Advian Oy/Tiedostot/Personal/Education/phd/backup_rscript/understory_reflectance/functions/error_bars_function.R', echo=FALSE)
#source("C:/Users/VincentMarkiet/OneDrive - Advian Oy/Tiedostot/Personal/Education/phd/backup_rscript/understory_reflectance/import_aisa_data_for_understory_brf.R", echo = FALSE)
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
#plots_fert_1_index <- which( plots$Kasvupaikka == 1 ) #&  c(plots$pine_fraction > 0.9 | plots$spruce_fraction > 0.9 | plots$birch_fraction > 0.9))  # excluded class 1 due to low number of plots (2)
#plots_fert_2_index <- which( plots$Kasvupaikka == 2 ) #&  c(plots$pine_fraction > 0.9 | plots$spruce_fraction > 0.9 | plots$birch_fraction > 0.9)) 
#plots_fert_3_index <- which( plots$Kasvupaikka == 3 ) #&  c(plots$pine_fraction > 0.9 | plots$spruce_fraction > 0.9 | plots$birch_fraction > 0.9)) 
#plots_fert_4_index <- which( plots$Kasvupaikka == 4 ) #&  c(plots$pine_fraction > 0.9 | plots$spruce_fraction > 0.9 | plots$birch_fraction > 0.9)) 
#plots_fert_5_index <- which( plots$Kasvupaikka == 5 ) #&  c(plots$pine_fraction > 0.9 | plots$spruce_fraction > 0.9 | plots$birch_fraction > 0.9)) 
# plots_fert_6_index <- which( plots$Kasvupaikka == 6 ) #&  c(plots$pine_fraction > 0.9 | plots$spruce_fraction > 0.9 | plots$birch_fraction > 0.9)) 


#------temporary to check understory for extreme LAI values ----------#
library("dplyr")

lai = 1.8 #4.9 # 1.8
plots_fert_2_index = which( plots$Kasvupaikka == 2 & plots$LAI < lai ) # selected the plots with smallest .25 procentile of LAI 
plots_fert_3_index = which( plots$Kasvupaikka == 3 & plots$LAI < lai )
plots_fert_4_index = which(plots$Kasvupaikka == 4 & plots$LAI < lai )
plots_fert_5_index = which(plots$Kasvupaikka == 5 & plots$LAI < lai )

# ---------------------------------------------------------------------  #

# create list with subsetplot indices 
indexlist <- list(  plots_fert_2_index,  plots_fert_3_index,  plots_fert_4_index, plots_fert_5_index ) #plots_fert_1_index, ,  plots_fert_6_index

##### 2) calculate BRF understory for subset plots

# create an empty list for results
plot_BRFu_list <- list()

# create empty list for bright pixel sample size tables
canopy_brf_list <- list()

# create list for albedo results
plot_bright_alb_est <- list()

# create list to store PRI for each plot by fertility class
pri_ru_fertility_list <- list()

# create empty list to store dataframes with canopy PRI for each fertility class
pri_canopy_fertility_list <- list()

# create list for leaf spectra PRI
pri_leaf_spectra_list <- list()
ndvi_ru_list <- list()
ndvi_rc_list <- list()
ndvi_can_list <- list()
red_edge_ru_list <- list()
red_edge_rc_list <- list()
red_edge_can_list <- list()
cl_green_ru_list <- list()
cl_green_rc_list     <- list()
cl_green_can_list  <- list()
sr_ru_list <-   list()
sr_rc_list <-   list()
sr_can_list <-  list()
car_ru_list <-  list()
car_rc_list <-  list()
car_can_list <- list()


# create list of fertility class names
fert_class_list <- c("moderately rich","moist upland","dryish upland","dry upland") # removed class ,"(1) "herb-rich" and (6) "nutrient poor upland"  due to lack of field plot data for these fertility classes

# set canopy gap fraction (CGF) input. Either from field measurements (LAI-2000) or predicted canopy gap fraction using Hadi's Betaregression between BRF in red and measured canopy gap fraction.

# predicted CGF
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


# set brightest pixel bright_treshold. e.g. 0.8 = Use 20% brightest pixels
bright_treshold <- 0.9

# create dataframe to store p, LAI, tree species fraction values for each plot 
p_alb_results <- data.frame()


# source script with vegetation indices: red edge, ndvi, pri, red edge chlorophyll
source("C:/git_projects/PHD_understory_paper/HSI_vegetation_indices.R")

for ( i_index in 1:length(indexlist)){ 
  
  i_fert <- indexlist[[i_index]] # loop over list of plot indices 
  
  # create empty dataframe for temporary results
  plot_result <- data.frame() # create empty dataframe
  bright_pix_table <- data.frame() # create empty dataframe
  albedo_est_plots <- data.frame() # create empty dataframe for calculated leaf albedos
  canopy_brf_plot <- data.frame() # ....... for calculated canopy BRF
  
  pri_ru.df <- data.frame() # empty list for understory pri
  pri_can.df <- data.frame() # empty list for canopy pri
  pri_leaf.df <- data.frame()
  
  ndvi_ru.df <- data.frame()
  ndvi_rc.df <- data.frame()
  ndvi_can.df <- data.frame()
  red_edge_ru.df <- data.frame()
  red_edge_rc.df <- data.frame()
  red_edge_can.df <- data.frame()
  cl_green_ru.df <- data.frame()
  cl_green_rc.df <- data.frame()
  cl_green_can.df <- data.frame()
  sr_ru.df <-       data.frame()
  sr_rc.df <-       data.frame()
  sr_can.df <-      data.frame()
  car_ru.df <-      data.frame()
  car_rc.df <-      data.frame()
  car_can.df <-     data.frame()
  
  
  
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
    # retrieve p and rho using prospect corrected leaf albedo
    
    yy <- unlist(bright_TOC_mean[69:82] / albedo_eff[69:82] ) #) #rc[69:82]     # estimate ratio of brf from leaf albedo using wavelength 710-790 nm
    
    
    bk <- lm(yy ~ unlist(bright_TOC_mean[69:82] ) ) # fit ordinary linear model between ratio of TOC BRF and leaf albedo.
    coefs_slf <- coef(bk) # linear model coeffecients
    slope_slf <- coefs_slf[2] # slope / p
    intercept_slf <- coefs_slf[1] # intercept / rho
    aa <- slope_slf * bright_TOC_mean + intercept_slf # Eq. 5 in manuscript Markiet,Mottus .2019
    rc <- as.data.frame(t(bright_TOC_mean  / aa )) # store estimated leaf albedo in dataframe and tranpose rows to columns.

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
    
    
    # ----------- Vegetation indices -------------------------------------------------- # 
    
    
    # 5) calculate understory, leaf and canopy PRI for each plot
    # use PRI equation (531nm - 570nm) / (531nm + 570nm)
    pri.ru <- (ru[30] - ru[39]) / (ru[30] + ru[39])
    pri.can <- (BRF_can[30] - BRF_can[39]) / (BRF_can[30] + BRF_can[39]) 
    pri.leaf <- (rc[30] - rc[39]) / (rc[30] + rc[39]) 
    
    # 6) calculate understory, leaf, and canopy NDVI 
    ndvi_ru <- (ru[97]- ru[64]) / (ru[97] + ru[64]) # NDVI= (R858-R687)/(R858+R687)
    ndvi_rc <- (rc[97]- rc[64]) / (rc[97] + rc[64])
    ndvi_can <- (BRF_can[97] - BRF_can[64]) / (BRF_can[97] + BRF_can[64])
    
    # 7) calculate Chlorophyll index (CI) 
    cl_green_ru <- cl_green.f(brf_tree_type = ru)
    cl_green_rc  <- cl_green.f(brf_tree_type = rc)
    cl_green_can  <- cl_green.f(brf_tree_type = BRF_can)
    
    # calculate red edge
    red_edge_ru <- r_edge.f(brf_input = ru)
    red_edge_rc <- r_edge.f(brf_input = rc)
    red_edge_can <- r_edge.f(brf_input = BRF_can)
    
    # calculate simple ratio
    sr_ru <- sr.f(brf_input = ru)
    sr_rc <- sr.f(brf_input = rc)
    sr_can <- sr.f(brf_input = BRF_can)
    
    # calculate carotenoid reflectance index 
    car_ru <- carotenoid.f(brf_input = ru)
    car_rc <- carotenoid.f(brf_input = rc)
    car_can <- carotenoid.f(brf_input = BRF_can)
    
    # 7)  store results
    
    ######## temp #############
    
    
    #p_alb_results <- rbind(p_alb_results,"p" = c(round(slope_slf,2),"LAI" = plots$LAI[i_plot], "pine_fraction" = plots$pine_fraction[i_plot],
    #                       "spruce_fraction" = plots$spruce_fraction[i_plot], "birch_fraction" = plots$birch_fraction[i_plot]))
    #rownumber = dim(plot_result)[1]
    #rownames(p_alb_results)[rownumber] <- paste(brf$fieldplot_buffers.ID[i_plot] )
    
    #########################
    
    # 5.1) store results, one row per plot
    if (dim(plot_result)[1] == 0 ) { # insert check to test if there are more than 0 rows.
      plot_result <- ru
      albedo_est_plots <- rc
      canopy_brf_plot <- BRF_can
      pri_ru.df <- pri.ru
      pri_can.df <- pri.can
      pri_leaf.df <- pri.leaf
      ndvi_ru.df <- ndvi_ru
      ndvi_rc.df <- ndvi_rc
      ndvi_can.df <- ndvi_can
      red_edge_ru.df <- red_edge_ru
      red_edge_rc.df <- red_edge_rc
      red_edge_can.df <- red_edge_can
      cl_green_ru.df <- cl_green_ru
      cl_green_rc.df  <- cl_green_rc
      cl_green_can.df  <- cl_green_can
      sr_ru.df        <- sr_ru
      sr_rc.df        <- sr_rc
      sr_can.df       <- sr_can
      car_ru.df       <- car_ru
      car_rc.df       <- car_rc
      car_can.df      <- car_can
        
      
    } else { # Add rows to existing result dataframe if number of rows is higher than 0.
      #temp_data <- data.frame(ru,Nr_bright_pix) # create temporary vector to store BRFus and Nr of bright pixels used for each plot. 
      # colnames(temp_data) <- colnames(plot_result) # assign result dataframe column names to temporary dataframe so they can be merged later.
      plot_result <- rbind( plot_result,ru) # add temporary dataframe to result dataframe
      albedo_est_plots <- rbind(albedo_est_plots,rc)
      canopy_brf_plot <- rbind(canopy_brf_plot, BRF_can)
      pri_ru.df <- rbind(pri_ru.df, pri.ru)
      pri_can.df <- rbind(pri_can.df, pri.can)
      pri_leaf.df <- rbind(pri_leaf.df, pri.leaf)
      ndvi_ru.df <- rbind(ndvi_ru.df, ndvi_ru)
      ndvi_rc.df <- rbind(ndvi_rc.df, ndvi_rc)
      ndvi_can.df <- rbind(ndvi_can.df,ndvi_can)
      red_edge_ru.df <- rbind(red_edge_ru.df,red_edge_ru)
      red_edge_rc.df <- rbind(red_edge_rc.df,red_edge_rc)
      red_edge_can.df <- rbind(red_edge_can.df,red_edge_can)
      cl_green_ru.df <- rbind(cl_green_ru.df,cl_green_ru)
      cl_green_rc.df  <- rbind(cl_green_rc.df, cl_green_rc)
      cl_green_can.df <- rbind(cl_green_can.df, cl_green_can)
      sr_ru.df        <- rbind(sr_ru.df, sr_ru)
      sr_rc.df        <- rbind(sr_rc.df,sr_rc)
      sr_can.df       <- rbind(sr_can.df,sr_can)
      car_ru.df       <- rbind(car_ru.df, car_ru)
      car_rc.df       <- rbind(car_rc.df,car_rc)
      car_can.df      <- rbind(car_can.df,car_can)
      
      
    }
    
    # 5.3) store plot IDs as row number
    rownumber = dim(plot_result)[1]
    rownames(plot_result)[rownumber] <-  paste(brf$fieldplot_buffers.ID[i_plot] )
    #
    rownumber = dim(albedo_est_plots)[1]
    rownames(albedo_est_plots)[rownumber] <-  paste(brf$fieldplot_buffers.ID[i_plot] )

    rownumber = dim(canopy_brf_plot)[1]
    rownames(canopy_brf_plot)[rownumber] <- paste(brf$fieldplot_buffers.ID[i_plot])

    rownumber = dim(pri_ru.df)[1]
    rownames(pri_ru.df)[rownumber] <- paste(brf$fieldplot_buffers.ID[i_plot])

    rownumber = dim(pri_can.df)[1]
    rownames(pri_can.df)[rownumber] <- paste(brf$fieldplot_buffers.ID[i_plot])

    rownumber = dim(pri_leaf.df)[1]
    rownames(pri_leaf.df)[rownumber] <- paste(brf$fieldplot_buffers.ID[i_plot])
    
    rownumber = dim(ndvi_ru.df)[1]
    rownames(ndvi_ru.df)[rownumber] <- paste(brf$fieldplot_buffers.ID[i_plot])
    
    rownumber = dim(ndvi_rc.df)[1]
    rownames(ndvi_rc.df)[rownumber] <- paste(brf$fieldplot_buffers.ID[i_plot])
    
    rownumber = dim(ndvi_can.df)[1]
    rownames(ndvi_can.df)[rownumber] <- paste(brf$fieldplot_buffers.ID[i_plot])

    rownumber = dim(red_edge_ru.df)[1]
    rownames(red_edge_ru.df)[rownumber] <- paste(brf$fieldplot_buffers.ID[i_plot])
    
    rownumber = dim(red_edge_rc.df)[1]
    rownames(red_edge_rc.df)[rownumber] <- paste(brf$fieldplot_buffers.ID[i_plot])

    rownumber = dim(cl_green_rc.df)[1]
    rownames(cl_green_rc.df)[rownumber] <- paste(brf$fieldplot_buffers.ID[i_plot])
    
    rownumber = dim(cl_green_can.df)[1]
    rownames(cl_green_can.df)[rownumber] <- paste(brf$fieldplot_buffers.ID[i_plot])
    
    rownumber = dim(sr_ru.df)[1]
    rownames(sr_ru.df)[rownumber] <- paste(brf$fieldplot_buffers.ID[i_plot])
    
    rownumber = dim(sr_rc.df)[1]
    rownames(sr_rc.df)[rownumber] <- paste(brf$fieldplot_buffers.ID[i_plot])
    
    rownumber = dim(sr_can.df)[1]
    rownames(sr_can.df)[rownumber] <- paste(brf$fieldplot_buffers.ID[i_plot])
    
    
    rownumber = dim(car_ru.df)[1]
    rownames(car_ru.df)[rownumber] <- paste(brf$fieldplot_buffers.ID[i_plot])
    
    rownumber = dim(car_rc.df)[1]
    rownames(car_rc.df)[rownumber] <- paste(brf$fieldplot_buffers.ID[i_plot])
    
    rownumber = dim(car_can.df)[1]
    rownames(car_can.df)[rownumber] <- paste(brf$fieldplot_buffers.ID[i_plot])
    
  }
  
  # 5.4) store dataframe results in a list 
  plot_BRFu_list[[i_index]] <- plot_result # store results in list.
  canopy_brf_list[[i_index]] <- bright_pix_table # store dataframe with number of bright pixels in a list.
  plot_bright_alb_est[[i_index]] <- albedo_est_plots # store leaf albedo dataframe in list
  canopy_brf_list[[i_index]] <- canopy_brf_plot # add each result dataframe to list
  pri_ru_fertility_list[[i_index]] <- pri_ru.df # store understory pri 
  pri_canopy_fertility_list[[i_index]] <- pri_can.df # store canopy pri
  pri_leaf_spectra_list[[i_index]] <- pri_leaf.df # store leaf pri
  ndvi_ru_list[[i_index]] <-  ndvi_ru.df
  ndvi_rc_list[[i_index]] <- ndvi_rc.df 
  ndvi_can_list[[i_index]] <- ndvi_can.df 
  red_edge_ru_list[[i_index]] <- red_edge_ru.df 
  red_edge_rc_list[[i_index]] <- red_edge_rc.df 
  red_edge_can_list[[i_index]] <- red_edge_can.df
  cl_green_ru_list[[i_index]] <- cl_green_ru.df 
  cl_green_rc_list[[i_index]] <- cl_green_rc.df 
  cl_green_can_list[[i_index]] <- cl_green_can.df
  sr_ru_list[[i_index]] <- sr_ru.df       
  sr_rc_list[[i_index]] <- sr_rc.df       
  sr_can_list[[i_index]] <- sr_can.df      
  car_ru_list[[i_index]] <- car_ru.df      
  car_rc_list[[i_index]] <- car_rc.df      
  car_can_list[[i_index]] <- car_can.df     
  
  
  # assign fertility class names to each column 
  #colnames(canopy_brf_list[[i_index]]) <- fert_class_list[i_index] 
  colnames(pri_ru_fertility_list[[i_index]]) <- fert_class_list[i_index] 
  colnames(pri_canopy_fertility_list[[i_index]]) <- fert_class_list[i_index] 
  colnames(pri_leaf_spectra_list[[i_index]]) <- fert_class_list[i_index] 
  colnames(ndvi_ru_list[[i_index]]) <- fert_class_list[i_index]
  colnames(ndvi_rc_list[[i_index]]) <- fert_class_list[i_index]
  colnames(ndvi_can_list[[i_index]]) <-fert_class_list[i_index]
  colnames(red_edge_ru_list[[i_index]]) <- fert_class_list[i_index]
  colnames(red_edge_rc_list[[i_index]]) <- fert_class_list[i_index]
  colnames(red_edge_can_list[[i_index]]) <-fert_class_list[i_index]
  colnames(cl_green_ru_list[[i_index]]) <- fert_class_list[i_index]
  colnames(cl_green_rc_list[[i_index]]) <- fert_class_list[i_index]
  colnames(cl_green_can_list[[i_index]]) <- fert_class_list[i_index]
  colnames(sr_ru_list[[i_index]]) <- fert_class_list[i_index]
  colnames(sr_rc_list[[i_index]]) <- fert_class_list[i_index]
  colnames(sr_can_list[[i_index]]) <- fert_class_list[i_index]
  colnames(car_ru_list[[i_index]] )<- fert_class_list[i_index]
  colnames(car_rc_list[[i_index]]) <- fert_class_list[i_index]
  colnames(car_can_list[[i_index]]) <- fert_class_list[i_index]
  
  
}








#-------------------------- plot results ------------------------------------------------------------------------------------------------------#

# set plot settings
#par(mfrow = c(1,2))

# plot linegraph with BRF understory for each plot as a function of wavelength.
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



# 1) ------------------------- create PRI boxplot of all plots for each fertility class -------------------------------------------------------------------

# 1.1 understory PRI
# melt PRI dataframes into  dataframe sorted by fertility.
pri_ru_df_melted <- melt(pri_ru_fertility_list) # merge all understory PRI dataframes into one dataframe

# 1.2 overstory PRI
# melt dataframe
pri_can_df_melted <- melt(pri_canopy_fertility_list)

#boxplot of PRI for different fertility classes
qplot(factor(variable), value, data = pri_can_df_melted, geom = "boxplot", xlab = "fertility class", ylab = "Overstory PRI", color = variable)

# 1.3 Aisa TOC PRI
# use PRI equation (531nm - 570nm) / (531nm + 570nm)
pri_aisa <- (brf[,30] - brf[,39]) / (brf[,30] + brf[,39])


# --------------------------- Scatter plots and t-tests of TOC PRI against understory, overstory and leaf PRI. ---------------------------------------------------------- #

par(mfrow=c(2,2))

# 1.4.1 scatter plot of TOC aisa PRI and modeled overstory pri 
pri_can_lm <- lm(pri_can_df_melted$value ~ pri_aisa )
col_list <- c("black","red","blue","green","orange")
plot(pri_aisa, pri_can_df_melted$value, xlab = "TOC pri", ylab = "overstory pri", col = "black", cex = 1, cex.axis = 1 ) #text(pri_aisa,pri_can_df_melted$value,labels = pri_can_df_melted$L1)

#legend("topleft", unique(pri_can_df_melted$variable), col = col_list, pch = 1 )
abline(pri_can_lm, col = "red") # fit a line 
r2 <- summary(pri_can_lm)$adj.r.squared # store R2 for legend display
legend("bottomleft", legend=bquote(atop(italic(R)^2 == .(format(r2, digits = 2)),"n" == .(nrow(pri_can_df_melted)))),bty="n",cex = 1)


# 1.4.2) plot leaf pri vs TOC aisa pri
pri_leaf_df_melted <- melt(pri_leaf_spectra_list) # merge all leaf PRI dataframes into one dataframe
pri_leaf_lm <- lm(pri_leaf_df_melted$value ~ pri_aisa )
plot(pri_aisa, pri_leaf_df_melted$value, xlab = "TOC pri", ylab = expression(paste("r"[C]," " ,"pri")), col = "black", cex = 1 ) 
abline(pri_leaf_lm, col = "red")
r2 <- summary(pri_leaf_lm)$adj.r.squared # store R2 for legend display
legend("bottomleft", legend=bquote(atop(italic(R)^2 == .(format(r2, digits = 2)))),bty="n",cex = 1)



# 1.4.3) plot understory pri vs TOC aisa pri
pri_ru_lm <- lm(pri_ru_df_melted$value ~ pri_aisa )
plot(pri_aisa, pri_ru_df_melted$value, xlab = "TOC pri", ylab = expression(paste("r"[italic(F)]," ","pri")), col = "black",  cex = 1 ) 
abline(pri_ru_lm, col = "red")
r2 <- summary(pri_leaf_lm)$adj.r.squared # store R2 for legend display
legend("bottomleft", legend=bquote(atop(italic(R)^2 == .(format(r2, digits = 2)))),bty="n",cex = 1)


# perform t-test between TOC PRI and overestory, understory, leaf PRI to test for significant differences ------------------------- 

# leaf PRI against TOC PRI
pri_leaf.t.test <- t.test(pri_leaf_df_melted$value, pri_aisa, "two.sided")

# overstory PRI against TOC PRI
pri_canopy.t.test <- t.test(pri_can_df_melted$value, pri_aisa, "two.sided")

# understory PRI against TOC PRI
pri_understory.t.test <- t.test(pri_ru_df_melted$value, pri_aisa, "two.sided")



# 2) -------------------- One plot with all Mean rU reflectance curves for each fertility class -----------------------------------------------------------------#

par(mfrow=c(1,1)) # set plot settings

# plot first Mean rU 
plot(type = "l", x = WL_noWV, y = lapply(plot_BRFu_list,colMeans)[[1]] ,  main = paste("rU", "slope slf p value test", p_test_slope_slf), #  "brightest %",1 - bright_treshold, "Mean rU AHS for each fertility class"
     xlab = "Wavelength (nm)", ylab = "BRF", lwd =2, ylim = c(0,0.7))

# plot other mean rU lines
for (i in (1:length(plot_BRFu_list))){
  lines(x = WL_noWV, y = colMeans(plot_BRFu_list[[i]]), col = i , lwd =2)
  legend("topleft",c(paste(1:length(fert_class_list),fert_class_list[1:length(fert_class_list)])), lwd = 2, col = 1:i , bty = "n", lty = 1) #,fert_class_list[i] length(plot_BRFu_list)
}

# add error bars to plot 
error_bar.f(plot_BRFu_list,2) 



# -------------------- plot mean leaf spectra for each fertility class --------------------------------------------------------------------------#

plot(type = "l", x = WL_noWV, y = colMeans(plot_bright_alb_est[[1]]), xlab = "Wavelength (nm)", ylab = "BRF" , lty = 1, lwd = 2,ylim=c(0,1)) #, 1 - bright_treshold # ,
#main = paste( "Mean leaf albedo for each fertility class", "brightest 10% pixels")

# add PROSPECT reference albedo to plot
lines(x = WL_noWV, y = albedo_eff , lwd = 2, col = "purple", lty = 6) # specref_int_aisa$y use this if you want to plot PROSPECT model

# create color index list for lines (skip the 5th color as light blue is not suitable for printing)
col_nr <- c(1,2,3,4,6)

# plot legend and leaf spectra lines for other fertility classes
for ( i in (2:length(plot_bright_alb_est))){
  lines(x = WL_noWV, y = colMeans(plot_bright_alb_est[[i]]), col = col_nr[i], lwd = 2,lty = i )  
  legend("topleft", bty = "n", c(paste(1:length(fert_class_list),fert_class_list[1:length(fert_class_list)]), "PROSPECT leaf albedo" ), col = c(1:4,6,"purple") , # we skip the 5th color as it is light blue which does not work well with printing.
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
  legend("topleft", bty = "n",c(paste(1:length(fert_class_list),fert_class_list[1:length(fert_class_list)])) , col = c(1:i), lwd = 2 ) # add legend
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
lines(x = WL_noWV, y = colMeans(plot_BRFu_list[[1]][which(rownames(plot_BRFu_list[[1]]) == "H3"),])  , col = "red" ,lwd = line.wd, lty = 2) # add Aisa derived understory BRF
# plot Aisa derived mean rU for this fertility class
lines(x = WL_noWV, y = colMeans(plot_BRFu_list[[1]]) , lwd = line.wd, lty =3, col = "blue")
# add legend
legend("topleft",lty = c(1,2,3),col = c("black","red","blue"), c("rU fieldplot","rU AHS fieldplot", "Mean AHS rU fertility class"),
       bty = "n", lwd = line.wd)

##### plot fert. class 3 Moist upland
plot(type="l",x = field_spectra$wavelength, y = field_spectra$u26 , ylim = c(0,1) , # plot first one line
     main = c("Moist upland"), lwd = line.wd,
     xlab = "Wavelength (nm)", ylab = "BRF")
# plot Aisa derived rU for plot U26
lines(x = WL_noWV, y = colMeans(plot_BRFu_list[[2]][which(rownames(plot_BRFu_list[[2]]) == "U26"),]) , col = "red" ,lwd = line.wd, lty = 2) # add Aisa derived understory BRF
# plot Aisa derived mean rU for this fertility class
lines(x = WL_noWV, y = colMeans(plot_BRFu_list[[2]]) , lwd = line.wd, lty =3, col = "blue")
# add legend
legend("topleft",lty = c(1,2,3),col = c("black","red","blue"), c("rU fieldplot","rU AHS fieldplot", "Mean AHS rU fertility class"),
       bty = "n", lwd = line.wd)


###### plot fert. class 4 Dryish upland
plot(type="l",x = field_spectra$wavelength, y = field_spectra$u18 , ylim = c(0,1) , # plot first one line
     main = c("Dryish upland"), lwd = line.wd,
     xlab = "Wavelength (nm)", ylab = "BRF")
# plot Aisa derived rU for plot U18
lines(x = WL_noWV, y = colMeans(plot_BRFu_list[[3]][which(rownames(plot_BRFu_list[[3]]) == "U18" ),]) , col = "red" ,lwd = line.wd, lty = 2) # add Aisa derived understory BRF
# plot Aisa derived mean rU for this fertility class
lines(x = WL_noWV, y = colMeans(plot_BRFu_list[[3]]) , lwd = line.wd, lty =3, col = "blue")
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
legend("topleft", bty = "n", lty = 1, c(fert_class_list[1],fert_class_list[2],fert_class_list[3],fert_class_list[4]), col = c(2,3,4,5), lwd = 2)     





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
 
# make sure to first test whether the sample vari ance  is equal or unequal. This determines what t-test to use. 
# source: https://en.wikipedia.org/wiki/Student%27s_t-test#Independent_two-sample_t-test

# e.g.  Welch's student t.test comparing mean rU of all plots of the first wavelength of one fertility class with all rU for all plots in the second fertility class.

# Set selection of wavelengths for statistical test: 552nm,645nm,682nm,739nm,820)
wvl_sel <- c(35,55,63,75,89)

# find all possible unique combinations
fert_class_combn <- permutations(n = length(fert_class_list), r = 2,v = 1:length(fert_class_list) )# unique combinations for the 5 site fertility rU dataframes 

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
      t.test_result <- t.test(subset1_wvl, subset2_wvl, alternative = "two.sided") # The t.test function performs a Welch t-test by default if var.equal = TRUE is used then the pooled variance is used to estimate the variance.
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
#setwd("C:/Users/vmvincent/PhD/data/understory_reflectance_project/") # set working directory to store p-value table 
setwd("C:/Users/VincentMarkiet/OneDrive - Advian Oy/Tiedostot/Personal/Education/phd/article_understory_reflectance/manuscript_plots/")
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
#results_t.test_BRFu_fert5 <- t.test.f(plot_BRFu_list[[5]])


# print p-values for each dataframe at selected wavelenghts
p_table <-xtable(results_t.test_BRFu_fert1[2,wvl_sel])
print.xtable(p_table, type="html", file="p_values_table.html")


# create table with number of plots for each fertility class
BRFu_plots_table <- as.table(c(length(plots_fert_2_index) ,length(plots_fert_3_index), length(plots_fert_4_index), length(plots_fert_5_index ) ) )
names(BRFu_plots_table) <- fert_class_list[1:length(fert_class_list)]











# plot ndvi scatter plots for all fertility classes 

# melt results into dataframe and assign new column names
ndvi_can_melt <- melt(ndvi_can_list)
colnames(ndvi_can_melt) <- c("fertility","ndvi", "code")

p_ndvi_can <- ggplot(data = ndvi_can_melt) +
        geom_boxplot(mapping = aes(x = fertility, y = ndvi, fill = fertility )) + coord_flip() + ggtitle("overstory ndvi") + scale_fill_manual( values = c("black","red","green","blue"))

ndvi_plotter <- function(ndvi_df, title, legend_p){
  
  # first melt all results into one dataframe
  ndvi_melt.df = melt(ndvi_df)
  # assign new column names
  colnames(ndvi_melt.df) <- c("fertility","ndvi", "code")
  # build boxplot framework
  plot_p = ggplot(data = ndvi_melt.df) +
                geom_boxplot(mapping = aes(x = fertility, y = ndvi, fill = fertility, shape = fertility )) +
                     ggtitle(paste0("ndvi") + scale_fill_manual( values = c("black","red","green","blue")))+
                    ylim(0,1) + ggtitle(title) + theme(legend.position = legend_p)
  
  # now plot boxploot
  plot(plot_p)  
}

# add all ndvi results to a list

p_ndvi_can <- ndvi_plotter(ndvi_can_list, title = "overstory", legend_p = "none")
p_ndvi_ru <- ndvi_plotter(ndvi_ru_list, title = expression(italic("r")[italic(F)~",X"]), legend_p = "none")
p_ndvi_rc <- ndvi_plotter(ndvi_rc_list, title =  expression("r"[C]), legend_p = "right")


# pri scatterplots using ggplot

# merge all pri results and TOC pri together into one dataframe

data.frame("pri_toc" = (brf[,30] - brf[,39]) / (brf[,30] + brf[,39]), "id" = brf$fieldplot_buffers.ID,  )

plot(pri_ru_df_melted$value,pri_aisa) 
pri_can_df_melted 
pri_leaf_df_melted


multiplot(p_ndvi_can,p_ndvi_ru,p_ndvi_rc, cols = 3 )


#qplot(factor(variable), value, data = pri_can_df_melted, geom = "boxplot", xlab = "fertility class", ylab = "Overstory PRI", color = variable)

# create boxplot of LAI for each fertility class
plots_c = plots # copy plots dataframe
# recode all fertility classes to fertility class names
plots_c = plots_c %>% mutate(fertility = recode(Kasvupaikka, '1' = 'herb-rich','2' = 'moderately-rich', '3' = 'moist-upland', '4' = 'dryish-upland', '5' = 'dry-upland', '6' = 'nutriend-poor-upland' ))