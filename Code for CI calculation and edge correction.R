#  *****************************************************************
#  ***                                                           ***
#  ***            Calculate Competition indices (CI)             *** 
#  ***          Written by Jingning Shi on Oct 20th 2021         ***
#  ***          Revised by Jingning Shi on Aug 27th 2022         ***
#  ***                                                           ***
#  *****************************************************************

# The latest version can be found at https://github.com/JingningShi

library(magrittr)
library(spatstat)

### Part 1 Calculate distance-dependent CI

## The function calc_CI is used to calculate distance-dependent CIs. Search
#  radius (SR) is used to identify competitors around the subject tree.
#  I here wrote a total of 18 SR, including 5 horizontal fixed SR (CI1-CI5), 
#  4 horizontal variable SR (CI6-CI9) and 9 vertical angle-based SR (CI10-18).
#  Principle and methodological details of CI1-5 can be found in Hegyi(1974) 
#  and Lorimer(1983), and CI6-CI18 can be found in Sharma et al.(2016).

## Parameter df is a data frame containing DBH, species, tree height, crown
#  width, plot size and XY coordinates of individual trees. Note that ALL 
#  trees in df should be located within the same plot. If df contains trees 
#  from different plots, this function can produce incorrect results.

## Parameter r is the horizontal fixed search radius of CI1-CI5. 
#  For example, r = 10 means that only trees located within a 10-m radius of
#  the subject tree are considered as competitor.

## Parameter var is a character string indicating which tree size variable 
#  (DBH or tree height) is to be used for the CI calculation. Its value must 
#  be "DBH" or "HT" (Sharma et al. 2016). Defaults to "DBH", which is more 
#  commonly used in forest growth and yield modeling.

## Parameter edge is a logical value indicating whether edge effects correction 
#  should be performed during the computation. Note that edge correction only 
#  applies to rectangular plot. See details below.

## Details of all cited sources in the code file can be found in the reference
#  list of the manuscript.

calc_CI <- function(df, r, var = 'DBH', edge = FALSE){
  
  set.seed(1020)
  # Check whether two trees have the same XY coordinate.
  tem_chr <- paste0(df$X, '-', df$Y)
  if(length(unique(tem_chr)) != nrow(df)){
    # The function runif is primarily to avoid coordinate clashes because two 
    # trees having the same XY coordinates will lead to incorrect results. 
    x <- df$X + runif(nrow(df), -0.01, 0.01) 
    y <- df$Y + runif(nrow(df), -0.01, 0.01) 
  } else {
    x <- df$X
    y <- df$Y
  }
  
  d <- df$D # D = DBH
  ht <- df$ht # ht = tree height
  
  if(!var %in% c('DBH', 'HT'))
    stop('Incorrect tree size variable')
  
  # Which tree size variable to use ?
  if(var == 'DBH'){
    v <- d
  } else {
    v <- ht
  }
  
  # ot = overstory tree. The purpose of this line of code is to obtain a
  # variable required for CI3 calculation (Lorimer, 1983). Here I select 
  # three tallest trees to represent overstory trees within a plot. You 
  # can opt for other methods depending on your datasets.
  ot_df <- df[order(df$ht, decreasing = T)[1:3],]
  
  # R is a coefficient of CI3 and R = b * mean crown radius of overstory trees. 
  # b is a constant. As suggested by Lorimer(1983), b is set to 3.5 
  # CWe is the crown width in the west side, and so on.
  R <- 3.5 * mean((ot_df$CWe + ot_df$CWw + ot_df$CWs + ot_df$CWn)/4)
  if(is.nan(R)){
    R <- 1
    warning(paste("The CI3 calculation produced incorrect results", 
                  "due to the lack of crown width information"))
  }
  
  # Calculate the minimum height angle (mha) of each tree, see Part C below.
  mha_df <- matrix(nrow = 0, ncol = 9) %>% as.data.frame
  for (i in 1:length(x)) {
    ha_dist <- ht[i]/tan(seq(30, 70, by = 5)/180*pi)
    mha_df <- rbind(mha_df, ha_dist)
  }
  colnames(mha_df) <- paste0('dist', seq(30, 70, by = 5))

  eval(parse(text = paste0('CI', 1:18, '_list = vector()', collapse = "; ")))
  for(i in 1:length(x)){
    # x1 and y1 are the x and y coordinates of the subject tree in the 
    # current cycle, respectively.
    x1 <- x[i]
    y1 <- y[i] 
    dist <- sqrt((x - x1)^2 + (y - y1)^2) 
    
    ## Part A: Calculate CI1-CI5 using the fixed search radius r
    
    # Check whether a competitor is located within the search radius. The 
    # not equal operator != 0 is used to exclude the subject tree itself.
    j <- which(dist <= r & dist != 0) 
    
    # For the previous line of code, if XY coordinates of a competitor is 
    # exactly the same with the subject tree, this competitor is ignored
    # due to the not equal operator. So I use the runif function to ensure 
    # that the position of all trees will not be exactly the same. 
    
    if(edge){
      # Calculate edge correction coefficient for CI1-CI5. 
      # Part 2 below describes the correct_edge function in detail.
      edge_coef <- correct_edge(df[i,], r)
      ec <- 1/edge_coef
    } else {
      ec <- 1
    }
    
	  if(length(j) == 0){
	    # Competitor doesn't exist
      ci1 <- 0
      ci2 <- 0
      ci3 <- 0
      ci4 <- 0
      ci5 <- 0
    } else {
		  ci1 <- sum((v[j]/v[i])/(dist[j]+1)*ec) # Hegyi (1974) index
		  ci2 <- sum((v[j]/v[i])^2/dist[j]*ec)   # Rouvinen and Kuuluvainen (1997)
		  ci3 <- sum((v[j]/v[i])/(dist[j]/R)*ec) # Lorimer (1983) index
		  ci4 <- sum(v[j]/dist[j]*ec)            # Rouvinen and Kuuluvainen (1997)
		  ci5 <- sum((v[j]/v[i])/dist[j]^2*ec)   # Rouvinen and Kuuluvainen (1997)
		}
	  CI1_list <- c(CI1_list, ci1)
	  CI2_list <- c(CI2_list, ci2)
	  CI3_list <- c(CI3_list, ci3)
	  CI4_list <- c(CI4_list, ci4)
	  CI5_list <- c(CI5_list, ci5)
	  
	  ## Part B: Calculate CI6-CI9 using the variable search radius r
	  # See Sharma et al. (2016) for details.
	  
	  # CI6-7
	  sr6 <- d[i]*0.25
	  sr7 <- d[i]*0.33
	  j6 <- which(dist <= sr6 & dist != 0) 
	  j7 <- which(dist <= sr7 & dist != 0) 
	  
	  if(edge){
	    edge_coef6 <- correct_edge(df[i,], sr6)
	    ec6 <- 1/edge_coef6
	    edge_coef7 <- correct_edge(df[i,], sr7)
	    ec7 <- 1/edge_coef7
	  } else {
	    ec6 <- 1
	    ec7 <- 1
	  }
	  
	  if(length(j6) == 0){
	    ci6 <- 0
	  } else {
	    ci6 <- sum((v[j6]/v[i])/(dist[j6]+1)*ec6) # Eq(1) in Sharma et al. (2016)
	  }
	  if(length(j7) == 0){
	    ci7 <- 0
	  } else {
	    ci7 <- sum((v[j7]/v[i])/(dist[j7]+1)*ec7) # Eq(1) in Sharma et al. (2016)
	  }
	  CI6_list <- c(CI6_list, ci6)
	  CI7_list <- c(CI7_list, ci7)
	  
	  # CI8-CI9 See Table 1 in Sharma et al. (2016)
	  sr8 <- (d[i] + d)/6
	  sr9 <- (d[i] + d)/8
	  j8 <- which(dist <= sr8 & dist != 0) 
	  j9 <- which(dist <= sr9 & dist != 0) 

	  if(length(j8) == 0){
	    ci8 <- 0
	  } else {
	    if(edge){
	      # For CI8-CI18, SR is set to the distance from subject tree to 
	      # outermost competitor.
	      edge_coef8 <- correct_edge(df[i,], max(dist[j8]))
	      ec8 <- 1/edge_coef8
	    } else {
	      ec8 <- rep(1, times = length(j8))
	    }
	    ci8 <- sum((v[j8]/v[i])/(dist[j8]+1)*ec8) # Eq(1) in Sharma et al. (2016)
	  }
	  if(length(j9) == 0){
	    ci9 <- 0
	  } else {
	    if(edge){
	      edge_coef9 <- correct_edge(df[i,], max(dist[j9]))
	      ec9 <- 1/edge_coef9
	    } else {
	      ec9 <- rep(1, times = length(j9))
	    }
	    ci9 <- sum((v[j9]/v[i])/(dist[j9]+1)*ec9) # Eq(1) in Sharma et al. (2016)
	  }
	  CI8_list <- c(CI8_list, ci8)
	  CI9_list <- c(CI9_list, ci9)
	  
	  ## Part C: Calculate CI10-CI18 using the vertical angle-based SR
	  
	  # Competitors are identified using height angle starting from the base of
	  # the subject tree, instead of from 1 m above the base. The height angle of 
	  # each potential competitor is calculated as the vertical angle from the 
	  # base of the subject tree to the top-center of the crown of a competitor. 
	  # Only neighbor trees with an height angle greater than the minimum height 
	  # angle θ are selected as competitors. See Richards et al. (2008) for a 
	  # graphical illustration of vertical angle-based SR. 
	  
	  # Here θ = 30°, 35°, 40°, 45°, 50°, 55°, 60°, 65°, 70°
	  
	  # Define a function to identify competitors of a subject tree at a given 
	  # height angle and calculate CI10-CI18 for that subject tree. 
	  f_mha <- function(df, v, i, dist, edge, mha_df, angle){
	    sr <- mha_df[[paste0('dist', angle)]]
	    j <- which(dist <= sr & dist != 0) 
	    if(length(j) == 0){
	      ci <- 0
	    } else {
	      if(edge){
	        edge_coef <- correct_edge(df[i,], max(dist[j]))
	        ec <- 1/edge_coef
	      } else {
	        ec <- 1
	      }
	      ci <- sum((v[j]/v[i])/(dist[j]+1)*ec) # Eq(1) in Sharma et al. (2016)
	    }
	    ci
	  }
	  
	  ci10 <- f_mha(df, v, i, dist, edge, mha_df, 30)
	  ci11 <- f_mha(df, v, i, dist, edge, mha_df, 35)
	  ci12 <- f_mha(df, v, i, dist, edge, mha_df, 40)
	  ci13 <- f_mha(df, v, i, dist, edge, mha_df, 45)
	  ci14 <- f_mha(df, v, i, dist, edge, mha_df, 50)
	  ci15 <- f_mha(df, v, i, dist, edge, mha_df, 55)
	  ci16 <- f_mha(df, v, i, dist, edge, mha_df, 60)
	  ci17 <- f_mha(df, v, i, dist, edge, mha_df, 65)
	  ci18 <- f_mha(df, v, i, dist, edge, mha_df, 70)
	  
	  CI10_list <- c(CI10_list, ci10)
	  CI11_list <- c(CI11_list, ci11)
	  CI12_list <- c(CI12_list, ci12)
	  CI13_list <- c(CI13_list, ci13)
	  CI14_list <- c(CI14_list, ci14)
	  CI15_list <- c(CI15_list, ci15)
	  CI16_list <- c(CI16_list, ci16)
	  CI17_list <- c(CI17_list, ci17)
	  CI18_list <- c(CI18_list, ci18)
  }

  # Combine CIs by columns 
  CI_df <- cbind(
    CI1_list, CI2_list, CI3_list, CI4_list, CI5_list,
    CI6_list, CI7_list, CI8_list, CI9_list, CI10_list, 
    CI11_list, CI12_list, CI13_list, CI14_list, 
    CI15_list, CI16_list, CI17_list, CI18_list) %>% as.data.frame
  colnames(CI_df) <- paste0('CI', 1:18)
  
  # Combine CI_df and df together  
  df <- cbind(df, CI_df)
}

### Part 2 Correct edge effects caused by off-plot competitors

# Off-plot competitor trees can affect subject trees near the plot boundary.
# To correct such edge effects, all competition indices were weighted using 
# an area-weighted method (Das et al. 2008). Now suppose r has the value 10m. 
# For trees that were located within 10 m of the plot edge, index values were 
# divided by the proportion of a 10-m radius circle centered on the tree that 
# would lie inside the plot boundaries. For example, if only 80% of the 10-m 
# radius circle centered around a given tree would be within the plot, the raw 
# index value for that tree was divided by 0.80. R package spatstat was used 
# to compute the area of the overlap between circle and rectangle.

# Run this function may take a few seconds because the parameter eps of 
# overlap.owin is set to 0.05. Such a setting is more time-consuming but 
# can significantly increase the computation accuracy of overlapping areas. 

# Parameter r is a numeric vector indicating search radius and the length of r
# can be 1 or larger. 

correct_edge <- function(df, r){
  # Position of the subject tree
  x1 <- df$X
  y1 <- df$Y
  # Xrange is the length of the rectangular plot.
  xrange <- df$Xrange
  yrange <- df$Yrange
  
  # Create an owin object representing the rectangle in two-dimensional plane.
  plot_window <- owin(c(0, xrange), c(0, yrange))
  
  # Calculate the edge length of inscribed square of circle with radius r.
  # This square is used to generate a circular owin object with radius r.
  edge_len <- r*sqrt(2)
  
  oa_vector <- vector() # oa = overlapping area
  for (i in 1:length(r)) {
    # Check whether edge correction is required. This aims to avoid unnecessary 
    # overlapping area calculation.  
    c1 <- x1 - r[i] >= 0
    c2 <- x1 + r[i] <= xrange
    c3 <- y1 - r[i] >= 0
    c4 <- y1 + r[i] <= yrange
    if(c1 & c2 & c3 & c4){
      oa_vector <- c(oa_vector, 1)
    } else {
      len1 <- edge_len[i]/2
      tree_window <- owin(c(x1 - len1, x1 + len1),
                          c(y1 - len1, y1 + len1)) %>% 
        boundingcircle.owin(eps = 0.05)
      oa <- overlap.owin(plot_window, tree_window) # get the overlapping area
      # Calculate the proportion of overlap
      oa_vector <- c(oa_vector, round(oa/(r[i]^2*pi), 2))
    }
  }
  oa_vector
}

### Part 3 An example of CI calculation

# Now let's see how to use this function. Don't forget to set working directory
# before reading the presentation data.
df <- read.csv('Presentation data for CI.csv') 

# Calculate CI without considering the edge effect correction.
df2 <- calc_CI(df, r = 10, var = 'DBH', edge = FALSE)
View(df2) # See what happened

# Edge effects correction is a very time-consuming calculation. Take the
# presentation data as an example. If you set edge to TRUE, running the
# code below took almost half a minute on my computer.
df3 <- calc_CI(df, r = 10, var = 'DBH', edge = TRUE)

