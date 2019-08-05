###Map code
#install.packages('udunits2', type = 'source', repo = 'cran.rstudio.com')
#install.packages("sf")
#devtools::install_github("tidyverse/ggplot2")
#devtools::install_github("r-lib/rlang", build_vignettes = TRUE)



# plot output from analysis_diff_in_diff

library(data.table)
#library(ggplot2)
library(viridis)
#install.packages("rgdal", repos = "http://cran.us.r-project.org", type = "source")
require(ggplot2)
library(sf)
library(scales)
#library(hyspdisp)

rm(list = ls())

#======================================================================#
# define function to read data and make plots
#======================================================================#
map_metrics <- function(zip_dataset_sf, 
                        metric, 
                        #year.map,
                        legend_lims = NULL,
                        longitude_LL = -180,
                        plot_xlim = c(-123, -69),
                        plot_ylim = c(24, 50),
                        outpath,
                        filename = NULL,
                        plot.title = NULL,
                        legend.text.angle = 0,
                        legend.title.text = NULL,
                        legend.title.size = 10,
                        diff.as.frac = NULL,
                        facility.loc = NULL,
                        show.legend = NULL){
  z_d_sf <- copy( zip_dataset_sf)
  s = c(metric, 'ZIP', 'geometry')
  z_d_sf <- na.omit(z_d_sf[Longitude >= longitude_LL,
                            s,
                            with = F])
  
#  if ( length( year.map) > 1){
#    color.option <- 'magma'
    
#    if( is.null( filename))
#      filename <- ifelse( is.null( diff.as.frac),
#                          paste('map_diff_',metric, '_',  
#                                paste( year.map, collapse = '-'), 
 #                               '.png', sep = ''),
 #                         paste('map_diff_frac',metric, '_',  
 #                               paste( year.map, collapse = '-'), 
 #                               '.png', sep = ''))
    
#    plot_sf_orig <- z_d_sf
#    plot_sf_targ <- z_d_sf
#    setnames(plot_sf_orig, metric, 'm1')
#    setnames(plot_sf_targ, metric, 'm2')
    
 #   plot_sf_merg <- data.table( merge( plot_sf_orig, 
#                                       plot_sf_targ, 
#                                       'ZIP'))
#    if( is.null( diff.as.frac)){
#      plot_sf_merg[ , metric := m1 - m2]
#    } else
#      plot_sf_merg[ , metric := (m1 - m2) / m1]
    
#    setnames(plot_sf_merg, 'geometry.x', 'geometry')
#  }else {
    color.option <- 'viridis'
     if (is.null( plot.title)) plot.title = paste(metric, 'in')
   if( is.null( filename))
      filename = paste('map_', metric,'.png', sep = '')
    
    plot_sf_merg <- z_d_sf
    setnames(plot_sf_merg, metric, 'metric')
#  }
  
  if( is.null(legend_lims))
    legend_lims <- c(0, quantile(plot_sf_merg[ , metric], 
                                 .95, 
                                 na.rm = T))
  if( !is.null( facility.loc)) {
    point_loc <- geom_point(
      x = facility.loc[1],
      y = facility.loc[2],
      shape = 1,
      colour = "forestgreen",
      # fill = "darkred",
      size = .75
    )
  } else
    point_loc <- NULL
  
  if( 'discrete' %in% legend_lims){
    colorscale <-
      scale_color_viridis(
        name = legend.title.text,
        discrete = T,
        option = color.option,
        direction = 1,
        na.value="white",
        guide = guide_legend( title.position = 'top',
                              title.hjust = 0.5,
                              title.vjust = 0 ,
                              label.vjust = 1)
      )
    fillscale <-
      scale_fill_viridis(
        name = legend.title.text,
        discrete = T,
        option = color.option,
        direction = 1,
        na.value = NA,
        guide = guide_legend( title.position = 'top',
                              title.hjust = 0.5,
                              title.vjust = 0 ,
                              label.vjust = 1)
      )
  }  else{
    colorscale <-
      scale_color_viridis(
        name = legend.title.text,
        discrete = F,
        option = color.option,
        limits = legend_lims,
        oob = squish,
        direction = 1,
        na.value = NA,
        guide = guide_colorbar( title.position = 'top',
                                title.hjust = 0.5,
                                title.vjust = 0 ,
                                label.vjust = 1)
      )
    fillscale <-
      scale_fill_viridis(
        name = legend.title.text,
        discrete = F,
        option = color.option,
        limits = legend_lims,
        oob = squish,
        direction = 1,
        na.value="white",
        guide = guide_colorbar( title.position = 'top',
                                title.hjust = 0.5,
                                title.vjust = 0 ,
                                label.vjust = 1)
      )
  }
  
  gg <- ggplot(data = plot_sf_merg, 
               aes(fill  = metric, 
                   color = metric)) +
    theme_bw() + 
    labs(title = plot.title) +
    geom_polygon(
      data = map_data("state"),
      aes(x = long, y = lat, group = group),
      fill = "white",
      colour = "white",
      size = .25
    ) +
    geom_sf(size = 0.01) +
    geom_polygon(
      data = map_data("state"),
      aes(x = long, y = lat, group = group),
      fill = NA,
      colour = "grey50",
      size = .25
    ) +
    coord_sf(
      xlim = plot_xlim,
      ylim = plot_ylim,
      datum = NA
    ) +
    point_loc +
    colorscale +
    fillscale +
    theme(
      plot.title = if( !is.null(plot.title)) {
        element_text(size = 16, hjust = 0.5)
      } else
        element_blank(),
      axis.title = element_text(size = 24),
      axis.text = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.title = if( !is.null(legend.title.text)) {
        element_text(size = legend.title.size)
      } else
        element_blank(),
      legend.title.align = 1,
      legend.position = if( !is.null(show.legend)) {
        show.legend
      } else
        c(.20, .12),  
      legend.text = element_text(size = 8, 
                                 angle = legend.text.angle),
      legend.background = element_rect(fill = 'transparent'),
      legend.key.size = unit(.05, 'npc'),
      legend.direction = 'horizontal',
      rect = element_blank() #( fill = 'transparent')
      # panel.background = element_rect( fill = 'transparent',
      #                                  color = 'white'),
      # plot.background = element_rect( fill = 'transparent',
      #                                 color = 'white')
    )
  
  invisible(ggsave(
    file.path(outpath, filename),
    gg,
    width = 13.5,
    height = 7.79,
    unit = 'cm',
    bg = "transparent",
    device="jpeg"
  ))
  
}


#======================================================================#
# load data
#======================================================================#
file_path = "/Users/shirleyliao/Dropbox/Shirley's Dissertation/Paper2/Figures/Maps/"
zcta_shapefile <- paste(file_path,"cb_2017_us_zcta510_500k.shp",sep="")
crosswalk_csv <- paste(file_path,"Zip_to_ZCTA_crosswalk_2015_JSI.csv",sep="") 
load(paste(file_path,"subset.data.frame.R",sep=""))
load(paste(file_path,"cut72_data.R",sep=""))

zip_dataset <- data.table(cbind(final.data,cut72.data[[4]]))

cw <- fread(crosswalk_csv)

cw$ZCTA <- formatC( cw$ZCTA,
                    width = 5,
                    format = "d",
                    flag = "0") # to merge on zcta ID
library(sf)
zips <- st_read(zcta_shapefile)
setnames( zips, 'ZCTA5CE10', 'ZCTA')
zips <- merge( zips, cw, by = "ZCTA", all = F, allow.cartesian = TRUE) # all.x = TRUE, all.y = FALSE, allow.cartesian = TRUE)
zips$ZIP <- formatC( zips$ZIP,
                     width = 5,
                     format = "d",
                     flag = "0") # to merge on zcta ID

zd_sf <- data.table( merge( zips, zip_dataset, by = c('ZIP'), all.y = T))

indi_cols = zd_sf[,79:82]
names(indi_cols) = c("Winter_expose","Spring_expose","Summer_expose","Fall_expose")
zd_sf = cbind(zd_sf,apply(indi_cols,2,as.character))

load(paste(file_path,"final_app_membership.R",sep=""))
ppta.once = apply(apply(save.membership,c(1,2),prod),1,max)

load(paste(file_path,"final_app_ipw_weights.R",sep=""))

load(paste(file_path,"final_app_stable_weights.R",sep=""))

load(paste(file_path,"final_app_overlap_weights.R",sep=""))
overlap.weights = apply(save.overlap,1,prod)

current_names = names(zd_sf)
zd_sf = cbind(zd_sf,as.character(ppta.once),log(ipw.weights),log(stabilized.weights),overlap.weights)
names(zd_sf) = c(current_names,"PPTA","log.IPW","log.Stabilized","Overlap")
#zd_sf$PPTA = as.numeric(zd_sf$PPTA)

#======================================================================#
# make plots
#======================================================================#
#outpath   <- file_path
#outpath.eval = file_path
#======================================================================#
# plots
#======================================================================#
###############################################################################
# FIGURE 1
map_metrics( zip_dataset_sf = zd_sf,
             metric = 'Winter_expose',
             #year.map = 2012,
             legend_lims = 'discrete',
             plot.title = 'Exposure in Winter 2012',
             legend.title.text = 'Exposure status',
             outpath = file_path,
             legend.text.angle = 0,
             filename = 'map_exposure_winter.png')

map_metrics( zip_dataset_sf = zd_sf,
             metric = 'Spring_expose',
             #year.map = 2012,
             legend_lims = 'discrete',
             plot.title = 'Exposure in Spring 2012',
             legend.title.text = 'Exposure status',
             outpath = file_path,
             legend.text.angle = 0,
             filename = 'map_exposure_spring.png')

map_metrics( zip_dataset_sf = zd_sf,
             metric = 'Summer_expose',
             #year.map = 2012,
             legend_lims = 'discrete',
             plot.title = 'Exposure in Summer 2012',
             legend.title.text = 'Exposure status',
             outpath = file_path,
             legend.text.angle = 0,
             filename = 'map_exposure_summer.png')

map_metrics( zip_dataset_sf = zd_sf,
             metric = 'Fall_expose',
             #year.map = 2012,
             legend_lims = 'discrete',
             plot.title = 'Exposure in Fall 2012',
             legend.title.text = 'Exposure status',
             outpath = file_path,
             legend.text.angle = 0,
             filename = 'map_exposure_fall.png')

#########################################################################
## FIGURES 2-5
library(viridis)
map_metrics( zip_dataset_sf = zd_sf,
             metric = 'PPTA',
             #year.map = 2012,
             legend_lims = 'discrete',
             plot.title = 'COP as identified by PPTA',
             legend.title.text = 'Contribute to COP',
             outpath = file_path,
             legend.text.angle = 0,
             filename = 'map_ppta.png')

map_metrics( zip_dataset_sf = zd_sf,
             metric = 'log.IPW',
             #year.map = 2012,
             #legend_lims = c( 0, max( zd_sf$log.IPW)) ,
             legend_lims = NULL,
             plot.title = 'Regional patterns of Logged IPW',
             legend.title.text = '2.5 to 97.5 quantile of Log(IPW)',
             outpath = file_path,
             legend.text.angle = 30,
             filename = 'map_ipw.png')

map_metrics( zip_dataset_sf = zd_sf,
             metric = 'log.Stabilized',
             #year.map = 2012,
             legend_lims = NULL,
             #legend_lims = c( 0, max( zd_sf$log.Stabilized)) ,
             plot.title = 'Regional patterns in Logged SW',
             legend.title.text = '2.5 to 97.5 quantile of Log(SW)',
             outpath = file_path,
             legend.text.angle = 30,
             filename = 'map_stable.png')

map_metrics( zip_dataset_sf = zd_sf,
             metric = 'Overlap',
             #year.map = 2012,
             #legend_lims = c( 0, max( zd_sf$Overlap)) ,
             legend_lims = NULL,
             plot.title = 'Regional patterns in OW',
             legend.title.text = '2.5 to 97.5 quantile of OW',
             outpath = file_path,
             legend.text.angle = 30,
             filename = 'map_overlap.png')

##############################################################################
