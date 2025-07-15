### Integration of the Global Water and Lake Sectors within the ISIMIP framework through scaling of streamflow inputs to lakes
# Ana I. Ayala. Limnology Unit, Department of Ecology and Genetics, Uppsala University, Uppsala, Sweden.
# José L. Hinostroza. Faculty of Civil Engineering, National University of Engineering, Lima, Perú.
# July 2025 

# libraries
library(data.table)
library(sf)
library(ncdf4)
library(dplyr)
library(ggplot2)
library(readxl)
library(writexl)
library(hydroGOF)
library(yarrr)

# working directory
setwd("Scaling of streamflow inputs to lakes")

# read functions written in another script
source("helper_functions_script.R")

### input data
## isimip3, global lake sector, representative lakes
# https://github.com/icra/ISIMIP_Lake_Sector
# one selected lake per pixel at a 0.5º resolution: the selected lake has a depth corresponding to the depth weighted median (weighted by area of the lakes) 
# for all the lakes contained in each pixel with a 0.5º resolution
lakes <- st_read(file.path("./inputs","isimip3_global_lake_sector","HL_selected_cent.shp")) # load shapefile 

## isimip 3a,  global water sector, watergap2-2e
# https://data.isimip.org/search/tree/ISIMIP3a/OutputData/water_global/watergap2-2e/gswp3-w5e5/obsclim/
# watergap2, gswp3-w5e5, obsclim_histsoc, default
# x[lon, lat, time], x: qtotal, qg, dis
# total runoff (qtotal)
nc_qtot <- nc_open(file.path("inputs","isimip3a_watergap2-2e","watergap2-2e_gswp3-w5e5_obsclim_histsoc_default_qtot_global_monthly_1901_2019.nc")) # open a netCDF file
names(nc_qtot$var) # var name
qtot <- ncvar_get(nc_qtot,"qtot") # total runoff (kg/m2/s)
qtot <- qtot/1000 # kg/m2/s -> m/s, kg/m2/s * 1l/1kg * 1m3/1000l 
lon <- ncvar_get(nc_qtot,"lon") # longitude (degree) 
lat <- ncvar_get(nc_qtot,"lat") # latitude (degree)
time <- ncvar_get(nc_qtot, "time") # time (days), start on 1901-01-01
date <- seq(as.Date("1901-01-01"), by="month", length.out=dim(time))
nc_close(nc_qtot) # close a netCDF file
# groundwater runoff (qg)
nc_qg <- nc_open(file.path("inputs","isimip3a_watergap2-2e","watergap2-2e_gswp3-w5e5_obsclim_histsoc_default_qg_global_monthly_1901_2019.nc")) # open a netCDF file
names(nc_qg$var) # var name
qg <- ncvar_get(nc_qg,"qg") # groundwater runoff (kg/m2/s)
qg <- qg/1000 # kg/m2/s -> m/s, kg/m2/s * 1l/1kg * 1m3/1000l 
# lon <- ncvar_get(nc_qg,"lon") # longitude (degree) 
# lat <- ncvar_get(nc_qg,"lat") # latitude (degree)
# time <- ncvar_get(nc_qg, "time") # time (days), start on 1901-01-01
nc_close(nc_qg) # close a netCDF file
# discharge (dis)
nc_dis <- nc_open(file.path("inputs","isimip3a_watergap2-2e","watergap2-2e_gswp3-w5e5_obsclim_histsoc_default_dis_global_monthly_1901_2019.nc")) # open a netCDF file
names(nc_dis$var) # var name
dis <- ncvar_get(nc_dis,"dis") # discharge (m3/s)
# lon <- ncvar_get(nc_dis,"lon") # longitude (degree) 
# lat <- ncvar_get(nc_dis,"lat") # latitude (degree)
# time <- ncvar_get(nc_dis, "time") # time (days), start on 1901-01-01
nc_close(nc_dis) # close a netCDF file

## isimip 3a, river routing 
# flow direction (fdir)
# fdir[lon, lat]
nc_fdir <- nc_open(file.path("inputs","isimip3a_flow_direction","ddm30_flowdir_cru_neva.nc")) # open a netCDF file
names(nc_fdir$var) # var name
fdir <- ncvar_get(nc_fdir,"flowdirection") # flow direction (1: east, 2: south east, 3: south, 4: south west,
#                                                            5: west, 6: north west; 7: north; 8: north east,
#                                                            0: sink or outlet to the ocean)
lon_fdir <- ncvar_get(nc_fdir,"lon") # longitude (degree)
lat_fdir <- ncvar_get(nc_fdir,"lat") # latitude (degree)
nc_close(nc_fdir) # close a netCDF file
# slope (slope)
# slopes[lon, lat]
nc_slope <- nc_open(file.path("inputs","isimip3a_slope","ddm30_slopes_cru_neva.nc")) # open a netCDF file
names(nc_slope$var) # var name
slope <- ncvar_get(nc_slope,"slope") # slope: m/m (elevation/horizontal distance)
# lon_slope <- ncvar_get(nc_slope,"lon") # longitude (degree)
# lat_slope <- ncvar_get(nc_slope,"lat") # latitude (degree)
nc_close(nc_slope) # close a netCDF file

### approach I(a/b): lake area <= grid area
## validation against hype model outputs: https://hypeweb.smhi.se/explore-water/historical-data/europe-time-series

# selected lakes (listed), lag (TRUE/FALSE) and method (slope/all)
df_selected_lakes <- data.frame(lake=c(102, 1052, 1057, 1061, 1068, 1072, 1079, 1081, 1092, 1097,
                                       1104, 1150, 1165, 11564, 11566, 11616, 11693, 11697, 11734, 11749,
                                       11794, 11842, 11864, 11877, 11973, 12000, 12058, 12247, 12316, 12423, 
                                       12533, 12623, 12681, 12727, 12729, 12791, 12809, 12936, 12965, 13020,  
                                       13030, 13057, 13120, 13126, 13175, 137598, 140992, 141844, 142240, 142660,  
                                       143213, 145192, 147907, 148104, 149126, 149804, 151109, 151471, 152117, 152977,
                                       153185, 157766, 158251, 158481, 159827, 160157, 1235053, 1287757),
                                lag=c(FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, 
                                      TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
                                      TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, 
                                      TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, 
                                      FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, 
                                      FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE,
                                      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE))
# method for selection of the upstream grids
method <- "slope" # slope: grid with steepest slope, all: all grids at the same level
dir.create("outputs_hype")
dir.create(file.path("./outputs_hype", method))
outputs_hype <- data.frame() # outputs
for (l in df_selected_lakes$lake) {
  # lag between simulations and observations: TRUE, without lag between simulations and observations: FALSE
  lag <- df_selected_lakes[df_selected_lakes$lake==l, ]$lag # TRUE or FALSE
  
  # create a directory for each lake
  dir.create(file.path("./outputs_hype", method, l))
  
  # grids occupied by the lake
  # hydrolakes + centroid inside polygon (lon,lat)
  lake <- lakes[lakes$Hylak_id==l, ]
  # grids occupied by the lake (0.5degree x 0.5degree)
  grids <- generate_grids(lake)
  
  # find the intersection between lake and grids and divide it into parts 
  int_grids_lake <- st_intersection(grids, lake)
  # centroid inside polygon (lat, long) of each part
  int_grids_lake_cent <- st_coordinates(st_point_on_surface(int_grids_lake)) 
  colnames(int_grids_lake_cent) <- c("lon", "lat")
  # nearest grids
  lake_grids <- data.frame()
  for (i in 1:nrow(int_grids_lake_cent)) {
    lake_grids <- rbind(lake_grids, t(c(nearest_f1(int_grids_lake_cent[i, "lon"])[1], nearest_f1(int_grids_lake_cent[i, "lat"])[1]))) # nearest grid (grid it belongs to; lon, lat)
  }
  
  # upstream grids of the lake
  upstream_grids <- upstream_grids_of_a_lake(fdir, lon_fdir, lat_fdir, lake_grids)
  
  # grid area where the centroid is located (km2)
  grid_area <- areakm2lat(nearest_f1(lake$Cent_Lat)[1])
  # watershed area (km2)
  wshd_area <- lake$Wshd_area
  # number of grids to pick up: ratio between watershed area and  grid area
  ratio_area <- wshd_area/grid_area
  ratio_area <- ifelse(ratio_area<1, 1, floor(ratio_area))
  
  # qin [m3/s] = qtot+qg [m/s] * wshd_area [m2]
  if (ratio_area<=nrow(lake_grids)) { # upstream grids <= grids occupied by the lake
    # qtot, qg, qin for each grid 
    qin_grids_lake <- list()
    for (i in 1:nrow(lake_grids)) {
      qin_grid_lake <- data.frame(date, qtot[which(lon==lake_grids[i, "lon"]), which(lat==lake_grids[i, "lat"]), 1:dim(qtot)[3]], qg[which(lon==lake_grids[i, "lon"]), which(lat==lake_grids[i, "lat"]), 1:dim(qg)[3]])
      setnames(qin_grid_lake, c("date","qtot","qg"))
      qin <- (qin_grid_lake$qtot+qin_grid_lake$qg)*lake$Wshd_area*1e+06 # qtot+qg [m/s] * wshd_area [m2] -> qin [m3/s]
      qin_grid_lake <- data.frame(qin_grid_lake, qin)
      qin_grids_lake <- append(qin_grids_lake, list(qin_grid_lake))
    }
    selected_grids <- lake_grids
    # average qin
    qin_sim <- 0
    for (i in 1:length(qin_grids_lake)) {
      qin_sim <- qin_sim+(qin_grids_lake[[i]]$qin/length(qin_grids_lake))
    }
    qin_sim <- data.frame(date, qin_sim)
    
  } else if (ratio_area>nrow(lake_grids)) { # upstream grids > grids occupied by the lake
    
    if(nrow(upstream_grids)<ratio_area) { # upstream grids < ratio (select all)
      # qtot, qg, qin for each grid 
      qin_grids_lake <- list()
      for (i in 1:nrow(upstream_grids)) {
        qin_grid_lake <- data.frame(date, qtot[which(lon==upstream_grids[i, "lon"]), which(lat==upstream_grids[i, "lat"]), 1:dim(qtot)[3]], qg[which(lon==upstream_grids[i, "lon"]), which(lat==upstream_grids[i, "lat"]), 1:dim(qg)[3]])
        setnames(qin_grid_lake, c("date","qtot","qg"))
        qin <- (qin_grid_lake$qtot+qin_grid_lake$qg)*lake$Wshd_area*1e+06 # qtot+qg [m/s] * wshd_area [m2] -> qin [m3/s]
        qin_grid_lake <- data.frame(qin_grid_lake, qin)
        qin_grids_lake <- append(qin_grids_lake, list(qin_grid_lake))
      }
      # average qin
      qin_sim <- 0
      for (i in 1:length(qin_grids_lake)) {
        qin_sim <- qin_sim+(qin_grids_lake[[i]]$qin/length(qin_grids_lake))
      }
      qin_sim <- data.frame(date, qin_sim)
      
    } else { # upstream grids > ratio 
      # selected grids
      if(method=="slope") { # for the last level to be selected, select the grid with the steepest slope.
        level <- 0
        selected_grids <- data.frame()
        while (ratio_area>0) {
          level <- level+1
          ngrids_level <- nrow(upstream_grids[upstream_grids$level==level, ])
          if (ngrids_level<=ratio_area) { # grids per level <= ratio: select all
            selected_grids <- rbind(selected_grids, upstream_grids[upstream_grids$level==level, ])
            ratio_area <- ratio_area - ngrids_level
          } else { # grids per level > ratio: select the steepest slope
            slope_grids_level <- data.frame()
            grids_level <- upstream_grids[upstream_grids$level==level, ]
            for (j in 1:nrow(grids_level)) {
              slope_grids_level0 <- data.frame(grids_level[j, ], slope[which(lon_fdir==grids_level[j, "lon"]), which(lat_fdir==grids_level[j, "lat"])]) 
              slope_grids_level <- rbind(slope_grids_level, slope_grids_level0)
            }
            colnames(slope_grids_level) <- c("lon", "lat", "level", "slope")
            slope_grids_level <- slope_grids_level[order(slope_grids_level$slope, decreasing=TRUE), ] 
            selected_grids <- rbind(selected_grids, slope_grids_level[1:ratio_area, c("lon","lat","level")])
            ratio_area <- ratio_area - nrow(slope_grids_level[1:ratio_area,])
          }
        }
      } else if (method=="all") { # for the last level to be selected, take all grids
        level <- 0
        selected_grids <- data.frame()
        while (ratio_area>0) {
          level <- level+1
          ngrids_level <- nrow(upstream_grids[upstream_grids$level==level, ])
          if (ngrids_level<=ratio_area) { # grids per level <= ratio: select all
            selected_grids <- rbind(selected_grids, upstream_grids[upstream_grids$level==level, ])
            ratio_area <- ratio_area - ngrids_level
          } else { # grids per level > ratio: select all grids at the same level
            grids_level <- upstream_grids[upstream_grids$level==level, ]
            selected_grids <- rbind(selected_grids, grids_level)
            ratio_area <- 0
          }
        }
      }
      # qtot, qg, qin for each grid 
      qin_grids_lake <- list()
      for (i in 1:nrow(selected_grids)) {
        qin_grid_lake <- data.frame(date, qtot[which(lon==selected_grids[i, "lon"]), which(lat==selected_grids[i, "lat"]), 1:dim(qtot)[3]], qg[which(lon==selected_grids[i, "lon"]), which(lat==selected_grids[i, "lat"]), 1:dim(qg)[3]])
        setnames(qin_grid_lake, c("date","qtot","qg"))
        qin <- (qin_grid_lake$qtot+qin_grid_lake$qg)*lake$Wshd_area*1e+06 # qtot+qg [m/s] * wshd_area [m2] -> qin [m3/s]
        qin_grid_lake <- data.frame(qin_grid_lake, qin)
        qin_grids_lake <- append(qin_grids_lake, list(qin_grid_lake))
      }
      # average qin
      qin_sim <- 0
      for (i in 1:length(qin_grids_lake)) {
        qin_sim <- qin_sim+(qin_grids_lake[[i]]$qin/length(qin_grids_lake))
      }
      qin_sim <- data.frame(date, qin_sim)
    }
  }
  qin_sim_yearly <- yearly_average(qin_sim) # yearly average
  colnames(qin_sim_yearly )[[2]] <- "qin_sim"
  
  # save outputs: .xlsx file 
  write_xlsx(qin_sim, file.path("outputs_hype", method, l, paste0(l,"_qin_sim_monthy",".xlsx")), col_names=TRUE)
  write_xlsx(qin_sim_yearly, file.path("outputs_hype", method, l, paste0(l,"_qin_sim_yearly",".xlsx")), col_names=TRUE)
  
  # validation
  # hype qin: https://hypeweb.smhi.se/explore-water/historical-data/europe-time-series/
  hydrolakes_hype_lakes_id <- read_excel(file.path("./inputs","hydrolakes_hype_lakes_id.xlsx"))
  hype_id <- as.numeric(hydrolakes_hype_lakes_id[hydrolakes_hype_lakes_id$hydrolakes==l, "hype"])
  qin_hype <- read_excel(file.path("./inputs/smhi_hype",paste0(hype_id,".xls")), sheet="River discharge")
  qin_hype_monthly <- monthly_average(qin_hype, lag=lag)
  colnames(qin_hype_monthly)[[2]] <- "qin_hype"
  qin_hype_yearly <- yearly_average(qin_hype)
  colnames(qin_hype_yearly)[[2]] <- "qin_hype"

  # merge: simulations and observations (grdc data, hype outputs)
  qin_sim_hype_monthly <- merge(qin_sim, qin_hype_monthly, by="date")
  qin_sim_hype_yearly <- merge(qin_sim_yearly, qin_hype_yearly, by="year")

  # plots
  # monthly
  png(file.path("outputs_hype", method, l, paste0(paste(l, "qin_sim_hype_monthly", sep="_"), ".png")), width=15, height=5, units="in", res=300)
  par(mar=c(4.1, 5, 4.1, 5)) # bottom, left, top, right
  layout(matrix(c(1,1,2), nrow=1, ncol=3, byrow=TRUE))
  # plot: simulations/observations over time
  plot(qin_sim_hype_monthly$date, qin_sim_hype_monthly$qin_sim, type="l", col="#003366", lwd=2,
       xlab="", ylab=expression("Q [m"^3*" s"^-1*"]"),
       ylim=c(0, ceiling(max(qin_sim_hype_monthly$qin_sim, qin_sim_hype_monthly$qin_hype))),
       cex.lab=1.5, cex.axis=1.5, cex=1.5)
  lines(qin_sim_hype_monthly$date, qin_sim_hype_monthly$qin_hype, col="black", lwd=2)
  legend("topleft", legend=c("Simulations", "HYPE"), col=c("#003366", "black"), lwd=2, cex=1.5, bty="n")
  # plot: simulations vs observations
  plot(qin_sim_hype_monthly$qin_sim, qin_sim_hype_monthly$qin_hype, type="p", col="#003366", pch=19,
       xlab=expression("Simulated Q [m"^3*" s"^-1*"]"),
       ylab=expression("HYPE Q [m"^3*" s"^-1*"]"),
       cex.lab=1.5, cex.axis=1.5, cex=1.5)
  abline(0, 1, col="black")
  # model performance
  mp_hype_monthly <- model_performance(qin_sim_hype_monthly$qin_hype, qin_sim_hype_monthly$qin_sim)
  mtext(paste0("KGE=", round(mp_hype_monthly$kge, 2), ", ",
               "KGEr=", round(mp_hype_monthly$kger, 2), ", ",
               "KGEb=", round(mp_hype_monthly$kgeb, 2), ", ",
               "KGEg=", round(mp_hype_monthly$kgeg, 2)),
        3, line=-2, col="black", cex=1)
  dev.off()
  # yearly
  png(file.path("outputs_hype", method, l, paste0(paste(l, "qin_sim_hype_yearly", sep="_"), ".png")), width=15, height=5, units="in", res=300)
  par(mar=c(4.1, 5, 4.1, 5)) # bottom, left, top, right
  layout(matrix(c(1,1,2), nrow=1, ncol=3, byrow=TRUE))
  # plot: simulations/observations over time
  plot(qin_sim_hype_yearly$year, qin_sim_hype_yearly$qin_sim, type="l", col="#003366", lwd=2,
       xlab="", ylab=expression("Q [m"^3*" s"^-1*"]"),
       ylim=c(0, ceiling(max(qin_sim_hype_yearly$qin_sim, qin_sim_hype_yearly$qin_hype))),
       cex.lab=1.5, cex.axis=1.5, cex=1.5)
  lines(qin_sim_hype_yearly$year, qin_sim_hype_yearly$qin_hype, col="black", lwd=2)
  legend("topleft", legend=c("Simulations", "HYPE"), col=c("#003366", "black"), lwd=2, cex=1.5, bty="n")
  # plot: simulations vs observations
  plot(qin_sim_hype_yearly$qin_sim, qin_sim_hype_yearly$qin_hype, type="p", col="#003366", pch=19,
       xlab=expression("Simulated Q [m"^3*" s"^-1*"]"),
       ylab=expression("HYPE Q [m"^3*" s"^-1*"]"),
       cex.lab=1.5, cex.axis=1.5, cex=1.5)
  abline(0, 1, col="black")
  # model performance
  mp_hype_yearly <- model_performance(qin_sim_hype_yearly$qin_hype, qin_sim_hype_yearly$qin_sim); mp_hype_yearly
  mtext(paste0("KGE=", round(mp_hype_yearly$kge, 2), ", ",
               "KGEr=", round(mp_hype_yearly$kger, 2), ", ",
               "KGEb=", round(mp_hype_yearly$kgeb, 2), ", ",
               "KGEg=", round(mp_hype_yearly$kgeg, 2)),
        3, line=-2, col="black", cex=1)
  dev.off()

  # save outputs
  outputs_hype0 <- rbind(data.table(l, nrow(lake_grids), nrow(upstream_grids), lake$Lake_area, wshd_area, grid_area, ifelse(floor(wshd_area/grid_area)<1, 1, floor(wshd_area/grid_area)), mp_hype_monthly, "monthly"),
                    data.table(l, nrow(lake_grids), nrow(upstream_grids), lake$Lake_area, wshd_area, grid_area, ifelse(floor(wshd_area/grid_area)<1, 1, floor(wshd_area/grid_area)), mp_hype_yearly, "yearly"))
  colnames(outputs_hype0) <- c("lake",
                          "lake_grids","upstream_grids","lake_area","wshd_area","grid_area","ratio_area",
                          "bias","mae","rmse","nrmse","r","nse", "kge", "kger", "kgeb", "kgeg","n","frequency")
  outputs_hype <- rbind(outputs_hype, outputs_hype0)
}
# save outputs: .xlsx file
write_xlsx(outputs_hype, file.path("outputs_hype", method,"outputs_hype_approachI.xlsx"), col_names=TRUE)

# plot KGE (KGEr, KGEb, KGEg) for the selected lakes
# monthly
df_kge_hype_monthly <- outputs_hype[frequency=="monthly", c("lake", "kge", "kger", "kgeb", "kgeg")]
dt_kge_hype_monthly_long <- melt(setDT(df_kge_hype_monthly), id.vars=c("lake"), variable.name=c("var"), measure.vars=c("kge", "kger", "kgeb", "kgeg"))
dt_kge_hype_monthly_long_mean <- dt_kge_hype_monthly_long[, lapply(.SD, mean, na.rm=TRUE), by=.(var), .SDcols=c("value")]
dt_kge_hype_monthly_long_sd <- dt_kge_hype_monthly_long[, lapply(.SD, sd, na.rm=TRUE), by=.(var), .SDcols=c("value")]
dt_kge_hype_monthly_long_mean_sd <- data.table(dt_kge_hype_monthly_long_mean$var, dt_kge_hype_monthly_long_mean$value, dt_kge_hype_monthly_long_sd$value)
setnames(dt_kge_hype_monthly_long_mean_sd, c("var","mean","sd"))
ggplot() + 
  geom_violin(data=dt_kge_hype_monthly_long, aes(x=var, y=value, fill=var), trim=FALSE, color=NA) +
  geom_jitter(data=dt_kge_hype_monthly_long, aes(x=var, y=value, color=var), alpha=0.5, size=2.5, position=position_jitter(0.2)) +
  scale_shape_manual(values=c(1, 16)) +
  scale_fill_manual(values=transparent(c("#003366","#003366","#003366","#003366"), 0.8)) +
  geom_point(data=dt_kge_hype_monthly_long_mean_sd, aes(x=var, y=mean, color=var), size=4) +
  geom_errorbar(data=dt_kge_hype_monthly_long_mean_sd, aes(x=var, ymin=mean-sd, ymax=mean+sd, color=var), size=1, width=0.1) +
  scale_color_manual(values=c("#003366","#003366","#003366","#003366")) +
  scale_y_continuous(name="", breaks=seq(from=-1,to=2.25,by=0.25)) +
  scale_x_discrete(name="", breaks=c("kge","kger","kgeb","kgeg"), labels=c("KGE",expression("KGE"[r]),expression("KGE"[b]),expression("KGE"[g]))) +
  geom_hline(yintercept=1, linetype="dashed") +
  theme(text=element_text(size=12,colour="black"),
        axis.text.y=element_text(size=12,colour="black"),
        axis.text.x=element_text(size=12,colour="black"),
        legend.position="none",
        panel.background=element_rect(fill="white"),
        panel.border=element_rect(linetype=1,fill=NA))
# yearly
df_kge_hype_yearly <- outputs_hype[frequency=="yearly", c("lake", "kge", "kger", "kgeb", "kgeg")]
dt_kge_hype_yearly_long <- melt(setDT(df_kge_hype_yearly), id.vars=c("lake"), variable.name=c("var"), measure.vars=c("kge", "kger", "kgeb", "kgeg"))
dt_kge_hype_yearly_long_mean <- dt_kge_hype_yearly_long[, lapply(.SD, mean, na.rm=TRUE), by=.(var), .SDcols=c("value")]
dt_kge_hype_yearly_long_sd <- dt_kge_hype_yearly_long[, lapply(.SD, sd, na.rm=TRUE), by=.(var), .SDcols=c("value")]
dt_kge_hype_yearly_long_mean_sd <- data.table(dt_kge_hype_yearly_long_mean$var, dt_kge_hype_yearly_long_mean$value, dt_kge_hype_yearly_long_sd$value)
setnames(dt_kge_hype_yearly_long_mean_sd, c("var","mean","sd"))
ggplot() +
  geom_violin(data=dt_kge_hype_yearly_long, aes(x=var, y=value, fill=var), trim=FALSE, color=NA) +
  geom_jitter(data=dt_kge_hype_yearly_long, aes(x=var, y=value, color=var), alpha=0.5, size=2.5, position=position_jitter(0.2)) +
  scale_shape_manual(values=c(1, 16)) +
  scale_fill_manual(values=transparent(c("#003366","#003366","#003366","#003366"), 0.8)) +
  geom_point(data=dt_kge_hype_yearly_long_mean_sd, aes(x=var, y=mean, color=var), size=4) +
  geom_errorbar(data=dt_kge_hype_yearly_long_mean_sd, aes(x=var, ymin=mean-sd, ymax=mean+sd, color=var), size=1, width=0.1) +
  scale_color_manual(values=c("#003366","#003366","#003366","#003366")) +
  scale_y_continuous(name="", breaks=seq(from=-1.25,to=3.5,by=0.5)) +
  scale_x_discrete(name="", breaks=c("kge","kger","kgeb","kgeg"), labels=c("KGE",expression("KGE"[r]),expression("KGE"[b]),expression("KGE"[g]))) +
  geom_hline(yintercept=1, linetype="dashed") +
  theme(text=element_text(size=12,colour="black"),
        axis.text.y=element_text(size=12,colour="black"),
        axis.text.x=element_text(size=12,colour="black"),
        legend.position="none",
        panel.background=element_rect(fill="white"),
        panel.border=element_rect(linetype=1,fill=NA))

# clean the environment
rm(list=ls()[!ls()%in%c("lakes", "qtot", "qg", "dis", "lat", "lon", "date", "fdir", "lat_fdir", "lon_fdir", "slope",
                        "nearest_f1", "nearest_f2", "nearest_f3", "nearest_f4", "grid_coords", "generate_grids", "generate_grids_from_center_coord",
                        "areakm2lat", "fdir_grid", "upstream_grids_of_a_grid", "upstream_grids_of_a_lake", "monthly_average", "yearly_average", "model_performance")]) 

### approach I(a/b): lake area <= grid area
## validation against observations: https://www.smhi.se/data/hydrologi/ladda-ner-hydrologiska-observationer#param=waterdischargeDaily,stations=core

# Mälaren 102, Siljan 1150, Erken 12809, 149288
df_selected_lakes <- data.frame(lake=c(102, 1150, 12809, 149288),
                                smhi=c("smhi-opendata_1_20040_20240814_141851",
                                       "smhi-opendata_1_896_20240924_121048",
                                       "smhi-opendata_1_1743_20240924_121102",
                                       "smhi-opendata_1_2152_20240423_145021"))
# method for selection of the upstream grids
method <- "slope" # slope: grid with steepest slope, all: all grids at the same level
dir.create("outputs_obs")
dir.create(file.path("./outputs_obs", method))
outputs_obs <- data.frame() # outputs
for (l in df_selected_lakes$lake) {
  # create a directory for each lake
  dir.create(file.path("outputs_obs", method, l)) 
  
  # grids occupied by the lake
  # hydrolakes + centroid inside polygon (lon,lat)
  lake <- lakes[lakes$Hylak_id==l, ] 
  # grids occupied by the lake (0.5degree x 0.5degree)
  grids <- generate_grids(lake)
  
  # find the intersection between lake and grids and divide it into parts 
  int_grids_lake <- st_intersection(grids, lake)
  # centroid inside polygon (lat, long) of each part
  int_grids_lake_cent <- st_coordinates(st_point_on_surface(int_grids_lake)) 
  colnames(int_grids_lake_cent) <- c("lon", "lat")
  # nearest grids
  lake_grids <- data.frame()
  for (i in 1:nrow(int_grids_lake_cent)) {
    lake_grids <- rbind(lake_grids, t(c(nearest_f1(int_grids_lake_cent[i, "lon"])[1], nearest_f1(int_grids_lake_cent[i, "lat"])[1]))) # nearest grid (grid it belongs to; lon, lat)
  }
  
  # upstream grids of the lake
  upstream_grids <- upstream_grids_of_a_lake(fdir, lon_fdir, lat_fdir, lake_grids)
  
  # grid area where the centroid is located (km2)
  grid_area <- areakm2lat(nearest_f1(lake$Cent_Lat)[1])
  # watershed area (km2)
  wshd_area <- lake$Wshd_area
  # number of grids to pick up: ratio between watershed area and  grid area
  ratio_area <- wshd_area/grid_area
  ratio_area <- ifelse(ratio_area<1, 1, floor(ratio_area))
  
  # qin [m3/s] = qtot+qg [m/s] * wshd_area [m2]
  if (ratio_area<=nrow(lake_grids)) { # upstream grids <= grids occupied by the lake
    # qtot, qg, qin for each grid 
    qin_grids_lake <- list()
    for (i in 1:nrow(lake_grids)) {
      qin_grid_lake <- data.frame(date, qtot[which(lon==lake_grids[i, "lon"]), which(lat==lake_grids[i, "lat"]), 1:dim(qtot)[3]], qg[which(lon==lake_grids[i, "lon"]), which(lat==lake_grids[i, "lat"]), 1:dim(qg)[3]])
      setnames(qin_grid_lake, c("date","qtot","qg"))
      qin <- (qin_grid_lake$qtot+qin_grid_lake$qg)*lake$Wshd_area*1e+06 # qtot+qg [m/s] * wshd_area [m2] -> qin [m3/s]
      qin_grid_lake <- data.frame(qin_grid_lake, qin)
      qin_grids_lake <- append(qin_grids_lake, list(qin_grid_lake))
    }
    selected_grids <- lake_grids
    # average qin
    qin_sim <- 0
    for (i in 1:length(qin_grids_lake)) {
      qin_sim <- qin_sim+(qin_grids_lake[[i]]$qin/length(qin_grids_lake))
    }
    qin_sim <- data.frame(date, qin_sim)
    
  } else if (ratio_area>nrow(lake_grids)) { # upstream grids > grids occupied by the lake
    
    if(nrow(upstream_grids)<ratio_area) { # upstream grids < ratio (select all)
      # qtot, qg, qin for each grid 
      qin_grids_lake <- list()
      for (i in 1:nrow(upstream_grids)) {
        qin_grid_lake <- data.frame(date, qtot[which(lon==upstream_grids[i, "lon"]), which(lat==upstream_grids[i, "lat"]), 1:dim(qtot)[3]], qg[which(lon==upstream_grids[i, "lon"]), which(lat==upstream_grids[i, "lat"]), 1:dim(qg)[3]])
        setnames(qin_grid_lake, c("date","qtot","qg"))
        qin <- (qin_grid_lake$qtot+qin_grid_lake$qg)*lake$Wshd_area*1e+06 # qtot+qg [m/s] * wshd_area [m2] -> qin [m3/s]
        qin_grid_lake <- data.frame(qin_grid_lake, qin)
        qin_grids_lake <- append(qin_grids_lake, list(qin_grid_lake))
      }
      # average qin
      qin_sim <- 0
      for (i in 1:length(qin_grids_lake)) {
        qin_sim <- qin_sim+(qin_grids_lake[[i]]$qin/length(qin_grids_lake))
      }
      qin_sim <- data.frame(date, qin_sim)
      
    } else { # upstream grids > ratio 
      # selected grids
      if(method=="slope") { # for the last level to be selected, select the grid with the steepest slope.
        level <- 0
        selected_grids <- data.frame()
        while (ratio_area>0) {
          level <- level+1
          ngrids_level <- nrow(upstream_grids[upstream_grids$level==level, ])
          if (ngrids_level<=ratio_area) { # grids per level <= ratio: select all
            selected_grids <- rbind(selected_grids, upstream_grids[upstream_grids$level==level, ])
            ratio_area <- ratio_area - ngrids_level
          } else { # grids per level > ratio: select the steepest slope
            slope_grids_level <- data.frame()
            grids_level <- upstream_grids[upstream_grids$level==level, ]
            for (j in 1:nrow(grids_level)) {
              slope_grids_level0 <- data.frame(grids_level[j, ], slope[which(lon_fdir==grids_level[j, "lon"]), which(lat_fdir==grids_level[j, "lat"])]) 
              slope_grids_level <- rbind(slope_grids_level, slope_grids_level0)
            }
            colnames(slope_grids_level) <- c("lon", "lat", "level", "slope")
            slope_grids_level <- slope_grids_level[order(slope_grids_level$slope, decreasing=TRUE), ] 
            selected_grids <- rbind(selected_grids, slope_grids_level[1:ratio_area, c("lon","lat","level")])
            ratio_area <- ratio_area - nrow(slope_grids_level[1:ratio_area,])
          }
        }
      } else if (method=="all") { # for the last level to be selected, take all grids
        level <- 0
        selected_grids <- data.frame()
        while (ratio_area>0) {
          level <- level+1
          ngrids_level <- nrow(upstream_grids[upstream_grids$level==level, ])
          if (ngrids_level<=ratio_area) { # grids per level <= ratio: select all
            selected_grids <- rbind(selected_grids, upstream_grids[upstream_grids$level==level, ])
            ratio_area <- ratio_area - ngrids_level
          } else { # grids per level > ratio: select all grids at the same level
            grids_level <- upstream_grids[upstream_grids$level==level, ]
            selected_grids <- rbind(selected_grids, grids_level)
            ratio_area <- 0
          }
        }
      }
      # qtot, qg, qin for each grid 
      qin_grids_lake <- list()
      for (i in 1:nrow(selected_grids)) {
        qin_grid_lake <- data.frame(date, qtot[which(lon==selected_grids[i, "lon"]), which(lat==selected_grids[i, "lat"]), 1:dim(qtot)[3]], qg[which(lon==selected_grids[i, "lon"]), which(lat==selected_grids[i, "lat"]), 1:dim(qg)[3]])
        setnames(qin_grid_lake, c("date","qtot","qg"))
        qin <- (qin_grid_lake$qtot+qin_grid_lake$qg)*lake$Wshd_area*1e+06 # qtot+qg [m/s] * wshd_area [m2] -> qin [m3/s]
        qin_grid_lake <- data.frame(qin_grid_lake, qin)
        qin_grids_lake <- append(qin_grids_lake, list(qin_grid_lake))
      }
      # average qin
      qin_sim <- 0
      for (i in 1:length(qin_grids_lake)) {
        qin_sim <- qin_sim+(qin_grids_lake[[i]]$qin/length(qin_grids_lake))
      }
      qin_sim <- data.frame(date, qin_sim)
    }
  }
  qin_sim_yearly <- yearly_average(qin_sim) # yearly average
  colnames(qin_sim_yearly )[[2]] <- "qin_sim"
  
  # save outputs: .xlsx file 
  write_xlsx(qin_sim, file.path("outputs_obs", method, l, paste0(l,"_qin_sim_monthy",".xlsx")), col_names=TRUE)
  write_xlsx(qin_sim_yearly, file.path("outputs_obs", method, l, paste0(l,"_qin_sim_yearly",".xlsx")), col_names=TRUE)
  
  # validation
  # SMHI discharge observations
  # https://www.smhi.se/data/hydrologi/ladda-ner-hydrologiska-observationer#param=waterdischargeDaily,stations=core
  # Mälaren - HydroLAKES: 102 - SMHI: 20040
  # Siljan - HydroLAKES: 1150 - SMHI: 896
  # Erken - HydroLAKES: 12809 - SMHI: 1743
  # HydroLAKES: 149288 - SMHI: 2152
  qin_obs <- fread(file.path("./inputs/smhi_obs", paste0(df_selected_lakes[df_selected_lakes$lake==l, ]$smhi,".csv")), skip=7, sep=";", fill=TRUE)
  qin_obs[, c("V3","V4","V5") := NULL]
  setnames(qin_obs, c("date", "qin"))
  qin_obs_monthly <- monthly_average(qin_obs, lag=FALSE) # lag: TRUE or FALSE
  colnames(qin_obs_monthly)[[2]] <- "qin_obs"
  qin_obs_yearly <- yearly_average(qin_obs)
  colnames(qin_obs_yearly)[[2]] <- "qin_obs"
  
  # merge: simulations and observations (grdc data, hype outputs)
  qin_sim_obs_monthly <- merge(qin_sim, qin_obs_monthly, by="date")
  qin_sim_obs_yearly <- merge(qin_sim_yearly, qin_obs_yearly, by="year")
  
  # plots
  # monthly
  png(file.path("outputs_obs", method, l, paste0(paste(l, "qin_sim_obs_monthly", sep="_"), ".png")), width=15, height=5, units="in", res=300)
  par(mar=c(4.1, 5, 4.1, 5)) # bottom, left, top, right
  layout(matrix(c(1,1,2), nrow=1, ncol=3, byrow=TRUE))
  # plot: simulations/observations over time
  plot(qin_sim_obs_monthly$date, qin_sim_obs_monthly$qin_sim, type="l", col="#990000", lwd=2,
       xlab="", ylab=expression("Q [m"^3*" s"^-1*"]"),
       ylim=c(0, ceiling(max(qin_sim_obs_monthly$qin_sim, qin_sim_obs_monthly$qin_obs))),
       cex.lab=1.5, cex.axis=1.5, cex=1.5)
  lines(qin_sim_obs_monthly$date, qin_sim_obs_monthly$qin_obs, col="black", lwd=2)
  legend("topleft", legend=c("Simulations", "Observations"), col=c("#990000", "black"), lwd=2, cex=1.5, bty="n")
  # plot: simulations vs observations
  plot(qin_sim_obs_monthly$qin_sim, qin_sim_obs_monthly$qin_obs, type="p", col="#990000", pch=19,
       xlab=expression("Simulated Q [m"^3*" s"^-1*"]"),
       ylab=expression("Observed Q [m"^3*" s"^-1*"]"),
       cex.lab=1.5, cex.axis=1.5, cex=1.5)
  abline(0, 1, col="black")
  # model performance
  mp_obs_monthly <- model_performance(qin_sim_obs_monthly$qin_obs, qin_sim_obs_monthly$qin_sim)
  mtext(paste0("KGE=", round(mp_obs_monthly$kge, 2), ", ",
               "KGEr=", round(mp_obs_monthly$kger, 2), ", ",
               "KGEb=", round(mp_obs_monthly$kgeb, 2), ", ",
               "KGEg=", round(mp_obs_monthly$kgeg, 2)),
        3, line=-2, col="black", cex=1)
  dev.off()
  # yearly
  png(file.path("outputs_obs", method, l, paste0(paste(l, "qin_sim_obs_yearly", sep="_"), ".png")), width=15, height=5, units="in", res=300)
  par(mar=c(4.1, 5, 4.1, 5)) # bottom, left, top, right
  layout(matrix(c(1,1,2), nrow=1, ncol=3, byrow=TRUE))
  # plot: simulations/observations over time
  plot(qin_sim_obs_yearly$year, qin_sim_obs_yearly$qin_sim, type="l", col="#990000", lwd=2,
       xlab="", ylab=expression("Q [m"^3*" s"^-1*"]"),
       ylim=c(0, ceiling(max(qin_sim_obs_yearly$qin_sim, qin_sim_obs_yearly$qin_obs))),
       cex.lab=1.5, cex.axis=1.5, cex=1.5)
  lines(qin_sim_obs_yearly$year, qin_sim_obs_yearly$qin_obs, col="black", lwd=2)
  legend("topleft", legend=c("Simulations", "Observations"), col=c("#990000", "black"), lwd=2, cex=1.5, bty="n")
  # plot: simulations vs observations
  plot(qin_sim_obs_yearly$qin_sim, qin_sim_obs_yearly$qin_obs, type="p", col="#990000", pch=19,
       xlab=expression("Simulated Q [m"^3*" s"^-1*"]"),
       ylab=expression("Observed Q [m"^3*" s"^-1*"]"),
       cex.lab=1.5, cex.axis=1.5, cex=1.5)
  abline(0, 1, col="black")
  # model performance
  mp_obs_yearly <- model_performance(qin_sim_obs_yearly$qin_obs, qin_sim_obs_yearly$qin_sim); mp_obs_yearly
  mtext(paste0("KGE=", round(mp_obs_yearly$kge, 2), ", ",
               "KGEr=", round(mp_obs_yearly$kger, 2), ", ",
               "KGEb=", round(mp_obs_yearly$kgeb, 2), ", ",
               "KGEg=", round(mp_obs_yearly$kgeg, 2)),
        3, line=-2, col="black", cex=1)
  dev.off()
  
  # save outputs
  outputs_obs0 <- rbind(data.table(l, nrow(lake_grids), nrow(upstream_grids), lake$Lake_area, wshd_area, grid_area, ifelse(floor(wshd_area/grid_area)<1, 1, floor(wshd_area/grid_area)), mp_obs_monthly, "monthly"),
                        data.table(l, nrow(lake_grids), nrow(upstream_grids), lake$Lake_area, wshd_area, grid_area, ifelse(floor(wshd_area/grid_area)<1, 1, floor(wshd_area/grid_area)), mp_obs_yearly, "yearly"))
  colnames(outputs_obs0) <- c("lake",
                              "lake_grids","upstream_grids","lake_area","wshd_area","grid_area","ratio_area",
                              "bias","mae","rmse","nrmse","r","nse", "kge", "kger", "kgeb", "kgeg","n","frequency")
  outputs_obs <- rbind(outputs_obs, outputs_obs0)  
}
# save outputs: .xlsx file
write_xlsx(outputs_obs, file.path("outputs_obs", method, "outputs_obs.xlsx"), col_names=TRUE)

# clean the environment
rm(list=ls()[!ls()%in%c("lakes", "qtot", "qg", "dis", "lat", "lon", "date", "fdir", "lat_fdir", "lon_fdir", "slope",
                        "nearest_f1", "nearest_f2", "nearest_f3", "nearest_f4", "grid_coords", "generate_grids", "generate_grids_from_center_coord",
                        "areakm2lat", "fdir_grid", "upstream_grids_of_a_grid", "upstream_grids_of_a_lake", "monthly_average", "yearly_average", "model_performance")]) 

### approach II: lake area > grid area
dir.create("outputs_big_lakes")

## 105 Vänern
lake <- lakes[lakes$Hylak_id==105, ]; 

# create a directory
dir.create(file.path("outputs_big_lakes", lake$Hylak_id))

# grids occupied by the lake (0.5degree x 0.5degree)
grids <- generate_grids(lake)

# find the intersection between lake and grids and divide it into parts 
int_grids_lake <- st_intersection(grids, lake)

# centroid inside polygon (lat, long) of each part
int_grids_lake_cent <- st_coordinates(st_point_on_surface(int_grids_lake)) 
colnames(int_grids_lake_cent) <- c("lon", "lat"); int_grids_lake_cent

# nearest grids
lake_grids <- data.frame()
for (i in 1:nrow(int_grids_lake_cent)) {
  lake_grids <- rbind(lake_grids, t(c(nearest_f1(int_grids_lake_cent[i, "lon"])[1], nearest_f1(int_grids_lake_cent[i, "lat"])[1]))) # nearest grid (grid it belongs to; lon, lat)
}

# upstream grids of the lake
upstream_grids <- upstream_grids_of_a_lake(fdir, lon_fdir, lat_fdir, lake_grids)

# grid area where the centroid is located (km2)
grid_area <- areakm2lat(nearest_f1(lake$Cent_Lat)[1])

# watershed area (km2)
wshd_area <- lake$Wshd_area

# number of grids to pick up: ratio between watershed area and  grid area
ratio_area <- wshd_area/grid_area
ratio_area <- ifelse(ratio_area<1, 1, floor(ratio_area))

# lake area of each grid [m2]
int_grids_lake_partial_lake_area <- as.numeric(st_area(int_grids_lake))

# area of each grid [m2]
grids_area <- unlist(lapply(lake_grids$lat, FUN=areakm2lat))*1e+06

# land area of each grid (grid area - lake area)  [m2]
int_grids_lake_partial_land_area <- grids_area-int_grids_lake_partial_lake_area 
int_grids_lake_partial_area <- data.frame(lake_grids, int_grids_lake_partial_lake_area, int_grids_lake_partial_land_area)
colnames(int_grids_lake_partial_area) <- c("lon","lat","water","land")

# selected grids
# boundary grids: # dis [m3/s]+ ((qtot+qg) [m/s] * land_area [m2]) -> qin [m3/s]
boundary_grids <- rbind(c(12.75, 59.75), c(13.25, 59.75), c(13.75, 59.75), c(14.25, 59.75), c(14.75, 59.75),
                        c(12.25, 59.25),
                        c(13.75, 58.25))
colnames(boundary_grids) <- c("lon","lat")

qin_boundary_grids_lake <- list()
for (i in 1:nrow(boundary_grids)) {
  qin_boundary_grid_lake <- data.frame(date, 
                                       dis[which(lon==boundary_grids[i, "lon"]), which(lat==boundary_grids[i, "lat"]), 1:dim(dis)[3]], 
                                       qtot[which(lon==boundary_grids[i, "lon"]), which(lat==boundary_grids[i, "lat"]), 1:dim(qtot)[3]],
                                       qg[which(lon==boundary_grids[i, "lon"]), which(lat==boundary_grids[i, "lat"]), 1:dim(qg)[3]])
  setnames(qin_boundary_grid_lake, c("date","dis","qtot","qg"))
  qin <- qin_boundary_grid_lake$dis
  qin_boundary_grid_lake <- data.frame(qin_boundary_grid_lake, qin)
  qin_boundary_grids_lake <- append(qin_boundary_grids_lake, list(qin_boundary_grid_lake))
  rm(qin_boundary_grid_lake, qin)
}

# inner grids: # (qtot+qg) [m/s] * land_area [m2] -> qin [m3/s]
inner_grids <- rbind(c(12.75, 59.25), c(13.25, 59.25), c(13.75, 59.25), c(14.25, 59.25),
                     c(12.25, 58.75), c(12.75, 58.75), c(13.25, 58.75), c(13.75, 58.75), c(14.25, 58.75),
                     c(12.75, 58.25), c(13.25, 58.25))
colnames(inner_grids) <- c("lon","lat")

qin_inner_grids_lake <- list()
for (i in 1:nrow(inner_grids)) {
  qin_inner_grid_lake <- data.frame(date, 
                                    qtot[which(lon==inner_grids[i, "lon"]), which(lat==inner_grids[i, "lat"]), 1:dim(qtot)[3]],
                                    qg[which(lon==inner_grids[i, "lon"]), which(lat==inner_grids[i, "lat"]), 1:dim(qg)[3]])
  setnames(qin_inner_grid_lake, c("date","qtot","qg"))
  qin <- (qin_inner_grid_lake$qtot+qin_inner_grid_lake$qg)*
    int_grids_lake_partial_area[which(int_grids_lake_partial_area$lon==inner_grids[i, "lon"] & int_grids_lake_partial_area$lat==inner_grids[i, "lat"]), "land"]
  qin_inner_grid_lake <- data.frame(qin_inner_grid_lake, qin)
  qin_inner_grids_lake <- append(qin_inner_grids_lake, list(qin_inner_grid_lake))
}

# total qin (sumo of qin)
qin_boundary_sim <- 0
for (i in 1:length(qin_boundary_grids_lake)) {
  qin_boundary_sim <- qin_boundary_sim+(qin_boundary_grids_lake[[i]]$qin)
}
qin_boundary_sim <- data.frame(date, qin_boundary_sim)
qin_inner_sim <- 0
for (i in 1:length(qin_inner_grids_lake)) {
  qin_inner_sim <- qin_inner_sim+(qin_inner_grids_lake[[i]]$qin)
}
qin_inner_sim <- data.frame(date, qin_inner_sim)
qin_sim <- data.frame(date, qin_sim=qin_boundary_sim$qin+qin_inner_sim$qin)
qin_sim_yearly <- yearly_average(qin_sim) # yearly average
colnames(qin_sim_yearly )[[2]] <- "qin_sim"

# save outputs: .xlsx file 
write_xlsx(qin_sim, file.path("outputs_big_lakes", lake$Hylak_id, paste0(lake$Hylak_id,"_qin_sim_monthy",".xlsx")), col_names=TRUE)
write_xlsx(qin_sim_yearly, file.path("outputs_big_lakes", lake$Hylak_id, paste0(lake$Hylak_id,"_qin_sim_yearly",".xlsx")), col_names=TRUE)

# validation
# hype qin: https://hypeweb.smhi.se/explore-water/historical-data/europe-time-series
hype_id <- 8300415
qin_hype <- read_excel(file.path("./inputs/smhi_hype",paste0(hype_id,".xls")), sheet="River discharge")
colnames(qin_hype) <- c("date","qin")
qin_hype_monthly <- monthly_average(qin_hype, lag=TRUE)
colnames(qin_hype_monthly)[[2]] <- "qin_hype"
qin_hype_yearly <- yearly_average(qin_hype)
colnames(qin_hype_yearly)[[2]] <- "qin_hype"

# observations
# SMHI discharge observations: https://www.smhi.se/data/hydrologi/ladda-ner-hydrologiska-observationer#param=waterdischargeDaily,stations=core
qin_obs <- fread(file.path("./inputs/smhi_obs","smhi-opendata_1_1954_20240814_124606.csv"), skip=7, sep=";", fill=TRUE)
qin_obs[, c("V3","V4","V5") := NULL]
setnames(qin_obs, c("date", "qin"))
qin_obs_monthly <- monthly_average(qin_obs, lag=FALSE)
colnames(qin_obs_monthly)[[2]] <- "qin_obs"
qin_obs_yearly <- yearly_average(qin_obs)
colnames(qin_obs_yearly)[[2]] <- "qin_obs"

# merge: simulations and observations
qin_sim_hype_monthly <- merge(qin_sim, qin_hype_monthly, by="date")
qin_sim_hype_yearly <- merge(qin_sim_yearly, qin_hype_yearly, by="year")
qin_sim_obs_monthly <- merge(qin_sim, qin_obs_monthly, by="date")
qin_sim_obs_yearly <- merge(qin_sim_yearly, qin_obs_yearly, by="year")

# plots
# monthly average
# hype
# plot: simulations/observations over time
png(file.path("outputs_big_lakes", lake$Hylak_id, paste0(paste(lake$Hylak_id, "qin_sim_hype_monthly", sep="_"), ".png")), width=15, height=5, units="in", res=300)
par(mar=c(4.1, 5, 4.1, 5)) # bottom, left, top, right
layout(matrix(c(1,1,2), nrow=1, ncol=3, byrow=TRUE))
plot(qin_sim_hype_monthly$date, qin_sim_hype_monthly$qin_sim, type="l", col="black", lwd=2,  
     xlab="", ylab=expression("Q [m"^3*" s"^-1*"]"), 
     ylim=c(0, ceiling(max(qin_sim_hype_monthly$qin_sim, qin_sim_hype_monthly$qin_hype))), 
     cex.lab=1.5, cex.axis=1.5, cex=1.5)
lines(qin_sim_hype_monthly$date, qin_sim_hype_monthly$qin_hype, col="#003366", lwd=2)
legend("topleft", legend=c("Simulations", "HYPE"), col=c("black","#003366"), lwd=2, cex=1.4, bty="n")
# plot: simulations vs observations 
plot(qin_sim_hype_monthly$qin_sim, qin_sim_hype_monthly$qin_hype, type="p", col="#003366", pch=19,
     xlab=expression("Simulated Q [m"^3*" s"^-1*"]"),
     ylab=expression("HYPE Q [m"^3*" s"^-1*"]"),
     cex.lab=1.5, cex.axis=1.5, cex=1.5)
abline(0, 1, col="black")
mp_hype_monthly <- model_performance(qin_sim_hype_monthly$qin_hype, qin_sim_hype_monthly$qin_sim); mp_hype_monthly
mtext(paste0("KGE=", round(mp_hype_monthly$kge, 2), ", ", 
             "KGEr=", round(mp_hype_monthly$kger, 2), ", ", 
             "KGEb=", round(mp_hype_monthly$kgeb, 2), ", ", 
             "KGEg=", round(mp_hype_monthly$kgeg, 2)), 
      3, line=-2, col="black", cex=1)
dev.off()
# obs
# plot: simulations/observations over time
png(file.path("outputs_big_lakes", lake$Hylak_id, paste0(paste(lake$Hylak_id, "qin_sim_obs_monthly", sep="_"), ".png")), width=15, height=5, units="in", res=300)
par(mar=c(4.1, 5, 4.1, 5)) # bottom, left, top, right
layout(matrix(c(1,1,2), nrow=1, ncol=3, byrow=TRUE))
plot(qin_sim_obs_monthly$date, qin_sim_obs_monthly$qin_sim, type="l", col="black", lwd=2,  
     xlab="", ylab=expression("Q [m"^3*" s"^-1*"]"), 
     ylim=c(0, ceiling(max(qin_sim_obs_monthly$qin_sim, qin_sim_obs_monthly$qin_obs))), 
     cex.lab=1.5, cex.axis=1.5, cex=1.5)
lines(qin_sim_obs_monthly$date, qin_sim_obs_monthly$qin_obs, col="#990000", lwd=2)
legend("topleft", legend=c("Simulations", "Observations"), col=c("black","#990000"), lwd=2, cex=1.4, bty="n")
# plot: simulations vs observations 
plot(qin_sim_obs_monthly$qin_sim, qin_sim_obs_monthly$qin_obs, type="p", col="#990000", pch=19,
     xlab=expression("Simulated Q [m"^3*" s"^-1*"]"),
     ylab=expression("Observed Q [m"^3*" s"^-1*"]"),
     cex.lab=1.5, cex.axis=1.5, cex=1.5)
abline(0, 1, col="black")
mp_obs_monthly <- model_performance(qin_sim_obs_monthly$qin_obs, qin_sim_obs_monthly$qin_sim); mp_obs_monthly
mtext(paste0("KGE=", round(mp_obs_monthly$kge, 2), ", ", 
             "KGEr=", round(mp_obs_monthly$kger, 2), ", ", 
             "KGEb=", round(mp_obs_monthly$kgeb, 2), ", ", 
             "KGEg=", round(mp_obs_monthly$kgeg, 2)), 
      3, line=-2, col="black", cex=1)
dev.off()
# yearly average
# hype
# plot: simulations/observations over time
png(file.path("outputs_big_lakes", lake$Hylak_id, paste0(paste(lake$Hylak_id, "qin_sim_hype_yearly", sep="_"), ".png")), width=15, height=5, units="in", res=300)
par(mar=c(4.1, 5, 4.1, 5)) # bottom, left, top, right
layout(matrix(c(1,1,2), nrow=1, ncol=3, byrow=TRUE))
plot(qin_sim_hype_yearly$year, qin_sim_hype_yearly$qin_sim, type="l", col="black", lwd=2,  
     xlab="", ylab=expression("Q [m"^3*" s"^-1*"]"), 
     ylim=c(0, ceiling(max(qin_sim_hype_yearly$qin_sim, qin_sim_hype_yearly$qin_hype))), 
     cex.lab=1.5, cex.axis=1.5, cex=1.5)
lines(qin_sim_hype_yearly$year, qin_sim_hype_yearly$qin_hype, col="#003366", lwd=2)
legend("topleft", legend=c("Simulations", "HYPE"), col=c("black","#003366"), lwd=2, cex=1.4, bty="n")
# plot: simulations vs observations 
plot(qin_sim_hype_yearly$qin_sim, qin_sim_hype_yearly$qin_hype, type="p", col="#003366", pch=19,
     xlab=expression("Simulated Q [m"^3*" s"^-1*"]"),
     ylab=expression("HYPE Q [m"^3*" s"^-1*"]"),
     cex.lab=1.5, cex.axis=1.5, cex=1.5)
abline(0, 1, col="black")
mp_hype_yearly <- model_performance(qin_sim_hype_yearly$qin_hype, qin_sim_hype_yearly$qin_sim); mp_hype_yearly
mtext(paste0("KGE=", round(mp_hype_yearly$kge, 2), ", ", 
             "KGEr=", round(mp_hype_yearly$kger, 2), ", ", 
             "KGEb=", round(mp_hype_yearly$kgeb, 2), ", ", 
             "KGEg=", round(mp_hype_yearly$kgeg, 2)), 
      3, line=-2, col="black", cex=1)
dev.off()
# obs
# plot: simulations/observations over time
png(file.path("outputs_big_lakes", lake$Hylak_id, paste0(paste(lake$Hylak_id, "qin_sim_obs_yearly", sep="_"), ".png")), width=15, height=5, units="in", res=300)
par(mar=c(4.1, 5, 4.1, 5)) # bottom, left, top, right
layout(matrix(c(1,1,2), nrow=1, ncol=3, byrow=TRUE))
plot(qin_sim_obs_yearly$year, qin_sim_obs_yearly$qin_sim, type="l", col="black", lwd=2,  
     xlab="", ylab=expression("Q [m"^3*" s"^-1*"]"), 
     ylim=c(0, ceiling(max(qin_sim_obs_yearly$qin_sim, qin_sim_obs_yearly$qin_obs))), 
     cex.lab=1.5, cex.axis=1.5, cex=1.5)
lines(qin_sim_obs_yearly$year, qin_sim_obs_yearly$qin_obs, col="#990000", lwd=2)
legend("topleft", legend=c("Simulations", "Observations"), col=c("black","#990000"), lwd=2, cex=1.4, bty="n")
# plot: simulations vs observations 
plot(qin_sim_obs_yearly$qin_sim, qin_sim_obs_yearly$qin_obs, type="p", col="#990000", pch=19,
     xlab=expression("Simulated Q [m"^3*" s"^-1*"]"),
     ylab=expression("Observed Q [m"^3*" s"^-1*"]"),
     cex.lab=1.5, cex.axis=1.5, cex=1.5)
abline(0, 1, col="black")
mp_obs_yearly <- model_performance(qin_sim_obs_yearly$qin_obs, qin_sim_obs_yearly$qin_sim); mp_obs_yearly
mtext(paste0("KGE=", round(mp_obs_yearly$kge, 2), ", ", 
             "KGEr=", round(mp_obs_yearly$kger, 2), ", ", 
             "KGEb=", round(mp_obs_yearly$kgeb, 2), ", ", 
             "KGEg=", round(mp_obs_yearly$kgeg, 2)), 
      3, line=-2, col="black", cex=1)
dev.off()

# clean the environment
rm(list=ls()[!ls()%in%c("lakes", "qtot", "qg", "dis", "lat", "lon", "date", "fdir", "lat_fdir", "lon_fdir", "slope",
                        "nearest_f1", "nearest_f2", "nearest_f3", "nearest_f4", "grid_coords", "generate_grids", "generate_grids_from_center_coord",
                        "areakm2lat", "fdir_grid", "upstream_grids_of_a_grid", "upstream_grids_of_a_lake", "monthly_average", "yearly_average", "model_performance")]) 

## 104 Vättern
lake <- lakes[lakes$Hylak_id==104, ]

# create a directory
dir.create(file.path("outputs_big_lakes", lake$Hylak_id))

# grids occupied by the lake (0.5degree x 0.5degree)
grids <- generate_grids(lake)

# find the intersection between lake and grids and divide it into parts 
int_grids_lake <- st_intersection(grids, lake)

# centroid inside polygon (lat, long) of each part
int_grids_lake_cent <- st_coordinates(st_point_on_surface(int_grids_lake)) 
colnames(int_grids_lake_cent) <- c("lon", "lat")

# nearest grids
lake_grids <- data.frame()
for (i in 1:nrow(int_grids_lake_cent)) {
  lake_grids <- rbind(lake_grids, t(c(nearest_f1(int_grids_lake_cent[i, "lon"])[1], nearest_f1(int_grids_lake_cent[i, "lat"])[1]))) # nearest grid (grid it belongs to; lon, lat)
}

# upstream grids of the lake
upstream_grids <- upstream_grids_of_a_lake(fdir, lon_fdir, lat_fdir, lake_grids)

# grid area where the centroid is located (km2)
grid_area <- areakm2lat(nearest_f1(lake$Cent_Lat)[1])

# watershed area (km2)
wshd_area <- lake$Wshd_area

# number of grids to pick up: ratio between watershed area and  grid area
ratio_area <- wshd_area/grid_area
ratio_area <- ifelse(ratio_area<1, 1, floor(ratio_area))

# lake area of each grid [m2]
int_grids_lake_partial_lake_area <- as.numeric(st_area(int_grids_lake))

# area of each grid [m2]
grids_area <- unlist(lapply(lake_grids$lat, FUN=areakm2lat))*1e+06

# land area of each grid (grid area - lake area)  [m2]
int_grids_lake_partial_land_area <- grids_area-int_grids_lake_partial_lake_area 
int_grids_lake_partial_area <- data.frame(lake_grids, int_grids_lake_partial_lake_area, int_grids_lake_partial_land_area)
colnames(int_grids_lake_partial_area) <- c("lon","lat","water","land")

# selected grids
# boundary grids: # dis [m3/s] -> qin [m3/s]
boundary_grids <- data.frame(14.75, 57.75)
colnames(boundary_grids) <- c("lon","lat")

qin_boundary_grids_lake <- list()
for (i in 1:nrow(boundary_grids)) {
  qin_boundary_grid_lake <- data.frame(date, 
                                       dis[which(lon==boundary_grids[i, "lon"]), which(lat==boundary_grids[i, "lat"]), 1:dim(dis)[3]])
  setnames(qin_boundary_grid_lake, c("date","dis"))
  qin <- qin_boundary_grid_lake$dis
  qin_boundary_grid_lake <- data.frame(qin_boundary_grid_lake, qin)
  qin_boundary_grids_lake <- append(qin_boundary_grids_lake, list(qin_boundary_grid_lake))
  rm(qin_boundary_grid_lake, qin)
}

# inner grids: # (qtot+qg) [m/s] * land_area [m2] -> qin [m3/s]
inner_grids <- rbind(c(14.75, 58.75),
                     c(14.25, 58.25),
                     c(14.25, 57.75),
                     c(14.75, 58.25))
colnames(inner_grids) <- c("lon","lat")

qin_inner_grids_lake <- list()
for (i in 1:nrow(inner_grids)) {
  qin_inner_grid_lake <- data.frame(date, 
                                    qtot[which(lon==inner_grids[i, "lon"]), which(lat==inner_grids[i, "lat"]), 1:dim(qtot)[3]],
                                    qg[which(lon==inner_grids[i, "lon"]), which(lat==inner_grids[i, "lat"]), 1:dim(qg)[3]])
  setnames(qin_inner_grid_lake, c("date","qtot","qg"))
  qin <- (qin_inner_grid_lake$qtot+qin_inner_grid_lake$qg)*
    int_grids_lake_partial_area[which(int_grids_lake_partial_area$lon==inner_grids[i, "lon"] & int_grids_lake_partial_area$lat==inner_grids[i, "lat"]), "land"]
  qin_inner_grid_lake <- data.frame(qin_inner_grid_lake, qin)
  qin_inner_grids_lake <- append(qin_inner_grids_lake, list(qin_inner_grid_lake))
}

# total qin (sumo of qin)
qin_boundary_sim <- 0
for (i in 1:length(qin_boundary_grids_lake)) {
  qin_boundary_sim <- qin_boundary_sim+(qin_boundary_grids_lake[[i]]$qin)
}
qin_boundary_sim <- data.frame(date, qin_boundary_sim)
qin_inner_sim <- 0
for (i in 1:length(qin_inner_grids_lake)) {
  qin_inner_sim <- qin_inner_sim+(qin_inner_grids_lake[[i]]$qin)
}
qin_inner_sim <- data.frame(date, qin_inner_sim)
qin_sim <- data.frame(date, qin_sim=qin_boundary_sim$qin+qin_inner_sim$qin)
qin_sim_yearly <- yearly_average(qin_sim) # yearly average
colnames(qin_sim_yearly )[[2]] <- "qin_sim"

# save outputs: .xlsx file 
write_xlsx(qin_sim, file.path("outputs_big_lakes", lake$Hylak_id, paste0(lake$Hylak_id,"_qin_sim_monthy",".xlsx")), col_names=TRUE)
write_xlsx(qin_sim_yearly, file.path("outputs_big_lakes", lake$Hylak_id, paste0(lake$Hylak_id,"_qin_sim_yearly",".xlsx")), col_names=TRUE)

# validation
# hype qin: https://hypeweb.smhi.se/explore-water/historical-data/europe-time-series
hype_id <- c(8322152, 8322006, 8322153, 8322101, 8302040, 8302027, 8302061, 8321638, 8321545, 8321647, 8321529, 8321587, 
             8302060, 8321534, 8302051, 8321639, 8000438, 8302044, 8321525, 8302033, 8000452, 8322100, 8302038)
qin <- 0
for (i in hype_id) {
  qin_hype0 <- read_excel(file.path("./inputs/smhi_hype",paste0(i,".xls")), sheet="River discharge")
  colnames(qin_hype0) <- c("date","qin")
  qin <- qin+(qin_hype0$qin)
}
qin_hype <- data.frame(qin_hype0$date, qin)
colnames(qin_hype) <- c("date", "qin")
qin_hype_monthly <- monthly_average(qin_hype, lag=FALSE)
colnames(qin_hype_monthly)[[2]] <- "qin_hype"
qin_hype_yearly <- yearly_average(qin_hype)
colnames(qin_hype_yearly)[[2]] <- "qin_hype"

# observations
# SMHI discharge observations: https://www.smhi.se/data/hydrologi/ladda-ner-hydrologiska-observationer#param=waterdischargeDaily,stations=core
qin_obs <- fread(file.path("./inputs/smhi_obs","smhi-opendata_1_1950_20240815_081417.csv"), skip=7, sep=";", fill=TRUE)
qin_obs[, c("V3","V4","V5") := NULL]
setnames(qin_obs, c("date", "qin"))
qin_obs_monthly <- monthly_average(qin_obs, lag=TRUE)
colnames(qin_obs_monthly)[[2]] <- "qin_obs"
qin_obs_yearly <- yearly_average(qin_obs)
colnames(qin_obs_yearly)[[2]] <- "qin_obs"

# merge: simulations and observations
qin_sim_hype_monthly <- merge(qin_sim, qin_hype_monthly, by="date")
qin_sim_hype_yearly <- merge(qin_sim_yearly, qin_hype_yearly, by="year")
qin_sim_obs_monthly <- merge(qin_sim, qin_obs_monthly, by="date")
qin_sim_obs_yearly <- merge(qin_sim_yearly, qin_obs_yearly, by="year")

# plots
# monthly average
# hype
# plot: simulations/observations over time
png(file.path("outputs_big_lakes", lake$Hylak_id, paste0(paste(lake$Hylak_id, "qin_sim_hype_monthly", sep="_"), ".png")), width=15, height=5, units="in", res=300)
par(mar=c(4.1, 5, 4.1, 5)) # bottom, left, top, right
layout(matrix(c(1,1,2), nrow=1, ncol=3, byrow=TRUE))
plot(qin_sim_hype_monthly$date, qin_sim_hype_monthly$qin_sim, type="l", col="black", lwd=2,  
     xlab="", ylab=expression("Q [m"^3*" s"^-1*"]"), 
     ylim=c(0, ceiling(max(qin_sim_hype_monthly$qin_sim, qin_sim_hype_monthly$qin_hype))), 
     cex.lab=1.5, cex.axis=1.5, cex=1.5)
lines(qin_sim_hype_monthly$date, qin_sim_hype_monthly$qin_hype, col="#003366", lwd=2)
legend("topleft", legend=c("Simulations", "HYPE"), col=c("black","#003366"), lwd=2, cex=1.4, bty="n")
# plot: simulations vs observations 
plot(qin_sim_hype_monthly$qin_sim, qin_sim_hype_monthly$qin_hype, type="p", col="#003366", pch=19,
     xlab=expression("Simulated Q [m"^3*" s"^-1*"]"),
     ylab=expression("HYPE Q [m"^3*" s"^-1*"]"),
     cex.lab=1.5, cex.axis=1.5, cex=1.5)
abline(0, 1, col="black")
mp_hype_monthly <- model_performance(qin_sim_hype_monthly$qin_hype, qin_sim_hype_monthly$qin_sim); mp_hype_monthly
mtext(paste0("KGE=", round(mp_hype_monthly$kge, 2), ", ", 
             "KGEr=", round(mp_hype_monthly$kger, 2), ", ", 
             "KGEb=", round(mp_hype_monthly$kgeb, 2), ", ", 
             "KGEg=", round(mp_hype_monthly$kgeg, 2)), 
      3, line=-2, col="black", cex=1)
dev.off()
# obs
# plot: simulations/observations over time
png(file.path("outputs_big_lakes", lake$Hylak_id, paste0(paste(lake$Hylak_id, "qin_sim_obs_monthly", sep="_"), ".png")), width=15, height=5, units="in", res=300)
par(mar=c(4.1, 5, 4.1, 5)) # bottom, left, top, right
layout(matrix(c(1,1,2), nrow=1, ncol=3, byrow=TRUE))
plot(qin_sim_obs_monthly$date, qin_sim_obs_monthly$qin_sim, type="l", col="black", lwd=2,  
     xlab="", ylab=expression("Q [m"^3*" s"^-1*"]"), 
     ylim=c(0, ceiling(max(qin_sim_obs_monthly$qin_sim, qin_sim_obs_monthly$qin_obs))), 
     cex.lab=1.5, cex.axis=1.5, cex=1.5)
lines(qin_sim_obs_monthly$date, qin_sim_obs_monthly$qin_obs, col="#990000", lwd=2)
legend("topleft", legend=c("Simulations", "Observations"), col=c("black","#990000"), lwd=2, cex=1.4, bty="n")
# plot: simulations vs observations 
plot(qin_sim_obs_monthly$qin_sim, qin_sim_obs_monthly$qin_obs, type="p", col="#990000", pch=19,
     xlab=expression("Simulated Q [m"^3*" s"^-1*"]"),
     ylab=expression("Observed Q [m"^3*" s"^-1*"]"),
     cex.lab=1.5, cex.axis=1.5, cex=1.5)
abline(0, 1, col="black")
mp_obs_monthly <- model_performance(qin_sim_obs_monthly$qin_obs, qin_sim_obs_monthly$qin_sim); mp_obs_monthly
mtext(paste0("KGE=", round(mp_obs_monthly$kge, 2), ", ", 
             "KGEr=", round(mp_obs_monthly$kger, 2), ", ", 
             "KGEb=", round(mp_obs_monthly$kgeb, 2), ", ", 
             "KGEg=", round(mp_obs_monthly$kgeg, 2)), 
      3, line=-2, col="black", cex=1)
dev.off()
# yearly average
# hype
# plot: simulations/observations over time
png(file.path("outputs_big_lakes", lake$Hylak_id, paste0(paste(lake$Hylak_id, "qin_sim_hype_yearly", sep="_"), ".png")), width=15, height=5, units="in", res=300)
par(mar=c(4.1, 5, 4.1, 5)) # bottom, left, top, right
layout(matrix(c(1,1,2), nrow=1, ncol=3, byrow=TRUE))
plot(qin_sim_hype_yearly$year, qin_sim_hype_yearly$qin_sim, type="l", col="black", lwd=2,  
     xlab="", ylab=expression("Q [m"^3*" s"^-1*"]"), 
     ylim=c(0, ceiling(max(qin_sim_hype_yearly$qin_sim, qin_sim_hype_yearly$qin_hype))), 
     cex.lab=1.5, cex.axis=1.5, cex=1.5)
lines(qin_sim_hype_yearly$year, qin_sim_hype_yearly$qin_hype, col="#003366", lwd=2)
legend("topleft", legend=c("Simulations", "HYPE"), col=c("black","#003366"), lwd=2, cex=1.4, bty="n")
# plot: simulations vs observations 
plot(qin_sim_hype_yearly$qin_sim, qin_sim_hype_yearly$qin_hype, type="p", col="#003366", pch=19,
     xlab=expression("Simulated Q [m"^3*" s"^-1*"]"),
     ylab=expression("HYPE Q [m"^3*" s"^-1*"]"),
     cex.lab=1.5, cex.axis=1.5, cex=1.5)
abline(0, 1, col="black")
mp_hype_yearly <- model_performance(qin_sim_hype_yearly$qin_hype, qin_sim_hype_yearly$qin_sim); mp_hype_yearly
mtext(paste0("KGE=", round(mp_hype_yearly$kge, 2), ", ", 
             "KGEr=", round(mp_hype_yearly$kger, 2), ", ", 
             "KGEb=", round(mp_hype_yearly$kgeb, 2), ", ", 
             "KGEg=", round(mp_hype_yearly$kgeg, 2)), 
      3, line=-2, col="black", cex=1)
dev.off()
# obs
# plot: simulations/observations over time
png(file.path("outputs_big_lakes", lake$Hylak_id, paste0(paste(lake$Hylak_id, "qin_sim_obs_yearly", sep="_"), ".png")), width=15, height=5, units="in", res=300)
par(mar=c(4.1, 5, 4.1, 5)) # bottom, left, top, right
layout(matrix(c(1,1,2), nrow=1, ncol=3, byrow=TRUE))
plot(qin_sim_obs_yearly$year, qin_sim_obs_yearly$qin_sim, type="l", col="black", lwd=2,  
     xlab="", ylab=expression("Q [m"^3*" s"^-1*"]"), 
     ylim=c(0, ceiling(max(qin_sim_obs_yearly$qin_sim, qin_sim_obs_yearly$qin_obs))), 
     cex.lab=1.5, cex.axis=1.5, cex=1.5)
lines(qin_sim_obs_yearly$year, qin_sim_obs_yearly$qin_obs, col="#990000", lwd=2)
legend("topleft", legend=c("Simulations", "Observations"), col=c("black","#990000"), lwd=2, cex=1.4, bty="n")
# plot: simulations vs observations 
plot(qin_sim_obs_yearly$qin_sim, qin_sim_obs_yearly$qin_obs, type="p", col="#990000", pch=19,
     xlab=expression("Simulated Q [m"^3*" s"^-1*"]"),
     ylab=expression("Observed Q [m"^3*" s"^-1*"]"),
     cex.lab=1.5, cex.axis=1.5, cex=1.5)
abline(0, 1, col="black")
mp_obs_yearly <- model_performance(qin_sim_obs_yearly$qin_obs, qin_sim_obs_yearly$qin_sim); mp_obs_yearly
mtext(paste0("KGE=", round(mp_obs_yearly$kge, 2), ", ", 
             "KGEr=", round(mp_obs_yearly$kger, 2), ", ", 
             "KGEb=", round(mp_obs_yearly$kgeb, 2), ", ", 
             "KGEg=", round(mp_obs_yearly$kgeg, 2)), 
      3, line=-2, col="black", cex=1)
dev.off()

# clean the environment
rm(list=ls()[!ls()%in%c("lakes", "qtot", "qg", "dis", "lat", "lon", "date", "fdir", "lat_fdir", "lon_fdir", "slope",
                        "nearest_f1", "nearest_f2", "nearest_f3", "nearest_f4", "grid_coords", "generate_grids", "generate_grids_from_center_coord",
                        "areakm2lat", "fdir_grid", "upstream_grids_of_a_grid", "upstream_grids_of_a_lake", "monthly_average", "yearly_average", "model_performance")]) 

## 102 Mälaren
lake <- lakes[lakes$Hylak_id==102, ]

# create a directory
dir.create(file.path("outputs_big_lakes", lake$Hylak_id))

# grids occupied by the lake (0.5degree x 0.5degree)
grids <- generate_grids(lake) 

# find the intersection between lake and grids and divide it into parts 
int_grids_lake <- st_intersection(grids, lake)

# centroid inside polygon (lat, long) of each part
int_grids_lake_cent <- st_coordinates(st_point_on_surface(int_grids_lake)) 
colnames(int_grids_lake_cent) <- c("lon", "lat")

# nearest grids
lake_grids <- data.frame()
for (i in 1:nrow(int_grids_lake_cent)) {
  lake_grids <- rbind(lake_grids, t(c(nearest_f1(int_grids_lake_cent[i, "lon"])[1], nearest_f1(int_grids_lake_cent[i, "lat"])[1]))) # nearest grid (grid it belongs to; lon, lat)
}

# upstream grids of the lake
upstream_grids <- upstream_grids_of_a_lake(fdir, lon_fdir, lat_fdir, lake_grids)

# grid area where the centroid is located (km2)
grid_area <- areakm2lat(nearest_f1(lake$Cent_Lat)[1])

# watershed area (km2)
wshd_area <- lake$Wshd_area

# number of grids to pick up: ratio between watershed area and  grid area
ratio_area <- wshd_area/grid_area
ratio_area <- ifelse(ratio_area<1, 1, floor(ratio_area))

# lake area of each grid [m2]
int_grids_lake_partial_lake_area <- as.numeric(st_area(int_grids_lake))

# area of each grid [m2]
grids_area <- unlist(lapply(lake_grids$lat, FUN=areakm2lat))*1e+06

# land area of each grid (grid area - lake area)  [m2]
int_grids_lake_partial_land_area <- grids_area-int_grids_lake_partial_lake_area 
int_grids_lake_partial_area <- data.frame(lake_grids, int_grids_lake_partial_lake_area, int_grids_lake_partial_land_area)
colnames(int_grids_lake_partial_area) <- c("lon","lat","water","land")

# selected grids
# boundary grids: # dis [m3/s] -> qin [m3/s]
boundary_grids <- rbind(c(15.75, 59.75), c(15.75, 59.25))
colnames(boundary_grids) <- c("lon","lat")

qin_boundary_grids_lake <- list()
for (i in 1:nrow(boundary_grids)) {
  qin_boundary_grid_lake <- data.frame(date, 
                                       dis[which(lon==boundary_grids[i, "lon"]), which(lat==boundary_grids[i, "lat"]), 1:dim(dis)[3]], 
                                       qtot[which(lon==boundary_grids[i, "lon"]), which(lat==boundary_grids[i, "lat"]), 1:dim(qtot)[3]],
                                       qg[which(lon==boundary_grids[i, "lon"]), which(lat==boundary_grids[i, "lat"]), 1:dim(qg)[3]])
  setnames(qin_boundary_grid_lake, c("date","dis","qtot","qg"))
  qin <- qin_boundary_grid_lake$dis
  qin_boundary_grid_lake <- data.frame(qin_boundary_grid_lake, qin)
  qin_boundary_grids_lake <- append(qin_boundary_grids_lake, list(qin_boundary_grid_lake))
  rm(qin_boundary_grid_lake, qin)
}

# inner grids: # (qtot+qg) [m/s] * land_area [m2] -> qin [m3/s]
inner_grids <- rbind(c(16.25, 59.75), c(16.75, 59.75), c(17.25, 59.75), c(17.75, 59.75),
                     c(16.25, 59.25), c(16.75, 59.25), c(17.25, 59.25))
colnames(inner_grids) <- c("lon","lat")

qin_inner_grids_lake <- list()
for (i in 1:nrow(inner_grids)) {
  qin_inner_grid_lake <- data.frame(date, 
                                    qtot[which(lon==inner_grids[i, "lon"]), which(lat==inner_grids[i, "lat"]), 1:dim(qtot)[3]],
                                    qg[which(lon==inner_grids[i, "lon"]), which(lat==inner_grids[i, "lat"]), 1:dim(qg)[3]])
  setnames(qin_inner_grid_lake, c("date","qtot","qg"))
  qin <- (qin_inner_grid_lake$qtot+qin_inner_grid_lake$qg)*
    int_grids_lake_partial_area[which(int_grids_lake_partial_area$lon==inner_grids[i, "lon"] & int_grids_lake_partial_area$lat==inner_grids[i, "lat"]), "land"]
  qin_inner_grid_lake <- data.frame(qin_inner_grid_lake, qin)
  qin_inner_grids_lake <- append(qin_inner_grids_lake, list(qin_inner_grid_lake))
}

# total qin (sumo of qin)
qin_boundary_sim <- 0
for (i in 1:length(qin_boundary_grids_lake)) {
  qin_boundary_sim <- qin_boundary_sim+(qin_boundary_grids_lake[[i]]$qin)
}
qin_boundary_sim <- data.frame(date, qin_boundary_sim)
qin_inner_sim <- 0
for (i in 1:length(qin_inner_grids_lake)) {
  qin_inner_sim <- qin_inner_sim+(qin_inner_grids_lake[[i]]$qin)
}
qin_inner_sim <- data.frame(date, qin_inner_sim)
qin_sim <- data.frame(date, qin_sim=qin_boundary_sim$qin+qin_inner_sim$qin)
qin_sim_yearly <- yearly_average(qin_sim) # yearly average
colnames(qin_sim_yearly )[[2]] <- "qin_sim"

# save outputs: .xlsx file 
write_xlsx(qin_sim, file.path("outputs_big_lakes", lake$Hylak_id, paste0(lake$Hylak_id,"_qin_sim_monthy",".xlsx")), col_names=TRUE)
write_xlsx(qin_sim_yearly, file.path("outputs_big_lakes", lake$Hylak_id, paste0(lake$Hylak_id,"_qin_sim_yearly",".xlsx")), col_names=TRUE)

# validation
# hype qin: https://hypeweb.smhi.se/explore-water/historical-data/europe-time-series
hype_id <- 8302852
qin_hype <- read_excel(file.path("./inputs/smhi_hype",paste0(hype_id,".xls")), sheet="River discharge")
colnames(qin_hype) <- c("date","qin")
qin_hype_monthly <- monthly_average(qin_hype, lag=TRUE)
colnames(qin_hype_monthly)[[2]] <- "qin_hype"
qin_hype_yearly <- yearly_average(qin_hype)
colnames(qin_hype_yearly)[[2]] <- "qin_hype"

# observations
# SMHI discharge observations: https://www.smhi.se/data/hydrologi/ladda-ner-hydrologiska-observationer#param=waterdischargeDaily,stations=core
qin_obs <- fread(file.path("./inputs/smhi_obs","smhi-opendata_1_20040_20240814_141851.csv"), skip=7, sep=";", fill=TRUE)
qin_obs[, c("V3","V4","V5") := NULL]
setnames(qin_obs, c("date", "qin"))
qin_obs_monthly <- monthly_average(qin_obs, lag=FALSE)
colnames(qin_obs_monthly)[[2]] <- "qin_obs"
qin_obs_yearly <- yearly_average(qin_obs)
colnames(qin_obs_yearly)[[2]] <- "qin_obs"

# merge: simulations and observations (grdc data, hype outputs)
qin_sim_hype_monthly <- merge(qin_sim, qin_hype_monthly, by="date")
qin_sim_hype_yearly <- merge(qin_sim_yearly, qin_hype_yearly, by="year")
qin_sim_obs_monthly <- merge(qin_sim, qin_obs_monthly, by="date")
qin_sim_obs_yearly <- merge(qin_sim_yearly, qin_obs_yearly, by="year")

# plots
# monthly average
# hype
# plot: simulations/observations over time
png(file.path("outputs_big_lakes", lake$Hylak_id, paste0(paste(lake$Hylak_id, "qin_sim_hype_monthly", sep="_"), ".png")), width=15, height=5, units="in", res=300)
par(mar=c(4.1, 5, 4.1, 5)) # bottom, left, top, right
layout(matrix(c(1,1,2), nrow=1, ncol=3, byrow=TRUE))
plot(qin_sim_hype_monthly$date, qin_sim_hype_monthly$qin_sim, type="l", col="black", lwd=2,  
     xlab="", ylab=expression("Q [m"^3*" s"^-1*"]"), 
     ylim=c(0, ceiling(max(qin_sim_hype_monthly$qin_sim, qin_sim_hype_monthly$qin_hype))), 
     cex.lab=1.5, cex.axis=1.5, cex=1.5)
lines(qin_sim_hype_monthly$date, qin_sim_hype_monthly$qin_hype, col="#003366", lwd=2)
legend("topleft", legend=c("Simulations", "HYPE"), col=c("black","#003366"), lwd=2, cex=1.4, bty="n")
# plot: simulations vs observations 
plot(qin_sim_hype_monthly$qin_sim, qin_sim_hype_monthly$qin_hype, type="p", col="#003366", pch=19,
     xlab=expression("Simulated Q [m"^3*" s"^-1*"]"),
     ylab=expression("HYPE Q [m"^3*" s"^-1*"]"),
     cex.lab=1.5, cex.axis=1.5, cex=1.5)
abline(0, 1, col="black")
mp_hype_monthly <- model_performance(qin_sim_hype_monthly$qin_hype, qin_sim_hype_monthly$qin_sim); mp_hype_monthly
mtext(paste0("KGE=", round(mp_hype_monthly$kge, 2), ", ", 
             "KGEr=", round(mp_hype_monthly$kger, 2), ", ", 
             "KGEb=", round(mp_hype_monthly$kgeb, 2), ", ", 
             "KGEg=", round(mp_hype_monthly$kgeg, 2)), 
      3, line=-2, col="black", cex=1)
dev.off()
# obs
# plot: simulations/observations over time
png(file.path("outputs_big_lakes", lake$Hylak_id, paste0(paste(lake$Hylak_id, "qin_sim_obs_monthly", sep="_"), ".png")), width=15, height=5, units="in", res=300)
par(mar=c(4.1, 5, 4.1, 5)) # bottom, left, top, right
layout(matrix(c(1,1,2), nrow=1, ncol=3, byrow=TRUE))
plot(qin_sim_obs_monthly$date, qin_sim_obs_monthly$qin_sim, type="l", col="black", lwd=2,  
     xlab="", ylab=expression("Q [m"^3*" s"^-1*"]"), 
     ylim=c(0, ceiling(max(qin_sim_obs_monthly$qin_sim, qin_sim_obs_monthly$qin_obs))), 
     cex.lab=1.5, cex.axis=1.5, cex=1.5)
lines(qin_sim_obs_monthly$date, qin_sim_obs_monthly$qin_obs, col="#990000", lwd=2)
legend("topleft", legend=c("Simulations", "Observations"), col=c("black","#990000"), lwd=2, cex=1.4, bty="n")
# plot: simulations vs observations 
plot(qin_sim_obs_monthly$qin_sim, qin_sim_obs_monthly$qin_obs, type="p", col="#990000", pch=19,
     xlab=expression("Simulated Q [m"^3*" s"^-1*"]"),
     ylab=expression("Observed Q [m"^3*" s"^-1*"]"),
     cex.lab=1.5, cex.axis=1.5, cex=1.5)
abline(0, 1, col="black")
mp_obs_monthly <- model_performance(qin_sim_obs_monthly$qin_obs, qin_sim_obs_monthly$qin_sim); mp_obs_monthly
mtext(paste0("KGE=", round(mp_obs_monthly$kge, 2), ", ", 
             "KGEr=", round(mp_obs_monthly$kger, 2), ", ", 
             "KGEb=", round(mp_obs_monthly$kgeb, 2), ", ", 
             "KGEg=", round(mp_obs_monthly$kgeg, 2)), 
      3, line=-2, col="black", cex=1)
dev.off()
# yearly average
# hype
# plot: simulations/observations over time
png(file.path("outputs_big_lakes", lake$Hylak_id, paste0(paste(lake$Hylak_id, "qin_sim_hype_yearly", sep="_"), ".png")), width=15, height=5, units="in", res=300)
par(mar=c(4.1, 5, 4.1, 5)) # bottom, left, top, right
layout(matrix(c(1,1,2), nrow=1, ncol=3, byrow=TRUE))
plot(qin_sim_hype_yearly$year, qin_sim_hype_yearly$qin_sim, type="l", col="black", lwd=2,  
     xlab="", ylab=expression("Q [m"^3*" s"^-1*"]"), 
     ylim=c(0, ceiling(max(qin_sim_hype_yearly$qin_sim, qin_sim_hype_yearly$qin_hype))), 
     cex.lab=1.5, cex.axis=1.5, cex=1.5)
lines(qin_sim_hype_yearly$year, qin_sim_hype_yearly$qin_hype, col="#003366", lwd=2)
legend("topleft", legend=c("Simulations", "HYPE"), col=c("black","#003366"), lwd=2, cex=1.4, bty="n")
# plot: simulations vs observations 
plot(qin_sim_hype_yearly$qin_sim, qin_sim_hype_yearly$qin_hype, type="p", col="#003366", pch=19,
     xlab=expression("Simulated Q [m"^3*" s"^-1*"]"),
     ylab=expression("HYPE Q [m"^3*" s"^-1*"]"),
     cex.lab=1.5, cex.axis=1.5, cex=1.5)
abline(0, 1, col="black")
mp_hype_yearly <- model_performance(qin_sim_hype_yearly$qin_hype, qin_sim_hype_yearly$qin_sim); mp_hype_yearly
mtext(paste0("KGE=", round(mp_hype_yearly$kge, 2), ", ", 
             "KGEr=", round(mp_hype_yearly$kger, 2), ", ", 
             "KGEb=", round(mp_hype_yearly$kgeb, 2), ", ", 
             "KGEg=", round(mp_hype_yearly$kgeg, 2)), 
      3, line=-2, col="black", cex=1)
dev.off()
# obs
# plot: simulations/observations over time
png(file.path("outputs_big_lakes", lake$Hylak_id, paste0(paste(lake$Hylak_id, "qin_sim_obs_yearly", sep="_"), ".png")), width=15, height=5, units="in", res=300)
par(mar=c(4.1, 5, 4.1, 5)) # bottom, left, top, right
layout(matrix(c(1,1,2), nrow=1, ncol=3, byrow=TRUE))
plot(qin_sim_obs_yearly$year, qin_sim_obs_yearly$qin_sim, type="l", col="black", lwd=2,  
     xlab="", ylab=expression("Q [m"^3*" s"^-1*"]"), 
     ylim=c(0, ceiling(max(qin_sim_obs_yearly$qin_sim, qin_sim_obs_yearly$qin_obs))), 
     cex.lab=1.5, cex.axis=1.5, cex=1.5)
lines(qin_sim_obs_yearly$year, qin_sim_obs_yearly$qin_obs, col="#990000", lwd=2)
legend("topleft", legend=c("Simulations", "Observations"), col=c("black","#990000"), lwd=2, cex=1.4, bty="n")
# plot: simulations vs observations 
plot(qin_sim_obs_yearly$qin_sim, qin_sim_obs_yearly$qin_obs, type="p", col="#990000", pch=19,
     xlab=expression("Simulated Q [m"^3*" s"^-1*"]"),
     ylab=expression("Observed Q [m"^3*" s"^-1*"]"),
     cex.lab=1.5, cex.axis=1.5, cex=1.5)
abline(0, 1, col="black")
mp_obs_yearly <- model_performance(qin_sim_obs_yearly$qin_obs, qin_sim_obs_yearly$qin_sim); mp_obs_yearly
mtext(paste0("KGE=", round(mp_obs_yearly$kge, 2), ", ", 
             "KGEr=", round(mp_obs_yearly$kger, 2), ", ", 
             "KGEb=", round(mp_obs_yearly$kgeb, 2), ", ", 
             "KGEg=", round(mp_obs_yearly$kgeg, 2)), 
      3, line=-2, col="black", cex=1)
dev.off()

# clean the environment
rm(list=ls()[!ls()%in%c("lakes", "qtot", "qg", "dis", "lat", "lon", "date", "fdir", "lat_fdir", "lon_fdir", "slope",
                        "nearest_f1", "nearest_f2", "nearest_f3", "nearest_f4", "grid_coords", "generate_grids", "generate_grids_from_center_coord",
                        "areakm2lat", "fdir_grid", "upstream_grids_of_a_grid", "upstream_grids_of_a_lake", "monthly_average", "yearly_average", "model_performance")]) 