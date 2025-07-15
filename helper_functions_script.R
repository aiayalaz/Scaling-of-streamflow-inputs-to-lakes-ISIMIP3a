### Integration of the Global Water and Lake Sectors within the ISIMIP framework through scaling of streamflow inputs to lakes
# Ana I. Ayala. Limnology Unit, Department of Ecology and Genetics, Uppsala University, Uppsala, Sweden.
# José L. Hinostroza. Faculty of Civil Engineering, National University of Engineering, Lima, Perú.
# July 2025 

# function: nearest_f1()
# given a number, calculate the nearest to .25 or .75 (nearest number)
# nearest_f1(18.6)[1]=18.75
nearest_f1 <- function(x){
  a <- (x+0.25)*2
  b <- floor(a)
  c <- b+1
  if (abs(a-b) <= abs(a-c)){
    d <- b
  } else {
    d <- c
  }
  lb <- b/2-0.25 # lower bound (bottom or left)
  ub <-c/2-0.25 # upper bound (top or right)
  nearest <- d/2-0.25 # nearest
  return(c(nearest,lb,ub))
}

# function: nearest_f2()
# given (x,y), calculate the nearest to .25 or .75 (nearest coordinate)
# nearest_f2(62.4,18.6)[2,]=(62.25,18.75)
nearest_f2 <- function(x,y){
  # nearest, lower bound, upper bound 
  vx <- nearest_f1(x) 
  vy <- nearest_f1(y) 
  # near (X,Y)
  tl <- c(vx[2],vy[3]) 
  tr <- c(vx[3],vy[3])
  br <- c(vx[3],vy[2])
  bl <- c(vx[2],vy[2])
  # distance between (x,y) and (X,Y)
  dist1 <- dist(t(matrix(c(x,y,tl), nrow=2))) 
  dist2 <- dist(t(matrix(c(x,y,tr), nrow=2)))
  dist3 <- dist(t(matrix(c(x,y,bl), nrow=2)))
  dist4 <- dist(t(matrix(c(x,y,br), nrow=2)))
  # minimum distance
  dist_min <- min(dist1,dist2,dist3,dist4)
  if (dist_min==dist1){
    nearest <- tl
  } else if (dist_min==dist2){
    nearest <- tr
  } else if (dist_min==dist3){
    nearest <- bl
  } else if (dist_min==dist4){
    nearest <- br
  } else{
    print("error")
  }
  return(t(matrix(c(x,y,nearest,tl,tr,bl,br), nrow=2)))
}

# function: nearest_f3()
# given a number, calculate the nearest to .0 or .5 (nearest number)
# nearest_f3(18.6)[1]=18.5
nearest_f3 <- function(x){
  a <- x*2
  b <- floor(a)
  if(a==b){
    c <- b
  } else {
    c <- b+1
  }
  if (abs(a-b) <= abs(a-c)){
    d <- b
  } else {
    d <- c
  }
  lb <- b/2# lower bound (bottom or left)
  ub <-c/2 # upper bound (top or right)
  nearest <- d/2 # nearest
  return(c(nearest,lb,ub))
}
# function: nearest_f4()
# given (x,y), calculate the nearest to .0 or .5 (nearest coordinate)
# nearest_f4(62.4,18.6)[2,]=(62.5,18.5)
nearest_f4 <- function(x,y){
  # nearest, lower bound, upper bound 
  vx <- nearest_f3(x)
  vy <- nearest_f3(y)
  # near (X,Y)
  tl <- c(vx[2],vy[3])
  tr <- c(vx[3],vy[3])
  bl <- c(vx[2],vy[2])
  br <- c(vx[3],vy[2])
  # distance between (x,y) and (X,Y)
  dist1 <- dist(t(matrix(c(x,y,tl), nrow=2)))
  dist2 <- dist(t(matrix(c(x,y,tr), nrow=2)))
  dist3 <- dist(t(matrix(c(x,y,bl), nrow=2)))
  dist4 <- dist(t(matrix(c(x,y,br), nrow=2)))
  # minimum distance
  dist_min <- min(dist1,dist2,dist3,dist4)
  if (dist_min==dist1){
    nearest <- tl
  } else if (dist_min==dist2){
    nearest <- tr
  } else if (dist_min==dist3){
    nearest <- bl
  } else if (dist_min==dist4){
    nearest <- br
  } else{
    print("error")
  }
  return(t(matrix(c(x,y,nearest,tl,tr,bl,br), nrow=2)))
}

# function: grid_coords
# generate a grid every 0.5 degrees (.shp file)
grid_coords <- function(xgrid_min, ygrid_min, res=0.5){
  m <- rbind(c(xgrid_min, ygrid_min), 
             c(xgrid_min+res, ygrid_min),
             c(xgrid_min+res, ygrid_min+res), 
             c(xgrid_min, ygrid_min+res),
             c(xgrid_min, ygrid_min))
  return(m)
}

# function: generate_grids()
# generate grids where the lake is located (from the polygonal coordinates)
generate_grids <- function(lake){ 
  
  coords <- as.vector(st_bbox(lake)) # minimum and maximum values of a polygon (xmin, ymin, xmax, ymax)
  
  grid_min <- nearest_f4(coords[1],coords[2])[5,] # min coord for the grids
  xgrid_min <- grid_min[1]
  ygrid_min <- grid_min[2]
  
  grid_max <- nearest_f4(coords[3],coords[4])[4,] # max coord for the grids
  xgrid_max <- grid_max[1]
  ygrid_max <- grid_max[2]
  
  res <- 0.5 # grid resolution: 0.5degree x 0.5degree
  n <- as.integer((xgrid_max-xgrid_min)/res) # number of grids on the x-axis
  m <- as.integer((ygrid_max-ygrid_min)/res) # number of grids on the y-axis
  
  polygons <- list() 
  for(i in 1:n){
    for(j in 1:m){
      polygons <- append(polygons, list(st_polygon(list(grid_coords(xgrid_min+(i-1)*res, ygrid_min+(j-1)*res)))))
    }
  }
  polygons <- polygons %>% st_sfc() %>% st_as_sf()
  
  st_crs(polygons) <- 4326
  
  return(polygons)
}

# function: generate_grids_from_center_coord()
# generate grids (from grid center coordinates)
generate_grids_from_center_coord <- function(upstream_grids){
  res <- 0.5 # grid resolution: 0.5degree x 0.5degree
  polygons <- list()
  for(i in 1:nrow(upstream_grids)){
    polygons <- append(polygons, list(st_polygon(list(grid_coords(upstream_grids[i, "lon"]-(res/2), upstream_grids[i, "lat"]-(res/2))))))
  }
  polygons <- polygons %>% st_sfc() %>% st_as_sf()
  
  st_crs(polygons) <- 4326
  
  return(polygons)
}

# function: areakm2lat()
# grid area (km2)
areakm2lat <- function(lat){
  res <- 0.5 # grid resolution: 0.5degree x 0.5degree
  R <- 6371007 # authalic earth radius at equator (m)
  height <- res*pi/180*R # height of the cells, same value for the whole grid
  width <- (sin((lat+res/2)*pi/180)-sin((lat-res/2)*pi/180))*R # cells width
  area_km2 <- width*height/1e6 # cells area depending on latitude
  return(area_km2)
}

# function: fdir_grid()
# grids draininig to a grid (8 grids around)
fdir_grid <- function(fdir, lon_fdir, lat_fdir, lon_grid, lat_grid) {
  # left
  fdir_left_grid <- fdir[which(lon_fdir==(lon_grid-0.5)), which(lat_fdir==lat_grid)] 
  fdir_left_grid <- ifelse(is.na(fdir_left_grid)==FALSE, ifelse(fdir_left_grid!=1, FALSE, TRUE), FALSE)
  # top-left
  fdir_top_left_grid <- fdir[which(lon_fdir==(lon_grid-0.5)), which(lat_fdir==(lat_grid+0.5))] 
  fdir_top_left_grid <- ifelse(is.na(fdir_top_left_grid)==FALSE, ifelse(fdir_top_left_grid!=2, FALSE, TRUE), FALSE)
  # top
  fdir_top_grid <- fdir[which(lon_fdir==lon_grid), which(lat_fdir==(lat_grid+0.5))] 
  fdir_top_grid <- ifelse(is.na(fdir_top_grid)==FALSE, ifelse(fdir_top_grid!=3, FALSE, TRUE), FALSE)
  # top-right
  fdir_top_right_grid <- fdir[which(lon_fdir==(lon_grid+0.5)), which(lat_fdir==(lat_grid+0.5))] 
  fdir_top_right_grid <- ifelse(is.na(fdir_top_right_grid)==FALSE, ifelse(fdir_top_right_grid!=4, FALSE, TRUE), FALSE)
  # right
  fdir_right_grid <- fdir[which(lon_fdir==(lon_grid+0.5)), which(lat_fdir==lat_grid)] 
  fdir_right_grid <- ifelse(is.na(fdir_right_grid)==FALSE, ifelse(fdir_right_grid!=5, FALSE, TRUE), FALSE)
  # bottom-right
  fdir_bot_right_grid <- fdir[which(lon_fdir==(lon_grid+0.5)), which(lat_fdir==(lat_grid-0.5))] 
  fdir_bot_right_grid <- ifelse(is.na(fdir_bot_right_grid)==FALSE, ifelse(fdir_bot_right_grid!=6, FALSE, TRUE), FALSE)
  # bottom
  fdir_bot_grid <- fdir[which(lon_fdir==lon_grid), which(lat_fdir==(lat_grid-0.5))] 
  fdir_bot_grid <- ifelse(is.na(fdir_bot_grid)==FALSE, ifelse(fdir_bot_grid!=7, FALSE, TRUE), FALSE)
  # bottom-left
  fdir_bot_left_grid <- fdir[which(lon_fdir==(lon_grid-0.5)), which(lat_fdir==(lat_grid-0.5))] 
  fdir_bot_left_grid <- ifelse(is.na(fdir_bot_left_grid)==FALSE, ifelse(fdir_bot_left_grid!=8, FALSE, TRUE), FALSE)
  
  fdir_grids <- data.frame(c(lon_grid-0.5, lon_grid-0.5, lon_grid, lon_grid+0.5, lon_grid+0.5, lon_grid+0.5, lon_grid, lon_grid-0.5),
                           c(lat_grid, lat_grid+0.5, lat_grid+0.5, lat_grid+0.5, lat_grid, lat_grid-0.5, lat_grid-0.5, lat_grid-0.5),
                           c(fdir_left_grid, fdir_top_left_grid, fdir_top_grid, fdir_top_right_grid, fdir_right_grid, fdir_bot_right_grid, fdir_bot_grid, fdir_bot_left_grid))
  colnames(fdir_grids) <- c("lon","lat","fdir")
  
  selected_fdir_grids <- fdir_grids[fdir_grids$fdir==TRUE, c("lon","lat")]
  
  return(selected_fdir_grids)
}

# function: upstream_grids_of_a_grid()
# upstream grids of a grid 
upstream_grids_of_a_grid <- function(fdir, lon_fdir, lat_fdir, lon_grid, lat_grid) {
  # lake grid
  level <- 1 # level 1
  coords_fdir_grid0 <- data.frame(lon_grid, lat_grid, level)
  colnames(coords_fdir_grid0) <- c("lon","lat","level")
  coords_fdir_grid0
  
  # grids around the lake grid
  level <- level+1 # level 2
  if(nrow(fdir_grid(fdir, lon_fdir, lat_fdir, lon_grid, lat_grid))!=0) {
    coords_fdir_grids0 <- data.frame(fdir_grid(fdir, lon_fdir, lat_fdir, lon_grid, lat_grid), level)
    colnames(coords_fdir_grids0) <- c("lon","lat","level")
    dum <- TRUE
  } else {
    coords_fdir_grids0 <- data.frame()
    dum <- FALSE
  }
  COORD_FDIR_GRIDS <- rbind(coords_fdir_grid0, coords_fdir_grids0)
  
  # grids around another grid
  while(dum==TRUE) {
    level <- level+1
    
    COORD_FDIR_GRIDS0 <- data.frame()
    for (i in 1:nrow(coords_fdir_grids0)) {
      if(nrow(fdir_grid(fdir, lon_fdir, lat_fdir, coords_fdir_grids0[i, "lon"], coords_fdir_grids0[i, "lat"]))!=0) {
        coords_fdir_grids <- data.frame(fdir_grid(fdir, lon_fdir, lat_fdir, coords_fdir_grids0[i, "lon"], coords_fdir_grids0[i, "lat"]), level)
        colnames(coords_fdir_grids) <- c("lon","lat","level")
        COORD_FDIR_GRIDS0 <- rbind(COORD_FDIR_GRIDS0, coords_fdir_grids)
      }
    }
    
    if (nrow(COORD_FDIR_GRIDS0)!=0) {
      COORD_FDIR_GRIDS <- rbind(COORD_FDIR_GRIDS, COORD_FDIR_GRIDS0)
      coords_fdir_grids0 <- COORD_FDIR_GRIDS0
    } else {
      dum <- FALSE
    }
    
  }
  
  return(COORD_FDIR_GRIDS)
}

# function: upstream_grids_of_a_lake()
# upstream grids of a lake (when the lake occupies more than one grid) 
upstream_grids_of_a_lake <- function(fdir, lon_fdir, lat_fdir, lake_grids) {
  # all lake grids
  COORD_FDIR_GRIDS0 <- data.frame()
  for (i in 1:nrow(lake_grids)) {
    coords_fdir_grids0 <- upstream_grids_of_a_grid(fdir, lon_fdir, lat_fdir, lake_grids[i, "lon"], lake_grids[i, "lat"])
    COORD_FDIR_GRIDS0 <- rbind(COORD_FDIR_GRIDS0, coords_fdir_grids0)
  }
  
  # remove duplicates
  COORD_FDIR_GRIDS <- data.frame()
  coords_fdir_grids <- unique(COORD_FDIR_GRIDS0[, c("lon","lat")])
  for (i in 1:nrow(coords_fdir_grids)) {
    COORD_FDIR_GRIDS <- rbind(COORD_FDIR_GRIDS, 
                              data.frame(coords_fdir_grids[i, ], 
                                         min(COORD_FDIR_GRIDS0[COORD_FDIR_GRIDS0$lon==coords_fdir_grids[i, "lon"] & COORD_FDIR_GRIDS0$lat==coords_fdir_grids[i, "lat"], "level"])))
    
  }
  colnames(COORD_FDIR_GRIDS) <- c("lon","lat","level")
  
  return(COORD_FDIR_GRIDS)
}  

# function: monthly_average()
# calculation of the monthly average from daily data
monthly_average <- function(qin, lag) {
  # monthly average
  dt_qin <- data.table(qin)
  setnames(dt_qin, c("date","qin"))
  dt_qin[, year := year(date)]
  dt_qin[, month := month(date)]
  dt_qin_monthly <- dt_qin[, lapply(.SD, mean, na.rm=TRUE), by=.(year, month), .SDcols=c("qin")] # monthly average per year
  dt_qin_monthly <- na.omit(dt_qin_monthly[, lapply(.SD, function(x) replace(x, is.nan(x), NA))])
  # set a date
  # sometimes there is a time lag of one month
  if (lag==TRUE) {
    dt_qin_monthly[, date := as.Date(paste(ifelse(month!=1, year, year-1), ifelse(month!=1, month-1, 12), "01", sep="-"), format="%Y-%m-%d")] 
  } else {
    dt_qin_monthly[, date := as.Date(paste(year, month, "01", sep="-"), format="%Y-%m-%d")]
  }
  return(dt_qin_monthly[, c("date","qin")])
}

# function: yearly_average()
# calculation of the yearly average from daily data
yearly_average <- function(qin) {
  dt_qin <- data.table(qin)
  setnames(dt_qin, c("date","qin"))
  dt_qin[, year := year(date)]
  dt_qin_yearly <- dt_qin[, lapply(.SD, mean, na.rm=TRUE), by=.(year), .SDcols=c("qin")] # yearly average per year
  dt_qin_yearly <- na.omit(dt_qin_yearly[, lapply(.SD, function(x) replace(x, is.nan(x), NA))])
  return(dt_qin_yearly)
}

# function: model_performance()
# bias, mae, rmse, nrmse, nse
model_performance <- function(obs, sim) {
  #  bias
  bias <- mean((sim-obs), na.rm=TRUE)
  # mean absolute error
  mae <- mean(abs(obs-sim), na.rm=TRUE)
  # root mean square error
  rmse <- sqrt(mean((obs-sim)^2, na.rm=TRUE))
  # normalized root mean square error
  nrmse <- rmse/(max(obs, na.rm=TRUE)-min(obs, na.rm=TRUE))
  # correlation coef
  r <- sum((obs - mean(obs, na.rm=TRUE))*(sim - mean(sim, na.rm=TRUE)), na.rm=TRUE)/sqrt(sum((obs - mean(obs, na.rm=TRUE))^2, na.rm=TRUE)*sum((sim - mean(sim, na.rm=TRUE))^2, na.rm=TRUE))
  # nash sutcliff efficiency
  nse <- 1 - sum((obs-sim)^2, na.rm=TRUE)/sum((obs - mean(obs, na.rm=TRUE))^2, na.rm=TRUE)
  # kling-gupta efficiency 
  kge <- KGE(sim, obs, method="2012", out.type="full")$KGE.value
  kger <- KGE(sim, obs, method="2012", out.type="full")$KGE.elements[[1]]
  kgeb <- KGE(sim, obs, method="2012", out.type="full")$KGE.elements[[2]]
  kgeg <- KGE(sim, obs, method="2012", out.type="full")$KGE.elements[[3]]
  
  df <- data.frame(bias=bias, mae=mae, rmse=rmse, nrmse=nrmse, r=r, nse=nse, kge=kge, kger=kger, kgeb=kgeb, kgeg=kgeg, n=length(obs))
  return(df)
}
