## Function to calculate the total area of shapefile within each cell #######
## Reference: https://stackoverflow.com/questions/61867321/finding-the-total-are-for-sf-polygons-intersecting-a-grid-cell

library(sf)

shp_area_prop = function(points_coords, shp_data, thres = 0.2){
  ## points_coords: data frame of grid centre points
  ## shp_data: shapefile
  options(warn=-1)
  PROJ_LATLON = '+proj=longlat +datum=WGS84'
 
  ##### 1) create point-centred grids ####
  
  # convert points to sf 
  names(points_coords) = c('lon', 'lat') # make sure have location names
  locs_sf = points_coords%>%st_as_sf(coords = c('lon', 'lat')) %>% 
    st_set_crs(4326)
  
  # get 0.025 degrees (2500km) buffer around the centre points
  box = st_buffer(locs_sf, 2500) %>%
    st_bbox() %>%  st_as_sfc()
  
  ## calculate grid dimension
  lon_length = length(seq(min(points_coords$lon),max(points_coords$lon),0.05 ))
  lat_length = length(seq(min(points_coords$lat),max(points_coords$lat),0.05 ))
  
  grid <- st_make_grid(box, n = c(lon_length,lat_length)) 
  
  ### 2) convert shapefile to sf 
  proj4string(shp_data) = CRS(PROJ_LATLON) # assign crs
  shp_sf = st_as_sf(shp_data, crs = PROJ_LATLON)  
  
  ## 3) Assign a ID to each grid (order in the coordinates)
  
  grid_sf = grid%>%st_as_sf()%>%mutate(CellID = row_number())
  
  ## 4) Subset the grids that have overlap with shapefile region
  
  intersect_info = st_intersection(x = grid_sf, y =  shp_sf)
  
  # (spatial feature is contained in intersect_info$x)
  
  # store ID that intersects
  inter_ID = intersect_info$CellID  
  
  ## 5) compute area of the entire grid that has intersection
  grid_sf_inter = grid_sf[inter_ID, ]
  grid_total_area = grid_sf_inter %>% mutate(area = st_area(grid_sf_inter))
  
  #### 6) Finally, Calculate the intersecting area within the grid ### 
  
  grid_inter_area = intersect_info %>%group_by(CellID)%>%   # grid-specific
    summarise(geom=st_union(x))%>%  ## collate all spatial feature in that intersection
    mutate(geom = st_sfc(geom), area = st_area(geom))   # calculate actual area 
  
  ## calculate the proportion of area in the grid ##
  grid_prop = as.vector(grid_inter_area$area/grid_total_area$area)
  finalID = inter_ID[grid_prop>thres]

  
   return(list(fullgrid = grid_sf, inter_area = grid_inter_area, 
               inter_ID=inter_ID, inter_prop =grid_prop, finalID =  finalID))
}
