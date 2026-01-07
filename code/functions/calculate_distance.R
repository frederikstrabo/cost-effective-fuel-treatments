##########################################    Calculate Distance    ##########################################

################# Purpose: Calculates distance from grid to different assets


calculate_distance <- function (SpatialRD_Grids, Assets){
  
  index <- st_nearest_feature(x = SpatialRD_Grids, y = Assets)
  
  # slice based on the index
  poly2 <- Assets %>% slice(index)
  
  # calculate distance between polygons
  grid_dist <- st_distance(x = SpatialRD_Grids, y = poly2, by_element = TRUE)

  return(as.numeric(grid_dist/1000))
  
}

calculate_distance_LAT <- function (SpatialRD_Grids, Assets){
  
  index <- st_nearest_feature(x = SpatialRD_Grids, y = Assets)
  
  # slice based on the index
  poly2 <- Assets %>% slice(index)
  
  # calculate distance between polygons
  grid_dist <- st_distance(x = SpatialRD_Grids, y = poly2, by_element = TRUE)
  
  return(as.numeric(grid_dist/1000))
  
}