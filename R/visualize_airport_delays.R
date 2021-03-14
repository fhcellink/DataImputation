
#'@title Plot mean delay
#'@description visualizes the mean delay of flights for different airports by longitude and latitude.
#'@field no argument
#'@examples
#' visualize_airport_delays()
#'@export 
visualize_airport_delays = function(){
  requireNamespace("nycflights13")
  requireNamespace("dplyr")
  library(dplyr)
  library(ggplot2)
  flights <- nycflights13::flights
  airports <- nycflights13::airports
  associatedData <- inner_join(flights,airports,by = c("dest" = "faa"))
  
  subsettedData <- dplyr::select(associatedData,dest,lat,lon,arr_delay)
  
  meansWithLatLon <- subsettedData %>%
    group_by(dest,lat,lon) %>%
    summarize(Mean = mean(arr_delay,na.rm = TRUE))
  
  p <- ggplot(meansWithLatLon, aes(x = lat, y = lon, label = dest)) +
    geom_point()  + theme_bw()
  
  return(p)
}

# visualize_airport_delays()
