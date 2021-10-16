library(RSQLite)
library(RPostgreSQL)
library(DBI)
library(yaml)
library(fasstr)
library(stringr)
library(tidyverse)
library(lubridate)
library(ggplot2)
library(reticulate)

station_file <- './data/stations.geojson'
river_file <- './data/river.geojson'



stations <-  rgdal::readOGR(station_file)
rivers <-  rgdal::readOGR(river_file)

HAZARD = c("EARTHQUAKE",
           "FLOOD",
           "WIND",
           "AVALANCHE",
           "LANDSLIDE",
           "ROCKFALL",
           "WILDFIRE")

TECH_HAZARD = c("TECHFIRE", "TECHBLAST")
print(getwd())
# use_python(python = "./venv/Scripts/python.exe", required = TRUE)
use_python(python = "./linuxvenv/bin/python", required = TRUE)
source_python("./evaluation_hazard.py")


result <- function(lon,lat) {
  res <- evaluation_hazard(lon, lat, name = 'Atlar', sector = 'Dolasiyor')
  return(res)
}

# initial.result <- result(29.1357,39.2663)
initial.result  <- read_rds('./data/res.rds')
# write_rds(initial.result,'./data/res.rds')