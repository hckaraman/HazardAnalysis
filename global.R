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

drv = dbDriver("PostgreSQL")

server_data <- yaml.load_file("./data/server.yaml")

con = dbConnect(
  RPostgres::Postgres(),
  dbname = server_data$database,
  host = server_data$host,
  port = server_data$port,
  user = server_data$user,
  password = server_data$password
)

options(digits = 2)



station_file <- './data/stations.geojson'
river_file <- './data/river.geojson'

station <- "03-13"
query = str_interp("SELECT * FROM discharge where station = '${station}';")
# query = "select * from discharge where station ='03-13'"
data = dbGetQuery(con, query)
data <- rapply(object = data, f = round, classes = "numeric", how = "replace", digits = 6) 

stations <-  rgdal::readOGR(station_file)
rivers <-  rgdal::readOGR(river_file)

df.st <- unique(stations$Station)


template <- readRDS('./data/template.rds')


stats <- function(data,station){
  
  data.1 <- as.tibble(data) %>% select(date,discharge)
  data.1$date <- as.Date(data.1$date)
  data.1$discharge <- as.numeric(data.1$discharge)
  data.1 <- na.omit(data.1)
  
  df=template[1:length(data.1$discharge),]
  
  df$STATION_NUMBER <- station
  df$Date <- data.1$date
  df$Value <- data.1$discharge
  return(df)
}

data.1 <- stats(data,station)
start_year <- year(data.1$Date[1])

flow_stats <- calc_daily_stats(data.1, 
                 start_year = start_year,ignore_missing = TRUE )

flow_plot <- plot_daily_stats(data.1, 
                 start_year = start_year,ignore_missing = TRUE)


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
