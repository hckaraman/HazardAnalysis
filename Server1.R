library(shiny)
library(ggplot2)
library(DT)
library(tidyverse)
library(reticulate)
library(plotly)
library(openxlsx)
library(vroom)
library(shinycssloaders)



HAZARD = c("EARTHQUAKE",
           "FLOOD",
           "WIND",
           "AVALANCHE",
           "LANDSLIDE",
           "ROCKFALL",
           "WILDFIRE")

TECH_HAZARD = c("TECHFIRE", "TECHBLAST")

use_python(python = "./venv/Scripts/python.exe", required = TRUE)
source_python("evaluation_hazard.py")

lon <- 27.37124
lat <- 39.11272
res <- evaluation_hazard(lon, lat, name = 'Atlar', sector = 'Dolasiyor')
res.earh <- res[res$Hazard == 'EARTHQUAKE',]
res.flood <- res[res$Hazard == 'FLOOD',]
res.wind <- res[res$Hazard == 'WIND',]
res.avalance <- res[res$Hazard == 'AVALANCHE',]
# length(res.avalance$IM)
res.landslide <- res[res$Hazard == 'LANDSLIDE',]
res.rackfall <- res[res$Hazard == 'ROCKFALL',]
res.wildfire <- res[res$Hazard == 'WILDFIRE',]

