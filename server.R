require("shiny")
library(leaflet)
library(vroom)
library(shinycssloaders)
library(highcharter)
library(xts)
library(leaflet.extras)
library(shinybusy)
library(leaflet.extras2)

options(digits=6)

function(input, output, session){

  reactive_objects=reactiveValues()
  # reactive_objects$result <-  v("03-13") 
  # reactive_objects$station <- "03-13"
  reactive_objects$lon <- 37
  reactive_objects$lat <- 38
  reactive_objects$df <- initial.result
 
  map = leaflet::createLeafletMap(session, 'map')
  
  session$onFlushed(once = T, function() {
    
    
    content <- paste(sep = "<br/>",
                     "<b><a href='http://www.samurainoodle.com'>Samurai Noodle</a></b>",
                     "606 5th Ave. S",
                     "Seattle, WA 98138"
    )
    
    
    output$map <- leaflet::renderLeaflet({
      # buildMap(sites=prof_sites, plot_polys=TRUE, au_poly=lake_aus)
      
      l  <- leaflet() %>% 
        addProviderTiles(providers$Esri.WorldImagery, group = "Esri Sattelite") %>%
        addProviderTiles(providers$OpenStreetMap.HOT, group = "OSM HOT") %>%
        addProviderTiles(providers$Stamen.TonerLite, group = "Toner Lite") %>%
        setView(lng = 35, lat = 35, zoom = 6) %>%
        setMaxBounds( lng1 = 25
                      , lat1 = 35
                      , lng2 = 45
                      , lat2 = 45
        ) %>%
        addPolylines(data=rivers,weight = 1,opacity = 0.5)
      
      l %>%
        addWMSTiles("https://ogcie.iblsoft.com/metocean/wms",
                    layers = "foreground-lines",
                    options = WMSTileOptions(format = "image/png", transparent = TRUE)
                    # group = "Lines"
        ) %>%
        addWMSTiles( "http://194.163.139.124:8600/geoserver/hazard/wms",
        # addWMSTiles( "http://172.17.0.6:8080/geoserver/hazard/wms",
                     layers = "Q500_Depth",
                     options = WMSTileOptions(format = "image/png", transparent = TRUE),
                     group = "Depth"
        ) %>%
      # addControl(paste0("<img src=","http://194.163.139.124:8600/geoserver/hazard/ows?service=WMS&request=GetLegendGraphic&format=image/png&width=20&height=20&layer=Q500_Depth&style=depth&", ">"),position = "bottomleft") %>%
        addWMSLegend(
          uri = paste0(
            "http://194.163.139.124:8600/geoserver/hazard/ows?service=WMS&request=GetLegendGraphic&format=image/png&width=20&height=20&layer=Q500_Depth&style=depth&"
          ),position = "bottomleft",
        ) %>%
      addMeasure(
          position = "topleft",
          primaryLengthUnit = "meters",
          primaryAreaUnit = "hectares",
          secondaryAreaUnit = "sqmeters",
        ) %>%
        leaflet.extras::addResetMapButton() %>%
        addLayersControl(
          baseGroups = c("Esri Sattelite", "OSM HOT", "Toner Lite"),
          overlayGroups = c("Depth"),
          options = layersControlOptions(collapsed = T),
          position = c( "bottomright"),
        ) %>%
        leaflet.extras:: addDrawToolbar(
          targetGroup = "draw",
          editOptions = editToolbarOptions(selectedPathOptions = selectedPathOptions()))  %>%
          leaflet.extras::addStyleEditor() %>%
        leaflet.extras2::addWMS(baseUrl = "http://194.163.139.124:8600/geoserver/hazard/wms",
               layers = c("Q500_Depth"),
               # layers = c("OSM-Overlay-WMS"),
               group = "wmsgroup",
               options = leaflet::WMSTileOptions(
                 transparent = TRUE,
                 format = "image/png",
                 info_format = "text/html",
                 tiled = FALSE
               ))
        # addLayersControl(baseGroups = "base",
        #                  # overlayGroups = c("wmsgroup"))
        #                  overlayGroups = c("Q500_Depth"))
        # 
  

    })
    
    output$earth <- renderHighchart({
      
      station <- reactive_objects$station
      # result <- v(station)
      result <- reactive_objects$df
      result <- result[result$Hazard == 'EARTHQUAKE',]
      
      
      highchart() %>%
        hc_add_series(data = result, name = "", type = "line", hcaes(x = IM, y = MAF)) %>%
        hc_title(text = paste("Earthquake HAZARD CURVE at: ",reactive_objects$lon,sep='')) %>%
        hc_yAxis(title = list(text = "MAF [1/TR] ")) %>%
        hc_xAxis(title = list(text = "MEAN IM")) 
    })
    
    output$flood <- renderHighchart({
      
      station <- reactive_objects$station
      # result <- v(station)
      result <- reactive_objects$df
      result <- result[result$Hazard == 'FLOOD',]
      
      
      highchart() %>%
        hc_add_series(data = result, name = "", type = "line", hcaes(x = IM, y = MAF)) %>%
        hc_title(text = paste("FLOOD HAZARD CURVE at: ",reactive_objects$lon,sep='')) %>%
        hc_yAxis(title = list(text = "MAF [1/TR] ")) %>%
        hc_xAxis(title = list(text = "MEAN IM")) 
    })
    
    output$wind <- renderHighchart({
      
      station <- reactive_objects$station
      # result <- v(station)
      result <- reactive_objects$df
      result <- result[result$Hazard == 'WIND',]
      
      
      highchart() %>%
        hc_add_series(data = result, name = "", type = "line", hcaes(x = IM, y = MAF)) %>%
        hc_title(text = paste("WIND HAZARD CURVE at: ",reactive_objects$lon,sep='')) %>%
        hc_yAxis(title = list(text = "MAF [1/TR] ")) %>%
        hc_xAxis(title = list(text = "MEAN IM")) 
    })
    
    
    output$avalance <- renderHighchart({
      
      station <- reactive_objects$station
      # result <- v(station)
      result <- reactive_objects$df
      result <- result[result$Hazard == 'AVALANCHE',]
      
      
      highchart() %>%
        hc_add_series(data = result, name = "", type = "line", hcaes(x = IM, y = MAF)) %>%
        hc_title(text = paste("AVALANCHE HAZARD CURVE at: ",reactive_objects$lon,sep='')) %>%
        hc_yAxis(title = list(text = "MAF [1/TR] ")) %>%
        hc_xAxis(title = list(text = "MEAN IM")) 
    })
    
    output$land <- renderHighchart({
      
      station <- reactive_objects$station
      # result <- v(station)
      result <- reactive_objects$df
      result <- result[result$Hazard == 'LANDSLIDE',]
      
      
      highchart() %>%
        hc_add_series(data = result, name = "", type = "line", hcaes(x = IM, y = MAF)) %>%
        hc_title(text = paste("LANDSLIDE HAZARD CURVE at: ",reactive_objects$lon,sep='')) %>%
        hc_yAxis(title = list(text = "MAF [1/TR] ")) %>%
        hc_xAxis(title = list(text = "MEAN IM")) 
    })
    
    output$rock <- renderHighchart({
      
      station <- reactive_objects$station
      # result <- v(station)
      result <- reactive_objects$df
      result <- result[result$Hazard == 'ROCKFALL',]
      
      
      highchart() %>%
        hc_add_series(data = result, name = "", type = "line", hcaes(x = IM, y = MAF)) %>%
        hc_title(text = paste("ROCKFALL HAZARD CURVE at: ",reactive_objects$lon,sep='')) %>%
        hc_yAxis(title = list(text = "MAF [1/TR] ")) %>%
        hc_xAxis(title = list(text = "MEAN IM")) 
    })
    
    output$fire <- renderHighchart({
      
      station <- reactive_objects$station
      # result <- v(station)
      result <- reactive_objects$df
      result <- result[result$Hazard == 'WILDFIRE',]
      
      
      highchart() %>%
        hc_add_series(data = result, name = "", type = "line", hcaes(x = IM, y = MAF)) %>%
        hc_title(text = paste("WILDFIRE HAZARD CURVE at: ",reactive_objects$lon,sep='')) %>%
        hc_yAxis(title = list(text = "MAF [1/TR] ")) %>%
        hc_xAxis(title = list(text = "MEAN IM")) 
    })
    
    
    output$caption <- renderText({ paste("lon : ",reactive_objects$lon," lat : ",reactive_objects$lat, sep='') })

    
    output$table <- DT::renderDataTable({
      df <- reactive_objects$df
      DT::datatable(df,  extensions = 'Buttons',options = list(dom = 'Blfrtip',
                                                               buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                                                               lengthMenu = list(c(10,25,50,-1),
                                                                                 c(10,25,50,"All"))))
    },
    options = 
      list(sPaginationType = "two_button")
    )
    
    
    
  })


  
  observeEvent(input$gobutton, {
    # showModal(modalDialog("Doing a function", footer=NULL))
    show_modal_spinner(spin ='semipolar', text = "Please Wait", color='#d11c19') # show the modal window
    print("started")
    lon <- as.numeric(reactive_objects$lon)
    lat <- as.numeric(reactive_objects$lat)
    # lon <- 29.1357
    # lat <- 39.2663
    res <- result(lon,lat)
    # df <- evaluation_hazard(reactive_objects$lon,reactive_objects$lat, name = 'Atlar', sector = 'Dolasiyor')
    reactive_objects$df <- res
    print("ended")
    remove_modal_spinner() # remove it when done
  })
  

  observeEvent(input$map_click, {
    click = input$map_click
    leafletProxy('map')%>%addMarkers(lng = click$lng, lat = click$lat, layerId = 'id')
    reactive_objects$lon <- click$lng
    reactive_objects$lat <- click$lat
  })
  
  observeEvent(input$map_marker_click, {
    leafletProxy("map") %>%
      removeMarker(input$map_marker_click$id)
  })
  
  observeEvent(input$auto1, {
    # map <- leafletProxy("map")
    df.1 <- stations[which(stations$Station == input$auto1),]
    print(input$auto1)
    
    leafletProxy("map") %>%
    setView(
      lng = df.1$lon,
      lat = df.1$lat,
      zoom = 16) 
  })
  
}

