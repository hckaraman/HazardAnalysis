require("shiny")
library(leaflet)
library(vroom)
library(shinycssloaders)
library(highcharter)
library(xts)
library(leaflet.extras)
library(shinybusy)
options(digits=6)




function(input, output, session){
  
  
  v <- function(station) {
    # df.1 <- df[which(df$Name == id),]
    query = str_interp("SELECT * FROM discharge where station = '${station}';")
    df = dbGetQuery(con, query)
    x <- df$discharge
    x1 <- replace(x,is.na(x),0)
    df$Dischage <- as.numeric(x1)
    df$date <- as.Date(df$date)
    # df$Month <- months(df$date)
    # df$Year <- format(df$date,format="%y")
    df <- rapply(object = df, f = round, classes = "numeric", how = "replace", digits = 6) 

    
    dfm <- df %>%
      select(date,Dischage) %>%
      group_by(lubridate::month(date),lubridate::year(date) ) %>%
      summarise(month_average = mean(Dischage))
    names(dfm) <- c("month","year","Dischage")
    dfm$date <- with(dfm, sprintf("%d-%02d-01", year, month))
    dfmx = xts(dfm$Dischage, order.by=as.Date(dfm$date))
    # dfm <- aggregate( X2 ~ Month + Year , df , mean )
    dfx = xts(df$Dischage, order.by=as.Date(df$date))
    
    
    x <-na.omit(df$discharge)
    x <- as.numeric(x)
    x.old <- x
    x <- sort(x)
    x.zero.index <- which(x==0)
    nzeros <- length(x.zero.index)
    ind <- match(x.old, x)
    n <- length(x)
    dc <- rep(NA, n)
    dc[1:n] <- sapply(1:n, function(j,y) {
      dc[j] <- length( which(y >= y[j]) )
    }, y = x)
    dc <- dc/n
    fdr<-data.frame(x=dc,y=x)
    
    
    data.1 <- stats(df,station)
    start_year <- year(data.1$Date[1])
    
    flow_stats <- calc_daily_stats(data.1, 
                                   start_year = start_year,ignore_missing = TRUE )
    
    flow_plot <- plot_daily_stats(data.1, 
                                  start_year = start_year,ignore_missing = TRUE)
    
    
    res <- list()
    res$df <- df
    res$dfx <- dfx
    res$dfmx <- dfmx
    res$fdr <- fdr
    res$flow_stats <- flow_stats
    res$flow_plot <- flow_plot
    res$lat <- 37
    res$lon <- 38
    # print(lat)
    return(res)
  }
  
  
  reactive_objects=reactiveValues()
  reactive_objects$result <-  v("03-13") 
  reactive_objects$station <- "03-13"
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
      # leaflet() %>%
        # addProviderTiles("Esri.WorldImagery") %>%
        setView(lng = 35, lat = 35, zoom = 6) %>%
        # addMarkers(
        #   data = stations,
        #   # label = paste0(pond_point$Name),
        #   popup = paste0(
        #     "<b>Name: </b>"
        #     , stations$Name
        #     , "<br>"
        #     ,"<b>Station No: </b>"
        #     , stations$Station
        #     , "<br>"
        #     ,"<b>Elevation : </b>"
        #     , stations$Elevation ," m"
        #     , "<br>"
        #     ,"<b>Basin Area: </b>"
        #     , stations$Basin_Area , "m"
        #   ),
        #   # labelOptions = labelOptions(noHide = F),
        #   layerId = ~Station,
        #   clusterOptions = markerClusterOptions()) %>%
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
      addMeasure(
          position = "topleft",
          primaryLengthUnit = "meters",
          primaryAreaUnit = "hectares",
          secondaryAreaUnit = "sqmeters",
        ) %>%
        leaflet.extras::addResetMapButton() %>%
        addLayersControl(
          baseGroups = c("Esri Sattelite", "OSM HOT", "Toner Lite"),
          options = layersControlOptions(collapsed = T),
          position = c( "bottomright"),
        ) %>%
        leaflet.extras:: addDrawToolbar(
          targetGroup = "draw",
          editOptions = editToolbarOptions(selectedPathOptions = selectedPathOptions()))  %>%
          leaflet.extras::addStyleEditor()
      
      
      
      
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
    
    
    output$flowdur <- renderHighchart({
      
      station <- reactive_objects$station
      # result <- v(station)
      fdr <- reactive_objects$result$fdr
      fdr$x <- fdr$x * 1e2
      
      
      
      # reactive_objects=reactiveValues()
      # reactive_objects$result <-  v("03-13") 
      # reactive_objects$station <- "03-13"
      
      
      # highchart(type = "stock") %>% 
      #   hc_yAxis_multiples(create_yaxis(2, height = c(1, 1), turnopposite = TRUE)) %>% 
      #   hc_title(text = paste("Observed discharge at station : ",station)) %>%
      #   hc_add_series(result$dfx, yAxis = 0,name = "Daily Mean") %>%
      
      #   hc_add_series(result$dfmx, yAxis = 1,name = "Monthly Mean") %>%
      #   # hc_add_yAxis(nid = 1L, title = list(text = "Discharge m3/s"), relative = 4) %>%
      #   hc_xAxis(
      #     type = 'datetime') %>%
      #   hc_legend(enabled = TRUE) %>%
      #   hc_tooltip(
      #     crosshairs = TRUE,
      #     backgroundColor = "#F0F0F0",
      #     shared = TRUE, 
      #     borderWidth = 5
      #   )
      highchart() %>%
        hc_add_series(data = fdr, name = "", type = "line", hcaes(x = x, y = y)) %>%
        hc_yAxis(type = "logarithmic",title = list(text = "Discharge, m3/s")) %>%
        hc_xAxis(labels = list(format = "{value}%"),title = list(text = "% Time flow equalled or exceeded")) 
      
    })
    
    output$plot2<-renderPlot({
      
      reactive_objects$result$flow_plot
      
    })
    
    output$caption <- renderText({ paste("lon : ",reactive_objects$lon," lat : ",reactive_objects$lat, sep='') })

    
    output$table <- DT::renderDataTable({
      
      # station <- reactive_objects$Station
      # # df <- data[which(data$Station == station),]
      # query = str_interp("SELECT * FROM discharge where station = '${station}';")
      # df = dbGetQuery(con, query)
      # df$Dischage <- as.numeric(df$dischage)
      # result <- reactive_objects$result
      
      df <- reactive_objects$df
      DT::datatable(df,  extensions = 'Buttons',options = list(dom = 'Blfrtip',
                                                               buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                                                               lengthMenu = list(c(10,25,50,-1),
                                                                                 c(10,25,50,"All"))))
    },
    options = 
      list(sPaginationType = "two_button")
    )
    
    
    output$flowstats <- DT::renderDataTable({
      
      # station <- reactive_objects$Station
      # # df <- data[which(data$Station == station),]
      # query = str_interp("SELECT * FROM discharge where station = '${station}';")
      # df = dbGetQuery(con, query)
      # df$Dischage <- as.numeric(df$dischage)
      # result <- reactive_objects$result
      
      # df <- reactive_objects$result$df 
      # df <- select(df,"date","Dischage")
      # row.names(df) <- NULL
      # reactive_objects$result$flow_stats
      DT::datatable(reactive_objects$result$flow_stats,  extensions = 'Buttons',options = list(dom = 'Blfrtip',
                                                                       scrollX=T,
                                                               buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                                                               lengthMenu = list(c(10,25,50,-1),
                                                                                 c(10,25,50,"All"))))
    },
    options = 
      list(sPaginationType = "two_button")
    )
    
  })
  # 
  # observe({
  #   click <- input$map_marker_click
  #   if (is.null(click)){return()}
  #   p <- input$map_marker_click$id
  #   # siteid=site_click$id
  #   # reactive_objects$sel_mlid=siteid
  #   reactive_objects$station=p
  #   reactive_objects$result=v(p)
  #   reactive_objects$lat <- click$lat
  #   reactive_objects$lon <- click$lon
  #   print(reactive_objects$lat)
  # })
  
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
    
    # removeModal()
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
    # 
    # reactive_objects$station=input$auto1
    # reactive_objects$result=v(p)
    # print(reactive_objects$station)
    
  })
  
}

