library(shiny)
library(leaflet)
library(tidyverse)
library(vroom)
library(shinycssloaders)
library(highcharter)
library(stringr)
library(xts)
library(dqshiny)
library(shinyBS)


fluidPage(
  list(
    tags$head(
      HTML('<link rel="icon" href="./logo.png"
                type="image/png" />'))),
  
navbarPage(
  # title=div(img(src="./logo"), "My Title in the Navbar"),
  "Face Tech", id="nav",
           tabPanel("Interactive map",
                    div(class="outer",
                        
                        tags$head(
                          # Include our custom CSS
                          includeCSS("styles.css"),
                          includeScript("gomap.js"),
                          # tags$head(includeHTML(("analytics.html"))),
                        ),
                        # If not using custom CSS, set height of leafletOutput to a number instead of percent
                        leafletOutput("map", width="100%", height="100%"),
                        
                        # Shiny versions prior to 0.11 should use class = "modal" instead.
                        
                        absolutePanel(
                          id = "controls",
                          class = "panel panel-default",
                          fixed = TRUE,
                          draggable = TRUE,
                          top = 80,
                          left = "auto",
                          right = 20,
                          width = 700,
                          height = "auto",
                          fluidRow(style='height:100%',

                          h2("Hazard Analysis"),
                          # h2("DSI Flow Data"),
                          # autocomplete_input("auto1", "Station No:", placeholder = "Search for Station", df.st, max_options = 10),
                          
                          fluidRow(column(
                            9,
                          h5(verbatimTextOutput("caption"))),
                          column(3,
                          )),
                          # textOutput("text"),
                          # fluidRow(textOutput("verb")),
                          
                          
                          bsCollapse(id = "collapseExample",
                                     multiple = FALSE,
                                     bsCollapsePanel(
                                       "Show Plot",
                                       "",
                                       tabsetPanel(
                                         id = "plot_tabs",
                                         # autocomplete_input("auto1", "Unnamed:", stations, max_options = 10),
                                         tabPanel("Earthquake",
                                                  fluidRow(column(
                                                    12, h4(""),
                                                    withSpinner(highchartOutput(outputId = "earth", height =
                                                                                  600)),
                                                  ))),
                                         tabPanel("Flood",
                                                  fluidRow(column(
                                                    12,
                                                    withSpinner(highchartOutput(outputId = "flood", height =
                                                                                  600)),
                                                  ))),
                                         tabPanel("Wind",
                                                  fluidRow(column(
                                                    12,
                                                    withSpinner(highchartOutput(outputId = "wind", height =
                                                                                  600)),
                                                  ))),
                                         tabPanel("Avalance",
                                                  fluidRow(column(
                                                    12, h4(""),
                                                    withSpinner(highchartOutput(outputId = "avalance", height =
                                                                                  600)),
                                                  ))),
                                         tabPanel("Rock Fall",
                                                  fluidRow(column(
                                                    12, h4(""),
                                                    withSpinner(highchartOutput(outputId = "rock", height =
                                                                                  600)),
                                                  ))),
                                         tabPanel("Wildfire",
                                                  fluidRow(column(
                                                    12, h4(""),
                                                    withSpinner(highchartOutput(outputId = "fire", height =
                                                                                  600)),
                                                  ))),
                                         tabPanel("Landslide",
                                                  fluidRow(column(
                                                    12,
                                                    withSpinner(highchartOutput(outputId = "land", height =
                                                                                  600)),
                                                  ))),
                                         tabPanel("Table",
                                                  fluidRow(column(
                                                    12,
                                                    withSpinner(DT::dataTableOutput(outputId = "table")),
                                                  )))
                                       ),
                                       actionButton("gobutton", "Run", style='padding:4px; height:100%')
                                       # bsButton("actTwo", label = "Click me if you dare!", icon = icon("ban"))
                                     ))
                          )
                        ),
                        
                        tags$div(id="cite",
                                 'Data prepared for Face Tech by ', tags$a(href="https://github.com/hckaraman", "Cagri Karaman")
                        )
                    )
           ),
  
  # tabPanel("Data explorer",
  #          shinycssloaders::withSpinner(DT::dataTableOutput("pond_stat", height="600px"),size=2, color="#0080b7")
  # ),
  # 
  conditionalPanel("false", icon("crosshair"))
))
