
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
# 
# http://www.rstudio.com/shiny/
#

library(shiny)

shinyUI(pageWithSidebar(
  
  # Application title
  headerPanel("Simon's Painful Data"),
  
  # Sidebar with a slider input for number of bins
  sidebarPanel(
    conditionalPanel(
      'input.tab == "Slice Series" || input.tab == "Single Slice"',
      sliderInput("low",
                  "Lower Threshold",
                  min = 0,
                  max = 15,
                  step=0.1,
                  value = 2)),
    conditionalPanel(
      'input.tab == "Slice Series" || input.tab == "Single Slice"',
      sliderInput("high",
                  "Upper Threshold",
                  min = 0,
                  max = 15,
                  step=0.1,
                  value = 20)),
    selectInput("statistic", "Statistic to display:", choices=c()),
                #choices=c("Neonatal" = "Neonatal", 
                #          "Sex" = "mouse.gender")),
    conditionalPanel(
      'input.tab == "Slice Series"',
      selectInput("dimension", "Dimension to display:",
                  choices=c("coronal" = 2, 
                            "sagittal" = 1, 
                            "axial" = 3))),
    conditionalPanel(
      'input.tab == "Slice Series"',
      sliderInput("begin",
                  "First slice",
                  min=1, max=340, value=20)),
    conditionalPanel(
      'input.tab == "Slice Series"',
      sliderInput("end",
                  "Last slice",
                  min=-340, max=-1, value=-20)),
    conditionalPanel(
      'input.tab == "Slice Series"',
      sliderInput("rows", "Number of rows",
                  min=1, max=10, value=5)),
    conditionalPanel(
      'input.tab == "Slice Series"',
      sliderInput("columns", "Number of columns",
                  min=1, max=10, value=4)),
    conditionalPanel(
      'input.tab == "Single Slice"',
      checkboxInput("updatePlot", "Update plot", TRUE)
    ),
    conditionalPanel(
      'input.tab == "Single Slice" || input.tab == "Volumes"',
      selectInput("xvar", "X axis on plot", choices=c("Adult",
                                              "Neonatal",
                                              "treatment",
                                              "sex"))
    ),
    conditionalPanel(
      'input.tab == "Single Slice" || input.tab == "Volumes"',
      selectInput("colour", "Colour on plot", choices=c("Adult",
                                              "Neonatal",
                                              "treatment",
                                              "sex"))
    ),
    conditionalPanel(
      'input.tab == "Single Slice" || input.tab == "Volumes"',
      selectInput("fill", "fill on plot", choices=c("Adult",
                                              "Neonatal",
                                              "treatment",
                                              "sex"))
    ),

    conditionalPanel(
      'input.tab == "Single Slice" || input.tab == "Volumes"',
      selectInput("graphType", "Graph type", 
                choices=c("boxplot", "point" = "jitter")))
    
  ),
  
  # Show a plot of the generated distribution
  mainPanel(
    tabsetPanel(
      id="tab",
      tabPanel("Slice Series", plotOutput("seriesPlot", height="800px")),
      tabPanel("Single Slice", 
               fluidRow(
                 column(7,
                        plotOutput("coronalPlot", height="500px", click="plot_click")),
                 column(5,
                        plotOutput("axialPlot", height="500px", click="click_axial"))),
               fluidRow(
                 column(8,
                        plotOutput("sagittalPlot", click="click_sagittal")),
                 column(4,
                        plotOutput("graphPlot"))
               ),
               #sliderInput("slice", "Slice:",
              #             min=1, max=340, value=120, width="100%"),
      fluidRow(

                 column(12,
                        verbatimTextOutput("summaryText")))),
      tabPanel("Volumes", 
               fluidRow(
                 DT::dataTableOutput(outputId="volumesTable")),
               fluidRow(
                 plotOutput("volumesPlot")
                 #column(6, 
                  #      plotOutput("volumesPlot2"))
               )
               )
               #),
      
    )
  )
))

