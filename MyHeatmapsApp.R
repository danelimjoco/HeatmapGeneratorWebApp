library(shiny)
library(shinythemes)
library(readxl)

######## UI ########

ui <- fluidPage(titlePanel(h1("RNA-Seq Shiny App!",align = "center")), theme = shinytheme("cerulean"),
                
                           fluidRow(
                             column(4, align="center", offset=4,
                                    h3(fileInput("countsFile", "1. Choose count data file (xlsx):",
                                              multiple = TRUE,
                                              accept = c("text/csv",
                                                         "text/comma-separated-values,text/plain",
                                                         ".csv",
                                                         ".xlsx")), 
                                       )
                             )
                           ),
                           conditionalPanel(condition = "output.fileUploaded", 
                                            h3("Preview of uploaded data"), 
                                            tableOutput('contents'), 
                                            align="center"
                           ),
                           br(), #some space 
                           conditionalPanel(condition = "output.fileUploaded", 
                                            h3(fileInput("metadataFile", "2. Choose metadata file (xlsx):",
                                                         multiple = TRUE,
                                                         accept = c("text/csv",
                                                                    "text/comma-separated-values,text/plain",
                                                                    ".csv",
                                                                    ".xlsx")),
                                             ),   
                                            align="center"
                           ),
                           br(),
                           br(), #some space
                           conditionalPanel(condition = "output.metadataUploaded", 
                                            h3("3. Select Filters"), 
                                            
                                            # Input: Simple integer interval ----
                                            sliderInput("integer", "Integer:",
                                                        min = 0, max = 1000,
                                                        value = 500),
                                            
                                            # Input: Decimal interval with step value ----
                                            sliderInput("decimal", "Decimal:",
                                                        min = 0, max = 1,
                                                        value = 0.5, step = 0.1),
    
                                            align="center",
                                            br(), #some space 
                                            br(), #some space 
                                            br(), #some space 
                                            br(), #some space 
                                            br(), #some space 
                                            br(), #some space 
                                            br(), #some space 
                                            br(), #some space 
                                            br(), #some space 
                                            br(), #some space 
                           ),
                          br(), #some space 
                          br(), #some space 
                  
) # fluidPage




######## SERVER ########

server <- function(input, output) {
  
  # Count data
  output$contents <- renderTable({
    inFile <- input$countsFile
    read_excel(inFile$datapath, range = "A1:D6")
  })
  
  # Metadata
  output$contents2 <- renderTable({
    inFile <- input$metadataFile
    read_excel(inFile$datapath)
  })
  
  
  output$fileUploaded <- reactive({
    return(!is.null(input$countsFile))
  })
  
  output$metadataUploaded <- reactive({
    return(!is.null(input$metadataFile))
  })
  
  outputOptions(output, 'fileUploaded', suspendWhenHidden=FALSE)
  outputOptions(output, 'metadataUploaded', suspendWhenHidden=FALSE)
  
} # server

######## CREATE APP ########
shinyApp(ui = ui, server = server)

