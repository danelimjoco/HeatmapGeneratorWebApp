library(shiny)
library(shinythemes)
library(readxl)

######## UI ########

ui <- fluidPage(titlePanel(h1("RNA-Seq Shiny App!",align = "center")), theme = shinytheme("cerulean"),
                
                           fluidRow(
                             column(4, align="center", offset=4,
                                    h3(fileInput("countsFile", "Choose xlsx File:",
                                              multiple = TRUE,
                                              accept = c("text/csv",
                                                         "text/comma-separated-values,text/plain",
                                                         ".csv",
                                                         ".xlsx")), 
                                       )
                             )#,
                             #column(12, h3("Preview of uploaded data", align="center")
                             #)
                           ),
                           conditionalPanel(condition = "output.fileUploaded", h3("Preview of uploaded data"), tableOutput('contents'), align="center")
                  
) # fluidPage




######## SERVER ########

server <- function(input, output) {
  
  output$contents <- renderTable({
    inFile <- input$countsFile
    read_excel(inFile$datapath, range = "A1:D6")
  })
  
  output$fileUploaded <- reactive({
    return(!is.null(input$countsFile))
  })
  outputOptions(output, 'fileUploaded', suspendWhenHidden=FALSE)
  
  

} # server

######## CREATE APP ########
shinyApp(ui = ui, server = server)

