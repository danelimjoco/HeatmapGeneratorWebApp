library(shiny)
library(shinythemes)
library(readxl)


# Can put stuff from other script here
# Remove genes with less than 10 counts across samples (double-check how to select threshold)
 
 #cervical<- read_excel(inFile$datapath, range = "A1:D6")
 #cervical <- cervical[rowSums(cervical)>10,]



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
                           br(), #some space
                           br(), #some space
                           conditionalPanel(condition = "output.metadataUploaded", 
                                            h3("3. Filter out low counts"), 
                                            
                                            # Input: Simple integer interval ----
                                            sliderInput("lowCount", textOutput("lowCounts"),
                                                        min = 0, max = 1000,
                                                        value = 100),
                                            
                                            actionButton("sendLowCount", "Done!", style="color: #	#87CEFA; background-color: #00FF7F; border-color: #2e6da4"),
                                            
                                            br(), #some space 
                                            br(), #some space 
                                            
                                            #h4(textOutput("lowCountMsg")),
                                            h4(textOutput("genesLeft")),
                                            
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
  
  # Low counts filter header
  output$lowCounts <- renderText({
    paste0("Drop genes with less than ", input$lowCount, " read counts across all samples")
  })
  
  # Low counts filter confirmation message
  #output$lowCountMsg <- renderText({
  #  if (input$sendLowCount>0) { 
  #    return(paste0("Genes with less than ", input$lowCount, " reads across samples dropped!"))
  #  }
  #})
  
  # Genes remaining message
  output$genesLeft <- renderText({
    if (input$sendLowCount>0) { 
      file <- input$countsFile
      req(file)
      countData <- read_excel(file$datapath)
      reducedCountData <- countData[rowSums(countData[,-1])>input$lowCount,]
      return(paste0(nrow(countData)-nrow(reducedCountData)," genes dropped! ", nrow(reducedCountData), " genes left for analysis."))
    }
  })

  outputOptions(output, 'fileUploaded', suspendWhenHidden=FALSE)
  outputOptions(output, 'metadataUploaded', suspendWhenHidden=FALSE)
  
} # server

######## CREATE APP ########
shinyApp(ui = ui, server = server)

