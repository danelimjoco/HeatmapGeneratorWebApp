library(shiny)
library(shinythemes)
library(readxl)


######## UI ########

ui <- fluidPage(titlePanel(h1("RNA-Seq Shiny App!",align = "center")), theme = shinytheme("cerulean"),
                 fluidRow(column(4, align="center", offset=4,
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
                                          
                           ), conditionalPanel(condition = "output.lowCountsSelected", 
                                            h3("4. How many genes do you want to analyze?"), 
                                                          
                                            sliderInput("numGene", textOutput("numGenes"),
                                                        min = 0, max = 100,
                                                        value = 50),
                                            
                                            actionButton("sendGeneCount", "Done!", style="color: #	#87CEFA; background-color: #00FF7F; border-color: #2e6da4"),
                                            plotOutput("heatmap"),
                                            align="center"
                           ),
                
                           
                  
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

  # Only show ones files have been uploaded
  output$fileUploaded <- reactive({
    return(!is.null(input$countsFile))
  })
  
  output$metadataUploaded <- reactive({
    return(!is.null(input$metadataFile))
  })
  
  output$lowCountsSelected <- reactive({
    return(input$sendLowCount>0)
  })
  
  outputOptions(output, 'fileUploaded', suspendWhenHidden=FALSE)
  outputOptions(output, 'metadataUploaded', suspendWhenHidden=FALSE)
  outputOptions(output, 'lowCountsSelected', suspendWhenHidden=FALSE)
  
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
  
  # Gene counts filter header
  output$numGenes <- renderText({
    paste0("Select the top ", input$numGene, " genes to analyze")
  })
  
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
  
  # Carry out DESEq analysis, get the top X genes, and make the heatmap upon submit
  output$heatmap <- renderPlot({
    if (input$sendGeneCount>0) { 
      
      # Count data
      file <- input$countsFile
      req(file)
      countData <- read_excel(file$datapath)
      countData <- countData[rowSums(countData[,-1])>input$lowCount,]
      
      # Metadata
      inFile <- input$metadataFile
      class <- read_excel(inFile$datapath)
      
      #DESeq
      dds <- DESeqDataSetFromMatrix(countData = countData[,-1], colData = class,
                                    design = formula(~condition)) 
      dds <- DESeq(dds)
      res <- results(dds)
      resDF <- as.data.frame(res)
      
      ######## FILTER GENES BY FDR-ADJUSTED P-VALUE ########
      
      # Eliminate genes with NA for FDR-adjusted p-value (read counts too low)
      resDF <-resDF %>% 
        rownames_to_column('gene') %>% 
        filter(!is.na(padj)) %>% 
        column_to_rownames('gene')
      
      # Filter for significant genes by FDR-adjusted p-value threshold
      resSig = resDF[resDF$padj < 0.1,] 
      
      # Store names of top 50 significant DEGs
      sigGenes <- rownames(resSig[1:input$numGene,])
      
      ######## RLOG TRANSFORM (NECESSARY FOR CLUSTERING) ########
      
      # Transform count data using rlog (could also use vst)
      # For high gene counts, transformation is basically log2. 
      # For low gene counts, values are shrunken towards genes' averages across samples
      ddsRLOG <- rlog(dds)
      
      # Convert DESeqTransform object to data frame
      ddsRLOG <- ddsRLOG %>% assay() %>% as.data.frame()
      
      # Only continue clustering with significant DEGs
      ddsRLOG <- ddsRLOG[rownames(ddsRLOG) %in% sigGenes,]
      
      ######## Z-SCORE SCALING ########
      
      # Subtract mean of each row from all values
      ddsRLOG <- ddsRLOG - rowMeans(ddsRLOG)
      
      # Divide all values by std dev of each row
      ddsRLOG$row_std = apply(ddsRLOG[,-1], 1, sd)
      ddsRLOG <- sweep(ddsRLOG, MARGIN=1, ddsRLOG$row_std, FUN="/")
      
      # Don't need row_std column anymore
      ddsRLOG$row_std <- NULL
      
      # Distance calculation and cluster
      distanceGene <- dist(ddsRLOG, method="manhattan") 
      clusterGene <- hclust(distanceGene, method="average")
      ddsRLOG$Gene <- rownames(ddsRLOG)
      ddsRLOG$Gene <- factor(ddsRLOG$Gene, levels=clusterGene$labels[clusterGene$order])
      
      
      ######## DRAW HEATMAP ########
      
      # Put it back in the melted format for ggplot
      ddsRLOG <- reshape2::melt(ddsRLOG, id.vars=c("Gene")) # converts to long version
      
      # Rename value column to "Z-score"
      ddsRLOG <- ddsRLOG %>% 
        rename(Z_score = value)
      
      # Draw heatmap
      heatmap <- ggplot(ddsRLOG, aes(x=variable, y=Gene, fill=Z_score)) +
        labs(x ="Tissue Sample", y = "miRNA gene") + 
        geom_raster() + 
        scale_fill_gradient2(midpoint = 0, low = "blue", mid = "white", high = "red") + 
        theme(aspect.ratio = 0.9, axis.text.x=element_text(size=10, angle=65, hjust=1), axis.text.y=element_text(size=10))
      
      return(heatmap)

    }
  })

  
} # server

######## CREATE APP ########
shinyApp(ui = ui, server = server)

