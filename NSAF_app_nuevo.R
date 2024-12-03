#########################################################
#
# NSAFF calculation app that uses Proteome discoverer 
# output files in order to perform the NSAFF calculation
#
#########################################################


if(!require(shiny)){install.packages("shiny")}
library(shiny)
if(!require(shinyWidgets)){install.packages("shinyWidgets")}
library(shinyWidgets)
if(!require(shinydashboard)){install.packages("shinydashboard")}
library(shinydashboard)
if(!require(shinyjs)){install.packages("shinyjs")}
library(shinyjs)
if(!require(DT)){install.packages("DT")}
library(DT)
if(!require(BiocManager)){install.packages("BiocManager")}
library(BiocManager)
if(!require(svglite)){install.packages("svglite")}
library(svglite)
options(repos = BiocManager::repositories())
options(shiny.maxRequestSize = 30*1024^2)


# Define UI ----
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "NSAF Calculator"),
  
  dashboardSidebar(
    fileInput(inputId = "file", 
              label = h3("File input"),
              multiple = FALSE, 
              accept = c("text/csv", ".csv")),
    selectInput(inputId = "comptplatform", 
                label = "Platform",
                choices = list("Proteome Discoverer" = 1, 
                               "MSFragger" = 2),
                selected = 1),
    radioButtons(inputId = "method",
                 label = "Choice",
                 choices = c("Aminoacid´s length" = 1,
                             "Molecular weight" = 2),
                 selected = 1),
    textInput(inputId = "filenamedownload",
              label = "Result filename",
              value = "data"),
    downloadButton(outputId = "downloaddata", icon("download"),
                   label = "Download",
                   style="display: block; margin: 0 auto; width: 200px; color:black;")
  ),
  
  dashboardBody(
    fluidPage(DTOutput('tabledata'))
  )
)

# Define server logic ----
server <- function(input, output) {
  
  dataset <- reactive({
    validate(need(!is.null(input$file$datapath),
                  "Please select a Proteome Discoverer or MSFragger output file"))
    
    if (input$comptplatform == 1){
      multi <- readxl::read_xlsx(input$file$datapath)
      psm_columns <- grep("PSM", colnames(multi), value = TRUE)
      
      if (input$method == 1){
        aa_column <- grep("AAs", colnames(multi), value = TRUE)
        
      } else if (input$method == 2){
        aa_column <- grep("MW", colnames(multi), value = TRUE)
      }
      
    } else if (input$comptplatform == 2){
      multi <- read.delim(input$file$datapath, sep = "\t", stringsAsFactors = FALSE, colClasses = "character") 
      psm_columns <- grep("Spectral.Count", colnames(multi), value = TRUE)
      aa_column <- grep("Length", colnames(multi), value = TRUE)
    }
  
    # Convert columns to numeric
    multi[c(psm_columns, aa_column)] <- sapply(multi[c(psm_columns, aa_column)], as.numeric)
    
    # Replace NAs in PSM columns with 0
    multi[, psm_columns] <- replace(multi[, psm_columns], is.na(multi[, psm_columns]), 0)
    
    if (input$comptplatform == 1){
      
      
      multi[paste0(gsub("PSM", "Cociente", psm_columns))] <- as.data.frame(lapply(multi[psm_columns], function(col) col / multi[[aa_column]]))
      
      cociente_columns <- paste0(gsub("PSM", "Cociente", psm_columns))
      cociente_sum <- colSums(multi[cociente_columns])
      
      # Divide the new columns by `cociente_sum`
      multi[paste0("NSAF_", gsub("PSM", "", psm_columns))] <- as.data.frame(
        lapply(seq_along(cociente_columns), function(i) multi[[cociente_columns[i]]] / cociente_sum[i])
      )
      
      
    } else if (input$comptplatform == 2){
      
      multi[paste0(gsub("Spectral.Count", "Cociente", psm_columns))] <- as.data.frame(lapply(multi[psm_columns], function(col) col / multi[[aa_column]]))
      
      cociente_columns <- paste0(gsub("Spectral.Count", "Cociente", psm_columns))
      cociente_sum <- colSums(multi[cociente_columns])
      
      # Divide the new columns by `cociente_sum`
      multi[paste0("NSAF_", gsub("Spectral.Count", "", psm_columns))] <- as.data.frame(
        lapply(seq_along(cociente_columns), function(i) multi[[cociente_columns[i]]] / cociente_sum[i])
      )
    }
    
    multi
  })
  
  output$tabledata <- DT::renderDataTable({
    
    
    datatable(dataset(), options = list(pageLength = 5,
                                    lengthMenu = c(5, 10, 15, 20), 
                                    scrollX = T,
                                    autoWidth = TRUE ))
    
    
  })
  
  output$downloaddata <- downloadHandler(
    
    filename = function() {
      paste(input$filenamedownload, Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(dataset(), file)
    }
  )
  
  
}

# Run the app ----
shinyApp(ui = ui, server = server)
