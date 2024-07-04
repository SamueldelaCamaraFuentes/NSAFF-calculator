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
  dashboardHeader(title = "NSAFF Calculator"),
  
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
    multi <- readxl::read_xlsx(input$file$datapath, sheet = "Proteins")
    
    if (input$comptplatform == 1){
      psm_columns <- grep("PSM", colnames(multi), value = TRUE)
      
      if (input$method == 1){
        aa_column <- grep("AAs", colnames(multi), value = TRUE)
        
      } else if (input$method == 2){
        aa_column <- grep("MW", colnames(multi), value = TRUE)
      }
      
    } else if (input$comptplatform == 2){
      psm_columns <- grep("Spectral Count", colnames(multi), value = TRUE)
      aa_column <- grep("Length", colnames(multi), value = TRUE)
    }
  
    # Convert columns to numeric
    multi[c(psm_columns, aa_column)] <- sapply(multi[c(psm_columns, aa_column)], as.numeric)
    
    # Replace NAs in PSM columns with 0
    multi[, psm_columns] <- replace(multi[, psm_columns], is.na(multi[, psm_columns]), 0)
    
    # Initialize a vector to store specified total lengths
    specified_total_lengths <- numeric(length(psm_columns))
    
    # Loop through PSM columns to calculate specified total lengths
    for (i in seq_along(psm_columns)) {
      psm_col <- psm_columns[i]
      
      # Calculate total length excluding rows with PSM value of 0
      specified_total_lengths[i] <- sum(multi[multi[[psm_col]] != 0, aa_column])
    }
    
    
    
    # Loop through PSM columns to calculate NSAF with specified total length
    for (i in seq_along(psm_columns)) {
      psm_col <- psm_columns[i]
      Total_PSMs <- sum(multi[, psm_col])
      
      # Calculate NSAF for the current PSM column with specified total length
      multi[[paste0("NSAF_", gsub("PSM", "", psm_col))]] <- (multi[[psm_col]] / multi[[aa_column]]) / (Total_PSMs / specified_total_lengths[i])
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
