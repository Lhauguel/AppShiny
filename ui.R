#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(DT)
library(ggplot2)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  navbarPage(
    "Enrichissement fonctionnelle avec Shiny",
    
    tabPanel(
      "Input Data",
      
      # Sidebar with a slider input for number of bins 
      sidebarLayout(
        sidebarPanel( width = 9,
          fileInput('fichier1', 'Sélectionner un fichier',
                    accept = c(
                      'text/csv',
                      'text/comma-separated-values',
                      'text/tab-separated-values',
                      'text/plain',
                      '.csv',
                      '.tsv'
                    ), placeholder = "Import you data (csv,tsv)",
                    width = "60%"),
          
          h2("Settings"),
          
          ## Choix base de données statistiques et organisme en colonne
          fluidRow(
            column(4, radioButtons("GeneID", label = "origine gene IDs", choices = list("Gene NCBI" = 1, "Ensembl" = 2), selected = 1)),
            column(4, radioButtons("Stat", label = "Statistics", choices = list("GSEA" = 1, "SEA" = 2), selected = 1)),
            column(4, radioButtons("NameOrga", label = "Organism Name", choices = list("biomaRt" = 1, "autre" = 2), selected = 1))
          ),
          
          h4("Choice p-value and q-value"),
          ## Choix p-value et q-value  
          numericInput("pValueID", label = "p-Value", value = 0.05, min = 0, max = 1, step = 0.01, width = "25%"),
          numericInput("qValueID", label = "q-Value", value = 0.05, min = 0, max = 1, step = 0.01, width = "25%")
    ),
       
        
        
        mainPanel(
        )
      )
      
    ),
     
    tabPanel(
      "Whole Data Inspection"
    ),
    tabPanel(
      "GO Term Enrichment"
    ),
    tabPanel(
      "Pathway Enrichment"
    ),
    tabPanel(
      "Protein enrichment"
    )
  )
)
)

