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
        sidebarPanel( width = 6,
          fileInput('fichier1', 'Sélectionner un fichier',
                    accept = c(
                      'text/csv',
                      'text/comma-separated-values',
                      'text/tab-separated-values',
                      'text/plain',
                      '.csv',
                      '.tsv'
                    ), placeholder = "Import you data (csv,tsv)"
          ),
          h2("Settings")
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
