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
        sidebarPanel(
          fileInput('fichier1', 'Sélectionner un fichier',
                    accept = c(
                      'text/csv',
                      'text/comma-separated-values',
                      'text/tab-separated-values',
                      'text/plain',
                      '.csv',
                      '.tsv'
                    )
          ),
          fileInput('fichier2', 'Sélectionner un fichier',
                    accept = c(
                      'text/csv',
                      'text/comma-separated-values',
                      'text/tab-separated-values',
                      'text/plain',
                      '.csv',
                      '.tsv'
                    )
          ),
          selectInput(
            "Gene",
            "Selection de genes :",
            choices = NULL
          )
        ), 
        # Show a plot of the generated distribution
        mainPanel(
          plotOutput("boxplot")
          #tableOutput("contents1"),
          #tableOutput("contents2")
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
