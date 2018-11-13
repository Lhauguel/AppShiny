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
        sidebarPanel( width = 11,
          h1("Shiny application for enrichment analysis"),
          fileInput('fichier1', 'Sélectionner un fichier',
                    accept = c(
                      'text/csv',
                      'text/comma-separated-values',
                      'text/tab-separated-values',
                      'text/plain',
                      '.csv',
                      '.tsv'
                    ), placeholder = "Import you data (csv,tsv)",
                    width = "50%"),
          
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
          numericInput("qValueID", label = "q-Value", value = 0.05, min = 0, max = 1, step = 0.01, width = "25%"),
          
          fluidRow(
            column(6, align="right", offset = 6,
                   actionButton("Start", "Start analysis")
            )
          )
    ),
       
        
        
        mainPanel(
        )
      )
      
    ),
     
    tabPanel(
      "Whole Data Inspection",
      sidebarLayout(
        sidebarPanel( width = 9,
          h1("Whole Data Inspection"),
          sliderInput("option1", label = h4("option 1"), min = 0, 
                      max = 100, value = 50, width = "60%"),
          sliderInput("option2", label = h4("option 2"), min = 0, 
                      max = 100, value = 50, width = "60%"),
          sliderInput("option3", label = h4("option 3"), min = 0, 
                      max = 100, value = 50, width = "60%"),
          sliderInput("option4", label = h4("option 4"), min = 0, 
                      max = 100, value = 50, width = "60%")
        ),
        
        mainPanel(
        )
      )
    ),
    tabPanel(
      "GO Term Enrichment",
      sidebarLayout(
        sidebarPanel( width = 9,
          h1("GO Term Enrichment"),
          sliderInput("option1", label = h4("option 1"), min = 0, 
                      max = 100, value = 50, width = "60%"),
          sliderInput("option2", label = h4("option 2"), min = 0, 
                      max = 100, value = 50, width = "60%"),
          sliderInput("option3", label = h4("option 3"), min = 0, 
                      max = 100, value = 50, width = "60%"),
          sliderInput("option4", label = h4("option 4"), min = 0, 
                      max = 100, value = 50, width = "60%")
        ),
        
        mainPanel(
        )
      )
    ),
    tabPanel(
      "Pathway Enrichment",
      sidebarLayout(
        sidebarPanel( width = 9,
          h1("Pathway Enrichment"),
          sliderInput("option1", label = h4("option 1"), min = 0, 
                      max = 100, value = 50, width = "60%"),
          sliderInput("option2", label = h4("option 2"), min = 0, 
                      max = 100, value = 50, width = "60%"),
          sliderInput("option3", label = h4("option 3"), min = 0, 
                      max = 100, value = 50, width = "60%"),
          sliderInput("option4", label = h4("option 4"), min = 0, 
                      max = 100, value = 50, width = "60%")
        ),
        
        mainPanel(
        )
      )
    ),
    tabPanel(
      "Protein Enrichment",
    sidebarLayout(
      sidebarPanel( width = 9,
        h1("Protein Enrichment"),
        sliderInput("option1", label = h4("option 1"), min = 0, 
                    max = 100, value = 50, width = "60%"),
        sliderInput("option2", label = h4("option 2"), min = 0, 
                    max = 100, value = 50, width = "60%"),
        sliderInput("option3", label = h4("option 3"), min = 0, 
                    max = 100, value = 50, width = "60%"),
        sliderInput("option4", label = h4("option 4"), min = 0, 
                    max = 100, value = 50, width = "60%")
      ),
      
      mainPanel(
      )
    )
    )
  )
)
)

