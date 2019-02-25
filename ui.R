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
shinyUI(
  fluidPage(
    tagList(
      shinythemes::themeSelector(),
      navbarPage(
      "Enrichissement fonctionnelle avec Shiny",
      
      ##################################################
      ############ Première page Input Data ############
      ##################################################
      
      tabPanel(
        "Input Data",
        
        # Sidebar with a slider input for number of bins 
        sidebarLayout(
          sidebarPanel( width = 11,
            h1("Shiny application for enrichment analysis"),
            fileInput('fichier1', 'Sélectionner un fichier',
                      accept = c(
                        'text/csv',
                        "text/comma-separated-values,text/plain",
                        '.csv'#,
                        #'.tsv'
                      ), placeholder = "Import you data (csv,tsv)",
                      width = "50%"),
            
            h2("Settings"),
            
            ## Choix base de données statistiques et organisme en colonne
            fluidRow(
              column(4, radioButtons("GeneID", label = "origine gene IDs", choices = list("Gene NCBI" = 1, "Ensembl" = 2), selected = 1)),
              column(4, radioButtons("Stat", label = "Statistics", choices = list("GSEA" = 1, "SEA" = 2), selected = 1)),
              column(4, uiOutput(("toCol")))
            ),
            
            h4("Choice p-value and q-value"),
            ## Choix p-value et q-value  
            numericInput("pValueID", label = "p-Value", value = 0.05, min = 0, max = 1, step = 0.01, width = "25%"),
            numericInput("qValueID", label = "q-Value", value = 0.05, min = 0, max = 1, step = 0.01, width = "25%"),
            numericInput("log2FCID", label = "log 2 Fold change", value = 0.05, min = 0, max = 1, step = 0.01, width = "25%"),
            
            fluidRow(
              column(6, align="right", offset = 6,
                     actionButton("Start", "Start analysis")
              )
            )
          ),
          mainPanel(
            dataTableOutput("contents")
          )
        )
      ),
      
      ###################################################
      ####### Deuxieme page Whole Data Inspection #######
      ###################################################
      
      tabPanel(
        "Whole Data Inspection", 
        sidebarLayout(
          sidebarPanel( width = 9,
            h1("Whole Data Inspection"),
            uiOutput("sliderQValue"),
            uiOutput("sliderFC"),
            ## Slider mis dans server.R et apparait ici, permet de mettre en dépendance 
            ## avec la valeur initiale mise dans la première page 
            h2("Settings"),
            
            ## Choix base de données statistiques et organisme en colonne
            fluidRow(
              column(4, radioButtons("Stat", label = "Statistics", choices = list("GSEA" = 1, "SEA" = 2), selected = 1))
            )
          ),
          mainPanel(
            tabsetPanel(
              tabPanel("VulcanoPlot", plotlyOutput("Vulcano")), 
              tabPanel("MAPlot", plotlyOutput("MAPlot"))
            )
          )
        )
      ),
  
      ###################################################
      ######## Troisieme page GO Term Enrichment ########
      ###################################################
      
      tabPanel(
        "GO Term Enrichment",
        sidebarLayout(
          sidebarPanel( width = 9,
            h1("GO Term Enrichment"),
            sliderInput("PathOption1", label = h4("option 1"), min = 0, 
                        max = 100, value = 50, width = "60%"),
            numericInput("level", label = h4("Go level"), value = 2),
            h2("Settings"),
            
            ## Choix base de données statistiques et organisme en colonne
            fluidRow(
              column(4, radioButtons("Stat", label = "Statistics", choices = list("GSEA" = 1, "SEA" = 2), selected = 1))
            )
          ),
          mainPanel(
            plotOutput("GroupGO"),
            tableOutput("GOID")
          )
        )
      ),
      
      ###################################################
      ######## Quatrième page Pathway Enrichment ########
      ###################################################
      
      tabPanel(
        "Pathway Enrichment",
      sidebarLayout(
        sidebarPanel( width = 9,
          h1("Pathway Enrichment"),
          sliderInput("PathOption1", label = h4("option 1"), min = 0, 
                      max = 100, value = 50, width = "60%"),
          h2("Settings"),
          
          ## Choix base de données statistiques et organisme en colonne
          fluidRow(
            column(4, radioButtons("Stat", label = "Statistics", choices = list("GSEA" = 1, "SEA" = 2), selected = 1))
          )
          ),
          mainPanel(
          )
        )
      ),
      
      ###################################################
      ######## Cinquième page Protein Enrichment ########
      ###################################################
      
      tabPanel(
        "Protein Enrichment",
        sidebarLayout(
          sidebarPanel( width = 9,
            h1("Protein Enrichment"),
            sliderInput("ProtOption1", label = h4("option 1"), min = 0, 
                        max = 100, value = 50, width = "60%"), h2("Settings"),
            
            ## Choix base de données statistiques et organisme en colonne
            fluidRow(
              column(4, radioButtons("Stat", label = "Statistics", choices = list("GSEA" = 1, "SEA" = 2), selected = 1))
            )
          ),
          mainPanel(
          )
          
        )
      )
    )
  ))
)
