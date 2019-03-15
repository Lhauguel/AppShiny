library(shiny)
library(DT)
library(ggplot2)
library(plotly)
library("shinythemes")
library(shinycssloaders)

# Define UI for application that draws a histogram
shinyUI(
  fluidPage(
    theme = shinytheme("sandstone"),
    tagList(
      themeSelector(),
      navbarPage(
        "Functional enhancement with Shiny",
        
        ##################################################
        ############ Première page Input Data ############
        ##################################################
        
        tabPanel(
          "Input Data",
          
          # Sidebar with a slider input for number of bins 
          sidebarLayout(
            sidebarPanel( 
              width = 11,
              h1("Shiny application for enrichment analysis"),
              
              fileInput(
                'fichier1', 'Select a file',
                accept = c(
                  'text/csv',
                  "text/comma-separated-values,text/plain",
                  '.csv',
                  '.tsv'
                ), 
                placeholder = "Import you data (csv,tsv)",
                width = "50%"
              ),
              
              h2("Settings"),
              
              ## Choix base de données statistiques et organisme en colonne
              fluidRow(
                column(4, radioButtons("GeneID", label = "origin gene IDs", choices = list("NCBI" = 1, "Ensembl" = 2), selected = 1)),
                column(4, radioButtons("Stat", label = "Statistics", choices = list("GSEA" = 1, "SEA" = 2), selected = 1)),
                column(4, uiOutput(("toCol")))
              ),
              
              h4("Choice p-value and q-value"),
              ## Choix p-value et q-value  
              fluidRow(
                column(4,
                       numericInput("pValueID", label = "p-Value", value = 0.05, min = 0, max = 1, step = 0.01, width = "40%")),
                column(4, 
                       numericInput("qValueID", label = "q-Value", value = 0.05, min = 0, max = 1, step = 0.01, width = "40%")),
                column(4,
                       numericInput("log2FCID", label = "log 2 Fold change", value = 0.05, min = 0, max = 1, step = 0.01, width = "60%"))
              ),
              fluidRow(
                column(6, align="right", offset = 6,
                       actionButton("Start", "Start analysis")
                )
              )
            ),
            mainPanel(
              conditionalPanel(
                condition = "input.Start == false",
                h4("DataFile example: "),
                # Return a list containing the filename and alt text
                img(src = "Fichier_titre.png")),
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
            sidebarPanel( 
              width = 4,
              h1("Whole Data Inspection"),
              uiOutput("sliderQValue"),
              uiOutput("sliderFC")
            ),
            mainPanel(
              withSpinner (
                type  =  getOption ( "spinner.type" , default =  sample(1:8,1)), 
                color =  getOption ( "spinner.color" , default =  "#333333" ), 
                color.background =  getOption ( "spinner.color.background" , default =  "#333333" ),
                tabsetPanel(
                  tabPanel("VulcanoPlot", plotlyOutput("Vulcano")), 
                  tabPanel("MAPlot", plotlyOutput("MAPlot"))
                )
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
            sidebarPanel( 
              width = 4,
              h1("GO Term Enrichment"),
              numericInput("level", label = h4("Go level"), value = 2, max = 8, min = 1, width = "60%"),
              numericInput("category", label = h4("Go category"), value = 25, max = 100, min = 1, width = "60%"),
              h2("Settings"),
              
              ## Choix base de données statistiques et organisme en colonne
              fluidRow(
                column(4, radioButtons("Stat1", label = "Statistics", choices = list("GSEA" = 1, "SEA" = 2))),
                column(6, 
                       conditionalPanel(
                         condition = "input.Stat1 == '2'",
                         numericInput("qValueIDSEA2", label = "q-Value", value = 0.05, min = 0, max = 1, step = 0.01)
                       )
                )
              ),
              fluidRow(
                column(12,align = "center",
                       actionButton("StartGO", "Start")
                )
              )
              
            ),
            
            mainPanel(
              downloadButton('GoTermPlot','Download', width = "30%", align = "right"),
              withSpinner (type  =  getOption ( "spinner.type" , default =  sample(1:8,1)), 
                           color =  getOption ( "spinner.color" , default =  "#333333" ), 
                           color.background =  getOption ( "spinner.color.background" , default =  "#333333" ), 
                           plotOutput("GroupGO")),
              
              dataTableOutput("GOID")
            )
          )
        ),
        
        ###################################################
        ######## Quatrième page Pathway Enrichment ########
        ###################################################
        
        tabPanel(
          "Pathway Enrichment",
          sidebarLayout(
            sidebarPanel( 
              width = 4,
              h1("Pathway Enrichment"),
              h2("Settings"),
              
              ## Choix base de données statistiques et organisme en colonne
              fluidRow(
                uiOutput("toPathway")
              )
            ),
            mainPanel(
              withSpinner (
                type  =  getOption ( "spinner.type" , default =  sample(1:8,1)), 
                color =  getOption ( "spinner.color" , default =  "#333333" ), 
                color.background =  getOption ( "spinner.color.background" , default =  "#333333" ),
                imageOutput("PathwayEnrichment"))
            )
          )
        ),
        
        ###################################################
        ######## Cinquième page Protein Enrichment ########
        ###################################################
        
        tabPanel(
          "Protein Enrichment",
          sidebarLayout(
            sidebarPanel( 
              width = 4,
              h1("Protein Enrichment"),
              h5("Please, choose the settings and launch 'Analysis Protein'"),
              uiOutput("sliderQValue2"),
              
              ## Choix base de données statistiques et organisme en colonne
              fluidRow(
                column(4, radioButtons("Ajust", label = "Ajustment", choices = list("holm" = 1, "hochberg" = 2, "hommel" = 3, "bonferroni" = 4, "BH" = 5, "BY" = 6, "fdr" = 7), selected = 4)),
                column(4, radioButtons("TestStat", label = "Statistics Test", choices = list("Chi2" = 1, "Fischer" = 2), selected = 1))
              ),
              
              uiOutput("sliderPadj"),
              
              fluidRow(
                column(4, radioButtons("Stat2", label = "Statistics", choices = list("GSEA" = 1, "SEA" = 2))),
                column(6, 
                       conditionalPanel(
                         condition = "input.Stat2 == '2'",
                         numericInput("qValueIDSEA", label = "q-Value", value = 0.05, min = 0, max = 1, step = 0.01)
                       )
                )
              ),
              
              fluidRow(
                column(12,align = "center",
                       actionButton("StartProteine", "Analysis Protein")
                )
              ),
              
              tags$head(
                tags$style(HTML('#info{background-color:gray}'))
              ),
              
              fluidRow(
                column(12,align = "center",
                       actionButton("info", "About"),
                       conditionalPanel(
                         condition = "input.info",
                         h5("The adjustment methods include the Bonferroni correction ('bonferroni') in which the p-values are multiplied by the number of comparisons."), 
                         h5("Less conservative corrections are also included by Holm (1979) ('holm'), Hochberg (1988) ('hochberg'), Hommel (1988) ('hommel'), Benjamini & Hochberg (1995) ('BH' or its alias 'fdr'), and Benjamini & Yekutieli (2001) ('BY'), respectively. A pass-through option ('none') is also included. The set of methods are contained in the p.adjust.methods vector for the benefit of methods that need to have the method as an option and pass it on to p.adjust."), 
                         h5("The first four methods are designed to give strong control of the family-wise error rate. There seems no reason to use the unmodified Bonferroni correction because it is dominated by Holm's method, which is also valid under arbitrary assumptions. Hochberg's and Hommel's methods are valid when the hypothesis tests are independent or when they are non-negatively associated (Sarkar, 1998; Sarkar and Chang, 1997)."), 
                         h5("Hommel's method is more powerful than Hochberg's, but the difference is usually small and the Hochberg p-values are faster to compute."),
                         h5("The 'BH' (aka 'fdr') and 'BY' method of Benjamini, Hochberg, and Yekutieli control the false discovery rate, the expected proportion of false discoveries amongst the rejected hypotheses. "),
                         h5("The false discovery rate is a less stringent condition than the family-wise error rate, so these methods are more powerful than the others. Note that you can set n larger than length(p) which means the unobserved p-values are assumed to be greater than all the observed p for 'bonferroni' and 'holm' methods and equal to 1 for the other methods.")
                       )
                )
              )
            ),
            mainPanel(
              withSpinner (
                type  =  getOption ( "spinner.type" , default =  sample(1:8,1)), 
                color =  getOption ( "spinner.color" , default =  "#333333" ), 
                color.background =  getOption ( "spinner.color.background" , default =  "#333333" ),
                dataTableOutput("domain_ID")),
                plotlyOutput("plotdomain_ID")
            )
          )
        )
      )
    )
  )
)