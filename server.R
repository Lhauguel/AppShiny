#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(DT)
library(ggplot2)

##source("http://bioconductor.org/biocLite.R")

## install biomart
##biocLite("biomaRt")
##library(biomaRt)

## exemple biomart
##install.packages(biomaRt)

##data<- read.table("Donnees_test.csv", sep = ",")
##GeneID <- data[1]

##ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
##ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", version=78)
##ensembl
##ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="uswest.ensembl.org")




## install clusterProfiler
##biocLite("clusterProfiler")
##vignette("clusterProfiler", package="clusterProfiler")  ## documentation
##library(clusterProfiler)

## install pathview
##biocLite("pathview")
##browseVignettes("pathview")  ## documentation
##library(pathview)


# Define server logic required to draw a histogram
shinyServer(
    function(session, input, output) {
      dataComplet = reactive({req(input$fichier1)
        data <- read.csv(input$fichier1$datapath, header = TRUE, sep = ",")
        data
      })
      

      
      ##################################################
      ## Première page Input Data
      ##################################################
      
      ## Premier fichier
      output$contents <- renderDataTable({
        dataComplet()
      })
      
      output$valueGeneID <- renderPrint({ input$GeneID })
      output$valueStat <- renderPrint({ input$Stat })
      output$valueOrga <- renderPrint({ input$NameOrga })
      
      output$valuePValueID <- renderText({ input$pValueID }) 
      output$valueQValueID <- renderText({ input$qValueID })
      
      output$valueStart <- renderPrint({ input$Start })
      
      ##################################################
      ## Deuxieme page Whole Data Inspection
      ##################################################

      ## Test de volcano plot mettre les limites de Foldchange et de -log de pvalue ajustée 
      
      output$Vulcano = renderPlot({
        data = dataComplet()
        data$Expression = as.factor(ifelse(
          (data$pvalue < input$qValueID1) & (data$log2FoldChange > 0),
          "Significativly Over-Expressed",
          ifelse(
            data$pvalue > input$qValueID1,
            "Not Significant",
            "Significativly Under-Expressed"
          )
        ))
        ggplot(data,
               aes(
                 x = log2FoldChange,
                 y = -log10(padj),
                 col = Expression
               )) +
          geom_point(alpha = 0.5) +
          xlim(c(-5, 5)) + ylim(c(0, 15)) +
          xlab("log2 fold change") +
          ylab("p-value")
      })
      
      ## Fait apparaitre le slider dans l'interface
      output$sliderQValue <- renderUI({ sliderInput("qValueID1", label = h4("q-value"), min = 0, 
                                                    max = 1, value = input$qValueID, width = "60%") 
                                      })
      output$valueDataOption2 <- renderPrint({ input$DataOption2 })
      output$valueDataOption3 <- renderPrint({ input$DataOption3 })
      output$valueDataOption4 <- renderPrint({ input$DataOption4 })
      
      ##################################################
      ## Troisieme page GO Term Enrichment
      ##################################################
      
      output$valueGOTermOption1 <- renderPrint({ input$GOTermOption1 })
      output$valueGOTermOption2 <- renderPrint({ input$GOTermOption2 })
      output$valueGOTermOption3 <- renderPrint({ input$GOTermOption3 })
      output$valueGOTermOption4 <- renderPrint({ input$GOTermOption4 })
      
      ##################################################
      ## Quatrième page Pathway Enrichment
      ##################################################
      
      output$valuePathOption1 <- renderPrint({ input$PathOption1 })
      output$valuePathOption2 <- renderPrint({ input$PathOption2 })
      output$valuePathOption3 <- renderPrint({ input$PathOption3 })
      output$valuePathOption4 <- renderPrint({ input$PathOption4 })
      
      ##################################################
      ## Cinquième page Protein Enrichment
      ##################################################

      
      output$valueProtOption1 <- renderPrint({ input$ProtOption1 })
      output$valueProtOption2 <- renderPrint({ input$ProtOption2 })
      output$valueProtOption3 <- renderPrint({ input$ProtOption3 })
      output$valueProtOption4 <- renderPrint({ input$ProtOption4 })
      
  }
)
  

