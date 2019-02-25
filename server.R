#################################################################################
###### This is the server logic of a Shiny web application. You can run the #####
#################### application by clicking 'Run App' above. ###################
#################################################################################
########## Find out more about building applications with Shiny here: ###########
########################### http://shiny.rstudio.com/ ###########################
#################################################################################

## A nstallé indépendemment du script  
## source("http://bioconductor.org/biocLite.R")
## biocLite("clusterProfiler")

library(shiny)
library(DT)
library(ggplot2)
library(clusterProfiler)
library(biomaRt)
library(org.Hs.eg.db)
library(biomaRt)

##source("http://bioconductor.org/biocLite.R")

## install biomart
## biocLite("biomaRt")

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
      
      dataGene = reactive({req(input$fichier1)
        data <- read.csv(input$fichier1$datapath, header = TRUE, sep = ",")
        data[1]

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
      
      ensembl = useEnsembl(biomart="ensembl")
      list_ensembl = listDatasets(ensembl)[2]
      output$toCol <- renderUI({
        selectInput("BiomaRtOrgo", "Organism Name", list_ensembl, selected = "Human genes (GRCh38.p12)")
      })
      
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
          scale_colour_discrete(drop=FALSE) +
          xlim(c(-1.5, 1.5)) + ylim(c(0, 3.2)) +
          xlab("log2 fold change") +
          ylab("-log10(p-value)")
      })
      
      ## Fait apparaitre le slider dans l'interface
      output$sliderQValue <- renderUI({ sliderInput("qValueID1", label = h4("q-value"), min = 0, 
                                                    max = 1, value = input$qValueID, width = "60%") 
                                      })
      
      #ensembl = useEnsembl(biomart="ensembl")
      #input$data_gene = dataComplet[1]
      #list_ensembl = listDatasets(ensembl)[2]
      #output$toCol <- renderUI({
      #  selectInput("BiomaRtOrgo", "Organism Name", list_ensembl, selected = "Human genes (GRCh38.p12)")
      #})
      output$data_gene <- renderPrint({ input$data_gene })
      output$valueDataOption2 <- renderPrint({ input$DataOption2 })
      output$valueDataOption3 <- renderPrint({ input$DataOption3 })
      output$valueDataOption4 <- renderPrint({ input$DataOption4 })
      
      ##################################################
      ## Troisieme page GO Term Enrichment
      ##################################################
      
      #########################
      ## Partie ClusterProfiler 
      #########################
      
      output$GroupGO = renderPlot({
        data = dataComplet()
        x = data[,1]
        
        ids = bitr(x, fromType="ENSEMBL", toType=c("UNIPROT", "PFAM"), OrgDb="org.Hs.eg.db")
        hg = ids[,2]
        eg2np = bitr_kegg(hg, fromType='uniprot', toType='kegg', organism='hsa')
        gse = eg2np[,2]
        
        ggo = groupGO(gene     = gse,
                      OrgDb    = org.Hs.eg.db,
                      level    = 3,
                      readable = TRUE)
        barplot(ggo, drop=TRUE, showCategory=12)
      })
      
      output$GOID = renderTable({
        ggo$ID
      })

      
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
      
      output$contents_gene <- renderDataTable({
        dataGene()
      })
      ensembl = useEnsembl(biomart="ensembl")
      h_sapiens = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
      output$pfam <- renderDataTable({
        getBM(attributes=c('ensembl_gene_id', 'pfam'), filters = 'ensembl_gene_id', values = dataGene(), mart = h_sapiens)
      })
      
      
      output$nb_gene_total = renderPrint({ nrow(dataGene()) })
      #nb_gene_total = nrow(dataGene())
      output$occurences = renderPrint({ 
        data_gene_pfam <- getBM(attributes=c('ensembl_gene_id', 'pfam'), filters = 'ensembl_gene_id', values = dataGene(), mart = h_sapiens )
        table(unlist(data_gene_pfam[2])) })
      
  }
)
