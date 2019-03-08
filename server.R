## A installé indépendemment du script  
## source("http://bioconductor.org/biocLite.R")
## biocLite("clusterProfiler")

library(shiny)
library(DT)
library(ggplot2)
library(clusterProfiler)
library(biomaRt)
library(org.Hs.eg.db)
library(plotly)
library(pathview)

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
    dataComplet = eventReactive(input$Start,
      {req(input$fichier1)
        data <- read.csv(input$fichier1$datapath, header = TRUE, sep = ",")
        X <-data[,"geneid"]
        BaseMean <- data[,"baseMean"]
        log2FC <- data[,"log2FoldChange"]
        Pvalue <- data[,"pvalue"]
        Qvalue <- data[,"padj"]
        
        newTable <- data.frame(id = X, basemean = BaseMean, log2FoldChange = log2FC, pvalue = Pvalue, padj = Qvalue)
      })
    
    originGeneID <- reactive({
      if (input$GeneID == 1) originID = "ncbiID"
      else originID = "ensemblID"
    })
    
    dataIDKegg = reactive({ 
      data = dataComplet()
      id <- originGeneID()
      
      if (id == "ncbiID") ID = "ENTREZID"
      else ID = "ENSEMBL"
      x = data[,1]
      
      ids = bitr(x, fromType=ID, toType=c("UNIPROT"), OrgDb="org.Hs.eg.db")
      hg = ids[,2]
      eg2np = bitr_kegg(hg, fromType='uniprot', toType='kegg', organism='hsa')
      gse = eg2np[,2]
      gse
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
    output$valuelog2FCIDD <- renderText({ input$log2FCID })
    
    output$valueStart <- renderPrint({ input$Start })
    
    ensembl = useEnsembl(biomart="ensembl")
    list_ensembl = listDatasets(ensembl)[2]
    output$toCol <- renderUI({
      selectInput("BiomaRtOrgo", "Organism Name", list_ensembl, selected = "Human genes (GRCh38.p12)")
    })
    
    ##################################################
    ## Deuxieme page Whole Data Inspection
    ##################################################
    
    ## Fait apparaitre le slider dans l'interface
    output$sliderQValue <- renderUI({ sliderInput("qValueID1", label = h4("q-value"), min = 0, 
                                                  max = 1, value = input$qValueID, width = "60%") 
    })
    output$sliderFC <- renderUI({ sliderInput("log2FCID1", label = h4("log 2 Fold change"), min = 0, 
                                              max = 2, step = 0.01, value = input$log2FCID, width = "60%") 
    })
    
    ## Test de volcano plot mettre les limites de Foldchange et de -log de pvalue ajustée 
    
    
    output$Vulcano = renderPlotly({
      data = dataComplet()
      log2FC = as.numeric(input$log2FCID1)
      data$Expression = as.factor(ifelse(
        (data$pvalue > input$qValueID1) | (-log2FC < data$log2FoldChange) & (data$log2FoldChange < log2FC),
        "Not Significant",
        ifelse(
          (data$pvalue < input$qValueID1) & (data$log2FoldChange < 0),
          "Under-Expressed",
          "Over-Expressed"
        )
      ))
      ggplotly(ggplot(data,
                      aes(
                        x = log2FoldChange,
                        y = -log10(padj),
                        col = Expression
                      )) +
                 ##  geom_point(alpha = 0.5) + 
                 geom_point(aes(text = paste("Symbol:", symbol)), size = 0.5) +
                 scale_colour_discrete(drop=FALSE) +
                 xlim(c(min(log2(data$log2FoldChange)), max(log2(data$log2FoldChange)))) + ylim(c(min(-log10(data$padj)), max(-log10(data$padj)))) +
                 xlab("log2 fold change") +
                 ylab("-log10(p-value)"))
    })
    
    
    output$MAPlot = renderPlotly({
      data = dataComplet()
      log2FC = as.numeric(input$log2FCID1)
      data$Expression = as.factor(ifelse(
        (data$pvalue > input$qValueID1) | (-log2FC < data$log2FoldChange) & (data$log2FoldChange < log2FC),
        "Not Significant",
        ifelse(
          (data$pvalue < input$qValueID1) & (data$log2FoldChange < 0),
          "Under-Expressed",
          "Over-Expressed"
        )
      ))
      ggplotly(ggplot(data,
                      aes(
                        x = log2(basemean),
                        y = log2FoldChange,
                        col = Expression
                      )) +
                 geom_point(aes(text = paste("Symbol:", symbol)), size = 0.5) +
                 scale_colour_discrete(drop=FALSE) +
                 xlim(c(min(log2(data$basemean)), max(log2(data$basemean)))) + ylim(c(min(data$log2FoldChange), max(data$log2FoldChange))) +
                 xlab("log2 basemean") +
                 ylab("log2 fold change"))
    })
    
    #ensembl = useEnsembl(biomart="ensembl")
    #list_ensembl = listDatasets(ensembl)[2]
    #output$toCol <- renderUI({
    #  selectInput("BiomaRtOrgo", "Organism Name", list_ensembl, selected = "Human genes (GRCh38.p12)")
    #})
    output$valueDataOption2 <- renderPrint({ input$DataOption2 })
    output$valueDataOption3 <- renderPrint({ input$DataOption3 })
    output$valueDataOption4 <- renderPrint({ input$DataOption4 })
    
    ##################################################
    ## Troisieme page GO Term Enrichment
    ##################################################
    
    #########################
    ## Partie ClusterProfiler 
    #########################
    
    
    
    tableGgo <- reactive({
      gse <- dataIDKegg() 
      
      
      ggo = groupGO(gene     = gse,
                    OrgDb    = org.Hs.eg.db,
                    level    = input$level,
                    readable = TRUE)
    })
    
    output$GroupGO = renderPlot({
      ggo <- tableGgo()
      barplot(ggo, drop=TRUE, showCategory = 25)
    })
    
    output$GOID = renderTable({
      ggo <- tableGgo()
      ggo[,1:4]
    })
    
    
    output$valueGOTermOption1 <- renderPrint({ input$GOTermOption1 })
    output$valueGOTermOption2 <- renderPrint({ input$GOTermOption2 })
    output$valueGOTermOption3 <- renderPrint({ input$GOTermOption3 })
    output$valueGOTermOption4 <- renderPrint({ input$GOTermOption4 })
    
    ##################################################
    ## Quatrième page Pathway Enrichment
    ##################################################
    
    
    ## Charge les pathways de l'homme
    data(paths.hsa)
    pathway <- as.list(paths.hsa)
    
    ## affiche une liste des pathways
    pathwayName <- unlist(pathway, use.names=FALSE)
    
    output$toPathway <- renderUI({
      selectInput("pathwayID", "Pathway Name", pathwayName, width = "60%")
    })
    
    ## Met les ID avec le nom des pathways 
    pathwayTable <- aggregate(.~ind,stack(pathway),paste,collapse=' ')
    
    ## Récupère uniquement l'ID sans "hsa"
    IDpathway <- reactive({
      gse <- dataIDKegg() 
      
      idpathway <- as.character(subset(pathwayTable, values==input$pathwayID)$ind)

      id <- substring(idpathway, 4)
      pathwaytOut <- pathview(gene.data = gse, pathway.id = id, species = "hsa")
    })
  
    output$PathwayEnrichment <- renderPrint({
      IDpathway()
      print("Téléchargement de la voie métabolique !")
    })

    
    
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
    
    output$occurences = renderPrint({ 
      #calcul du nombre de gene dans le fichier#
      nb_gene_total <- nrow(dataGene())
      #recherche des ids PFAM en fonction des genes id#
      data_gene_pfam <- getBM(attributes=c('ensembl_gene_id', 'pfam'), filters = 'ensembl_gene_id', values = dataGene(), mart = h_sapiens )
      #calcul des occurences#
      occ_data <- table(unlist(data_gene_pfam[2]))
      #cherche la ligne qui est égale à "" #
      mauvaise_ligne <- which(rownames(occ_data) == "")
      #supprime la ligne qui nous sert pas#
      tableau_pfam <- occ_data[-mauvaise_ligne]
      #calcul de la frequence#
      freq_pfamVSdata <- tableau_pfam / nb_gene_total      #recuperation de la liste de tous les genes du genome#
      liste_gene_genome <- getBM(attributes = c('ensembl_gene_id'), mart = h_sapiens)
      #calcul du nombre de gene dans le genome#
      nb_gene_genome <- nrow(liste_gene_genome)
      #recherche des ids PFAM en fonction des genes id#
      data_genome_pfam <- getBM(attributes=c('ensembl_gene_id', 'pfam'), filters = 'ensembl_gene_id', values = liste_gene_genome, mart = h_sapiens)
      #calcul des occurences#
      occ_genome <- table(unlist(data_genome_pfam[2]))
      occ_genome
      #cherche la ligne degueu#
      mauvaise_ligne_genome <- which(rownames(occ_genome) == "")
      #supprime la ligne#
      tableau_pfam_genome <- occ_genome[-mauvaise_ligne_genome]
      #calcul de la frequence#
      freq_pfamVSgenome = tableau_pfam_genome / nb_gene_genome
      freq_pfamVSgenome
    })
    
  }
)