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
        
        newTable <- data.frame(symbol = X, basemean = BaseMean, log2FoldChange = log2FC, pvalue = Pvalue, padj = Qvalue)
      })
    
    dataIDKegg = reactive({ 
      data = dataComplet()
      x = data[,1]
      
      ids = bitr(x, fromType="ENSEMBL", toType=c("UNIPROT", "PFAM"), OrgDb="org.Hs.eg.db")
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
      pathview(gene.data = gse, pathway.id = id, species = "hsa")
      
    })
  
    output$PathwayEnrichment <- renderPrint({
      pathway <- IDpathway()
      id <- pathway$id
      print("Téléchargement de la voie métabolique hsa",id)
    })
    
    
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

    output$domain_ID <- renderDataTable({
      data <- dataComplet()
      
      #requete avec le jeu de données
      data_gene_id <- data[1]
      requete_data <- getBM(attributes=c('ensembl_gene_id', 'pfam'), filters = 'ensembl_gene_id', values = data_gene_id, mart = h_sapiens)
      nb_gene_data <- nrow(data_gene_id)
      occurences <- table(requete_data[2])
      mauvaise_ligne <- which(rownames(occurences) == "")
      occurences <- occurences[-mauvaise_ligne] #supprime la ligne avec "" comme nom de domaine
      tableau_data <- as.data.frame(occurences)
      
      #requete avec le genome
      data_genome_id <- getBM(attributes = c('ensembl_gene_id'), mart = h_sapiens)
      nb_gene_genome <- nrow(data_genome_id)
      requete_genome <- getBM(attributes=c('ensembl_gene_id', 'pfam'), filters = 'ensembl_gene_id', values = data_genome_id, mart = h_sapiens)
      occ_genome <- table(requete_genome[2])
      mauvaise_ligne_genome <- which(rownames(occ_genome) == "")
      occ_genome <- occ_genome[-mauvaise_ligne_genome]
      tableau_genome <- as.data.frame(occ_genome)
      tableau_genome
      
      cpt <- 1
      tableau_final <- data.frame(matrix(NA,ncol=5,nrow=nrow(tableau_data)))
      for (x in 1:nrow(tableau_data)){
        ligne_genome <- rownames(subset(tableau_genome, Var1 == as.character(tableau_data[cpt,1])))
        m <- matrix(c(tableau_data[cpt,2],nb_gene_data,tableau_genome[ligne_genome,2],nb_gene_genome), nrow=2, ncol=2)
        # Test du Chi 2 #
        res_chisq <- chisq.test(m)
        pvalue_chisq <- res_chisq$p.value
        # Test de Fisher #
        res_fish <- fisher.test(m)
        pvalue_fisher <- res_fish$p.value

        # Construction du tableau final #
        tableau_final[cpt,1]= tableau_data[cpt,1]
        tableau_final[cpt,2]= pvalue_chisq
        tableau_final[cpt,3]= 0
        tableau_final[cpt,4]= pvalue_fisher
        tableau_final[cpt,5]= 0
        
        cpt = cpt + 1
        
        
      }
      tableau_final
    })
    
  }
)