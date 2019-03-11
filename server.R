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
library("png")

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
    
    requeteGenome <- reactive ({
      h_sapiens = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
      #requete avec le genome
      data_genome_id <- getBM(attributes = c('ensembl_gene_id'), mart = h_sapiens)
      getBM(attributes=c('ensembl_gene_id', 'interpro'), filters = 'ensembl_gene_id', values = data_genome_id, mart = h_sapiens)
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
    output$sliderQValue <- renderUI({ numericInput("qValueID1", label = h4("q-value"), 
                                                   value = input$qValueID, min = 0, max = 1, step = 0.01, width = "60%")
    })
    output$sliderFC <- renderUI({ numericInput("log2FCID1", label = h4("log 2 Fold change"), min = 0, 
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
                 geom_point(aes(text = paste("Symbol:", id)), size = 0.5) +
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
                 geom_point(aes(text = paste("Symbol:", id)), size = 0.5) +
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
    
    output$sliderPValue <- renderUI({ numericInput("pValueID1", label = h4("p-value"), 
                                                    value = input$pValueID, min = 0, max = 1, step = 0.01, width = "60%")
    })
    
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
      readPNG("hsa00010.pathview.png")
    })
  
    output$PathwayEnrichment <- renderImage({
      IDpathway()
      
    })

    
    
    ##################################################
    ## Cinquième page Protein Enrichment
    ##################################################
    
    output$sliderQValue2 <- renderUI({ numericInput("qValueID2", label = h4("q-value"), 
                                                   value = input$qValueID, min = 0, max = 1, step = 0.01, width = "60%")
    })
    
    ajustement <- reactive({
      if (input$Ajust == 1) ajust = "holm"
      if (input$Ajust == 2) ajust = "hochberg"
      if (input$Ajust == 3) ajust = "hommel"
      if (input$Ajust == 4) ajust = "bonferroni"
      if (input$Ajust == 5) ajust = "BH"
      if (input$Ajust == 6) ajust = "BY"
      else ajust = "fdr"
    })
    
    testStat <- reactive({
      if (input$TestStat == 1) test = "Chi2"
      else test = "Fischer"
    })
    
    output$contents_gene <- renderDataTable({
      dataGene()
    })

    output$domain_ID <- renderDataTable({
      data <- dataComplet()
      ajustement <- ajustement()
      requete_genome <- requeteGenome()
      testStatistique <- testStat()
      h_sapiens = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
      
      #requete avec le jeu de données
      data_gene_id <- data[1]
      requete_data <- getBM(attributes=c('ensembl_gene_id', 'interpro'), filters = 'ensembl_gene_id', values = data_gene_id, mart = h_sapiens)
      nb_gene_data <- nrow(data_gene_id)
      occurences <- table(requete_data[2])
      mauvaise_ligne <- which(rownames(occurences) == "")
      occurences <- occurences[-mauvaise_ligne] #supprime la ligne avec "" comme nom de domaine
      tableau_data <- as.data.frame(occurences)
      
      #requete genome 
      data_genome_id <- getBM(attributes = c('ensembl_gene_id'), mart = h_sapiens)
      nb_gene_genome <- nrow(data_genome_id)
      occ_genome <- table(requete_genome[2])
      mauvaise_ligne_genome <- which(rownames(occ_genome) == "")
      occ_genome <- occ_genome[-mauvaise_ligne_genome]
      tableau_genome <- as.data.frame(occ_genome)
      
      # Construction du tableau final
      cpt <- 1
      tableau_final <- matrix(data="NA", nrow=nrow(tableau_data), ncol=3, dimnames=list(c(), c("Domain ID", "pvalue", "padj")), byrow = TRUE)
      
      liste_test <- list()
      for (x in 1:nrow(tableau_data)){
        ligne_genome <- rownames(subset(tableau_genome, Var1 == as.character(tableau_data[cpt,1])))
        m <- matrix(c(tableau_data[cpt,2],nb_gene_data,tableau_genome[ligne_genome,2],nb_gene_genome), nrow=2, ncol=2)
        
        if (testStatistique == "Chi2") {
          # Test du Chi 2 #
          res_chisq <- chisq.test(m)
          pvalue_test <- res_chisq$p.value
          liste_test <- append(liste_test, pvalue_test)
        }
        
        else {
          # Test de Fisher #
          res_fish <- fisher.test(m)
          pvalue_test <- res_fish$p.value
          liste_test <- append(liste_test, pvalue_test)
        }
      
        # Construction du tableau final #
        tableau_final[cpt,1]= as.character(tableau_data[cpt,1])
        tableau_final[cpt,2]= pvalue_test
        tableau_final[cpt,3]= 0
        
        cpt = cpt + 1
      }
      

      # Calcul de la pvalue adj pour le test du chi2
      cpt = 1
      liste_adj_test = p.adjust(liste_test, method = ajustement)
      for (x in 1:length(liste_adj_test)){
        tableau_final[cpt, 3] = liste_adj_test[cpt]
        cpt = cpt + 1
      }
      tableau_final
    })

  }
)