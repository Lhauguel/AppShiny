## A installé indépendemment du script  
# source("http://bioconductor.org/biocLite.R")
# biocLite("clusterProfiler")
# biocLite("biomaRt")
# biocLite("pathview")
# biocLite("org.Hs.eg.db")

library(shiny)
library(DT)
library(ggplot2)
library(clusterProfiler)
library(biomaRt)
library(org.Hs.eg.db)
library(plotly)
library(pathview)
library("png")
library(stringr)

# Define server logic required to draw a histogram
shinyServer(
  function(session, input, output) {
    dataComplet = eventReactive(input$Start,
                                {req(input$fichier1)
                                  data <- read.csv(input$fichier1$datapath, header = TRUE, sep = ",")
                                  if(ncol(data) < 3){
                                    data <- read.csv(input$fichier1$datapath, header = TRUE, sep = "\t")
                                  }
                                  Name <-data[,"GeneName"]
                                  X <-data[,"geneid"]
                                  BaseMean <- data[,"baseMean"]
                                  log2FC <- data[,"log2FoldChange"]
                                  Pvalue <- data[,"pvalue"]
                                  Qvalue <- data[,"padj"]
                                  
                                  newTable <- data.frame(id = X, basemean = BaseMean, log2FoldChange = log2FC, pvalue = Pvalue, padj = Qvalue, GeneName = Name)
                                })
    
    requeteGenome <- reactive ({
      h_sapiens = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
      #requete avec le genome
      id <- originGeneID()
      if (id == "ncbiID") ID = "entrezgene"
      else ID = "ensembl_gene_id"
      data_genome_id <- getBM(attributes = c(ID), mart = h_sapiens)
      getBM(attributes=c(ID, 'interpro'), filters = ID, values = data_genome_id, mart = h_sapiens)
    })
    
    originGeneID <- reactive({
      if (input$GeneID == 1) originID = "ncbiID"
      else originID = "ensemblID"
    })
    
    SEA_stat1 <- reactive({
      if (input$Stat1 == 1) stat = "GSEA"
      else stat = "SEA"
    })
    
    dataIDKegg = reactive({ 
      data <- dataComplet()
      stat1 <-SEA_stat1()
      if (stat1 == "SEA"){
        data <- data[which(data[,5] >= input$qValueIDSEA2),]
      }
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
    
    output$contents <- renderDataTable({
      dataComplet()
    })
    
    output$valueStart <- renderPrint({ 
      input$Start 
    })
    
    ensembl = useEnsembl(biomart="ensembl")
    list_ensembl = listDatasets(ensembl)[2]
    
    output$toCol <- renderUI({
      selectInput(
        "BiomaRtOrgo", 
        "Organism Name", 
        list_ensembl, 
        selected = "Human genes (GRCh38.p12)")
    })
    
    ##################################################
    ## Deuxieme page Whole Data Inspection
    ##################################################
    
    output$sliderQValue <- renderUI({ 
      numericInput(
        "qValueID1", 
        label = h4("q-value"), 
        value = input$qValueID, 
        min = 0, 
        max = 1, 
        step = 0.01, 
        width = "60%")
    })
    output$sliderFC <- renderUI({ 
      numericInput(
        "log2FCID1", 
        label = h4("log 2 Fold change"), 
        min = 0, 
        max = 2, 
        step = 0.01, 
        value = input$log2FCID, 
        width = "60%") 
    })
    
    
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
                 geom_point(aes(text = paste("Symbol:", id)), size = 0.5) +
                 scale_colour_discrete(drop=FALSE) +
                 xlim(c(min(log2(data$log2FoldChange)), max(log2(data$log2FoldChange)))) + 
                 ylim(c(min(-log10(data$padj)), max(-log10(data$padj)))) +
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
                 xlim(c(min(log2(data$basemean)), max(log2(data$basemean)))) + 
                 ylim(c(min(data$log2FoldChange), max(data$log2FoldChange))) +
                 xlab("log2 basemean expression") +
                 ylab("log2 fold change"))
    })
    
    ##################################################
    ## Troisieme page GO Term Enrichment
    ##################################################
    
    #########################
    ## Partie ClusterProfiler 
    #########################
    
    output$sliderPValue <- renderUI({ 
      numericInput("pValueID1", label = h4("p-value"), 
                   value = input$pValueID, min = 0, max = 1, step = 0.01, width = "60%")
    })
    
    output$ButtonStat1 <- renderUI({
      radioButtons("Stat1", label = "Statistics", choices = list("GSEA" = 1, "SEA" = 2), selected = input$Stat)
    })
    
    tableGgo <- eventReactive(input$StartGO, {
      gse <- dataIDKegg() 
      
      ggo = groupGO(gene     = gse,
                    OrgDb    = org.Hs.eg.db,
                    level    = input$level,
                    readable = TRUE)
    })
    
    testplot <- function(){
      ggo <- tableGgo()
      barplot(
        ggo, 
        order=TRUE, 
        drop=TRUE, 
        showCategory = input$category, 
        options = list(dom = 'Bfrtip', buttons = c('csv', 'excel', 'pdf'))
      )
    }
    
    output$GroupGO = renderPlot({
      testplot()
    })
    
    output$GoTermPlot <- downloadHandler(
      filename=function(){
        paste("barplot","png",sep=".")
      },content=function(file){
        png(file)
        print(testplot())
        dev.off() 
      }
    )
    
    output$GOID = renderDataTable({
        ggo <- tableGgo()
        newggo <- ggo[order(-ggo$Count)]
        datatable(
          newggo[,1:4],
          rownames = F,
          options = list(
            lengthMenu = list(
              c(10, 25, 50, 100,-1),
              list('10', '25', '50', '100', 'All')
            ),
            dom = 'Bfrtip',
            buttons = c('csv', 'excel', 'pdf')
          ),
          extensions = 'Buttons',
          escape = F
        )
      })
    
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
    
    output$ButtonStat2 <- renderUI({
      radioButtons("Stat2", label = "Statistics", choices = list("GSEA" = 1, "SEA" = 2), selected = input$Stat)
    })
    
    ## Met les ID avec le nom des pathways
    pathwayTable <- aggregate(.~ind,stack(pathway),paste,collapse=' ')
    
    ## Récupère uniquement l'ID sans "hsa"
    IDpathway <- reactive({
      data = dataComplet()
      id <- originGeneID()
      idpathway <- as.character(subset(pathwayTable, values==input$pathwayID)$ind)
      id <- substring(idpathway, 4)
      if (input$GeneID == 1) {
        matrixFC <- matrix(data=data[,3],ncol=1, dimnames=list(c(data[,1]), c()), byrow = TRUE)
        pathwaytOut <- pathview(gene.data = matrixFC, pathway.id = id, species = "hsa")
      }
      if (input$GeneID == 2) {
        dataEnsembl = data[,1]
        h_sapiens = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
        idEnsemblNcbi<- getBM(attributes=c('ensembl_gene_id', 'entrezgene'), filters = 'ensembl_gene_id', values = dataEnsembl, mart = h_sapiens)
        idEnsemblNcbi <-na.exclude(idEnsemblNcbi)
        for (i in 1:nrow(data)) {
          for (j in 1:nrow(idEnsemblNcbi)) {
            if (data[i,1] == idEnsemblNcbi[j,1]) {
              idEnsemblNcbi[j,3] = data[i,3]
            }
          }
        }
        matrixFC <- matrix(data=idEnsemblNcbi[,3], ncol=1, dimnames=list(c(idEnsemblNcbi[,2]), c()), byrow = TRUE)
        pathwaytOut <- pathview(gene.data = matrixFC, pathway.id = id, species = "hsa")
      }
    })
    
    PathwayImage <- reactive ({
      data = dataComplet()
      idpathway <- as.character(subset(pathwayTable, values==input$pathwayID)$ind)
      id <- substring(idpathway, 4)
      if (input$GeneID == 2) {
        dataEnsembl = data[,1]
        h_sapiens = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
        idEnsemblNcbi<- getBM(attributes=c('ensembl_gene_id', 'entrezgene'), filters = 'ensembl_gene_id', values = dataEnsembl, mart = h_sapiens)
        idEnsemblNcbi <-na.exclude(idEnsemblNcbi)
        for (i in 1:nrow(data)) {
          for (j in 1:nrow(idEnsemblNcbi)) {
            if (data[i,1] == idEnsemblNcbi[j,1]) {
              idEnsemblNcbi[j,3] = data[i,3]
            }
          }
        }
        matrixFC <- matrix(data=idEnsemblNcbi[,3], ncol=1, dimnames=list(c(idEnsemblNcbi[,2]), c()), byrow = TRUE)
        file <- paste("hsa",id,".pathview.png",sep="")
        list(src = file)
      }
      else {
        matrixFC <- matrix(data=data[,3],ncol=1, dimnames=list(c(data[,1]), c()), byrow = TRUE)
        file <- paste("hsa",id,".pathview.png",sep="")
        list(src = file)
      }
    })
    
    output$PathwayEnrichment <- renderImage({
      IDpathway()
      PathwayImage()
    })
    
    
    ##################################################
    ## Cinquième page Protein Enrichment
    ##################################################
    
    output$sliderPadj <- renderUI({ 
      numericInput("pAdjID1", label = h4("q-value (filter table)"), value = input$qValueID, min = 0, max = 1, step = 0.01, width = "60%")
    })
    
    SEA_stat2 <- reactive({
      if (input$Stat1 == 1) stat = "GSEA"
      else stat = "SEA"
    })
    

    ajustement <- eventReactive(input$StartProteine, {
      if (input$Ajust == 1) ajust = "holm"
      if (input$Ajust == 2) ajust = "hochberg"
      if (input$Ajust == 3) ajust = "hommel"
      if (input$Ajust == 4) ajust = "bonferroni"
      if (input$Ajust == 5) ajust = "BH"
      if (input$Ajust == 6) ajust = "BY"
      else ajust = "fdr"
    })
    
    testStat <- eventReactive(input$StartProteine, {
      if (input$TestStat == 1) test = "Chi2"
      else test = "Fischer"
    })
    
    output$contents_gene <- renderDataTable({
      dataGene()
    })
    
    tableauprot <- reactive({
      data <- dataComplet()
      stat2 <-SEA_stat2()
      if (stat2 == "SEA"){
        data <- data[which(data[,5] >= input$qValueIDSEA),]
      }
      ajustement <- ajustement()
      requete_genome <- requeteGenome()
      id <- originGeneID()
      if (id == "ncbiID") ID = "entrezgene"
      else ID = "ensembl_gene_id"
      testStatistique <- testStat()
      h_sapiens = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
      
      #requete avec le jeu de données
      data_gene_id <- data[1]
      requete_data <- getBM(attributes=c(ID, 'interpro'), filters = ID, values = data_gene_id, mart = h_sapiens)
      nb_gene_data <- nrow(data_gene_id)
      occurences <- table(requete_data[2])
      mauvaise_ligne <- which(rownames(occurences) == "")
      occurences <- occurences[-mauvaise_ligne] #supprime la ligne avec "" comme nom de domaine
      tableau_data <- as.data.frame(occurences)
      
      #requete genome 
      data_genome_id <- getBM(attributes = c(ID), mart = h_sapiens)
      nb_gene_genome <- nrow(data_genome_id)
      occ_genome <- table(requete_genome[2])
      mauvaise_ligne_genome <- which(rownames(occ_genome) == "")
      occ_genome <- occ_genome[-mauvaise_ligne_genome]
      tableau_genome <- as.data.frame(occ_genome)
      
      # requete interpro vers description
      requete_description <- getBM(attributes=c("interpro","interpro_description"), filters = 'interpro', values = tableau_data[,1], mart = h_sapiens)
      
      # Construction du tableau final
      cpt <- 1
      tableau_final <- matrix(data="NA", nrow=nrow(tableau_data), ncol=5, dimnames=list(c(), c("Domain ID", "Description", "Effectif", "pvalue", "padj")), byrow = TRUE)
      
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
        tableau_final[cpt,2] = requete_description[cpt,2]
        tableau_final[cpt,3] = as.character(tableau_data[cpt,2])
        tableau_final[cpt,4]= round(pvalue_test,3)
        tableau_final[cpt,5]= 0
        
        cpt = cpt + 1
      }
      
      # Calcul de la pvalue adj pour le test du chi2
      cpt = 1
      liste_adj_test = round(p.adjust(liste_test, method = ajustement),3)
      for (x in 1:length(liste_adj_test)){
        tableau_final[cpt, 5] = liste_adj_test[cpt]
        cpt = cpt + 1
      }
      
      taille <- length(which(tableau_final[,5] >= input$pAdjID1))
      tableau_final_final <- matrix(
        data="NA", 
        nrow=taille,  
        ncol=5, 
        dimnames=list(c(), c("Domain ID", "Description", "Effectif", "pvalue", "padj")), byrow = TRUE)
      
      cpt = 1
      for (x in 1:nrow(tableau_final)){
        padj <- tableau_final[x,5]
        if (padj >= input$pAdjID1 ) {
          tableau_final_final[cpt,1] = as.character(tableau_final[x,1])
          tableau_final_final[cpt,2] = as.character(tableau_final[x,2])
          tableau_final_final[cpt,3] = as.character(tableau_final[x,3])
          tableau_final_final[cpt,4]= as.character(tableau_final[x,4])
          tableau_final_final[cpt,5]= as.character(tableau_final[x,5])
          cpt = cpt + 1
        }
      }
      tableau_final_final
    })
    
    output$domain_ID <- renderDataTable({
      tableauprot()
    })
    
    output$plotdomain_ID <- renderPlotly({
      tableau = tableauprot()
      df = data.frame(domainID=tableau[,1], effectif=tableau[,3])
      ggplot(data=df, aes(x=domainID, y=effectif, fill=domainID)) + geom_bar(stat="identity") + coord_flip() #+ theme_minimal()
    })
  }
)