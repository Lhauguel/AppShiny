#####################
###### BiomaRt ######
#####################

library(biomaRt)

listMarts()
#mart = useMart("ensembl")
#mart=useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
listEnsembl() # list des bdd
ensembl = useEnsembl(biomart="ensembl")

#head(listDatasets(ensembl))
listDatasets(ensembl) #toutes les datasets d'ensemble

h_sapiens = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl") #si on choisit par exemple l'humain

listFilters(h_sapiens)  #si on veut voir ce qu'il y a dans humain

#exemple de resultat
## attributes ce que l'on veut trouver
## filters = recherche par ce mot clé
## valeurs = ce que l'on donne 
## mart = la bdd sur laquelle on va travailler 
getBM(attributes=c('ensembl_gene_id','ensembl_transcript_id', 'ensembl_peptide_id', 'kegg_enzyme', 'pdb', 'pfam'),filters = 'ensembl_gene_id', values = 'ENSG00000223972', mart = h_sapiens)
getBM(attributes=c('ensembl_gene_id', 'ensembl_transcript_id', 'pfam'),filters = 'ensembl_gene_id', values = 'ENSG00000223972', mart = h_sapiens)

getBM(attributes=c('source'),filters = 'source', values='source', mart = h_sapiens)

data <- read.csv("court_data.csv", header = TRUE)
#data <- read.csv("S3_Table.csv", header = TRUE)
gene_id<-data[1]
gene_id #liste des genes du fichier

list_id <- getBM(attributes=c('pdb', 'pfam', 'family', 'superfamily'),filters = 'ensembl_gene_id', values = gene_id, mart = h_sapiens) # listes des id pdb et pfam (de chaque genes du fichier)
list_id
id_pfam <- list_id[2]
id_pfam
id_family <- list_id[3]
#getBM(attributes=c('ensembl_gene_id', 'pfam', 'pdb'), filters = 'with_pfam', values = id_pfam, mart = h_sapiens)

test <- getBM(attributes=c('ensembl_gene_id', 'pfam'), filters = 'ensembl_gene_id', values = gene_id, mart = h_sapiens)
test
nb_gene_data <- nrow(gene_id)
nb_gene_data
occurences <- table(test[2])
occurences
mauvaise_ligne <- which(rownames(occurences) == "")
occurences <- occurences[-mauvaise_ligne]
tableau_data <- as.data.frame(occurences)
tableau_data
#chisq.test(occurences)

liste_gene_genome <- getBM(attributes = c('ensembl_gene_id'), mart = h_sapiens)
#liste_gene_genome
nb_gene_genome <- nrow(liste_gene_genome)
nb_gene_genome
data_genome_pfam <- getBM(attributes=c('ensembl_gene_id', 'pfam'), filters = 'ensembl_gene_id', values = liste_gene_genome, mart = h_sapiens)
data_genome_pfam
occ_genome <- table(data_genome_pfam[2]) #unlist()
mauvaise_ligne_genome <- which(rownames(occ_genome) == "")
occ_genome <- occ_genome[-mauvaise_ligne_genome]
occ_genome
tableau_genome <- as.data.frame(occ_genome)
tableau_genome



cpt <- 0
for (x in 1:nrow(tableau_data)){
  print(paste(" ===> ID: ", as.character(tableau_data[cpt,1])))
  print(paste(as.character(tableau_data[cpt,1]), as.character(tableau_data[cpt,2])))
  ligne_genome <- rownames(subset(tableau_genome, Var1 == as.character(tableau_data[cpt,1])))
  print(paste(as.character(tableau_genome[ligne_genome,1]), as.character(tableau_genome[ligne_genome,2])))

  m <- matrix(c(tableau_data[cpt,2],nb_gene_data,tableau_genome[ligne_genome,2],nb_gene_genome), nrow=2, ncol=2)
  m
  
  res_chisq <- chisq.test(m)
  pvalue_chisq <- res_chisq$p.value
  print(paste("test de chisq: ", pvalue_chisq))
  
  res_fish <- fisher.test(m)
  pvalue_fisher <- res_fish$p.value
  print(paste("test de Fischer: ", pvalue_fisher))
  
  cpt = cpt + 1
}

