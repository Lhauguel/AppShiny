#####################
###### BiomaRt ######
#####################

library(biomaRt)

listMarts()
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

data <- read.csv("court_data.csv", header = TRUE)
gene_id<-data[1]
gene_id #liste des genes du fichier

list_id <- getBM(attributes=c('pdb', 'pfam', 'family', 'superfamily'),filters = 'ensembl_gene_id', values = gene_id, mart = h_sapiens) # listes des id pdb et pfam (de chaque genes du fichier)
list_id
id_pfam <- list_id[2]
id_pfam
id_family <- list_id[3]
getBM(attributes=c('ensembl_gene_id', 'pfam', 'pdb'), filters = 'with_pfam', values = id_pfam, mart = h_sapiens)

test <- getBM(attributes=c('ensembl_gene_id', 'pfam'), filters = 'ensembl_gene_id', values = gene_id, mart = h_sapiens)
test
(test)

nb_gene_total <- nrow(gene_id)
nb_gene_total
occurences <- table(unlist(gene_id))
occurences

for (i in table(unlist(gene_id))){
  print(i)
}