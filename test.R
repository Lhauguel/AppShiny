#####################
###### BiomaRt ######
#####################

library(biomaRt)

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
getBM(attributes=c('ensembl_gene_id','ensembl_transcript_id','hgnc_symbol','hgnc_id'),filters = 'ensembl_gene_id', values = 'ENSG00000223972', mart = h_sapiens)

