# AppShiny
Authors : Benjamin , Lysiane, Pauline, Thomas

AppShiny : Functionnal Enrichment Analysis in RNA-Seq

**Différentes versions utilisées et packages :**

* R 
    * Version : 3.4.4
    * Documentation : 
* Shiny 
    * Version : 1.2.0
    * Documentation : 
* DT 
    * Version : 0.5
    * Documentation : 
* ggplot2
    * Version : 3.1.0
    * Documentation : 
* biomaRt
    * Version : 2.34.2
    * Documentation : https://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/biomaRt.html
* clusterProfiler
    * Version : 3.6.0
    * Documentation : `vignette("clusterProfiler", package="clusterProfiler")`
* pathview
    * Version : 1.18.2
    * Documentation : `browseVignettes("pathview")
    `
Pour les packages : "biomaRt", "clusterProfiler", "pathview" il faut effectuer quelques lignes :

```
source("http://bioconductor.org/biocLite.R")
biocLite("NomPackage")
library(NomPackage)
```


package ("Cairo") à installer pour le vulcano plot

