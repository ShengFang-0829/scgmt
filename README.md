## scgmt
# Description
This is an R pack exercise work, aimed at integrating multiple single cell rank-based gene set enrichment analysis methods and integrating the results into a standardized format for storage at metedata at Seurat object.
# Installation
Before installing the scgmt package, you need to install the following dependent R packages.
~~~
  VISION,
  Seurat,
  BiocParallel,
  AUCell,
  GSEABase,
  GSVA,
  UCell,
  singscore,
  magrittr,
  dplyr,
  tidyr,
  clusterProfiler,
  scales,
  methods,
  purrr,
  reshape2,
  stringr,
  ggplot2,
  pheatmap,
  ggridges,
  viridis
~~~
Install scgmt package:
~~~
devtools::install_github("ShengFang-0829/scgmt")
~~~
## Usage
Here, we use the pbmc3k dataset that comes with SeuratData.

# Prepare Data
~~~
library(SeuratData)
InstallData("pbmc3k")
data("pbmc3k")
rds <- pbmc3k.final
~~~

# Main calculation process
~~~
library(scgmt)
# Just provide the gmt file path for storing the gene list for each signaling pathway.
# The gmt file in the example was obtained from this website address:https://www.gsea-msigdb.org/gsea/msigdb/index.jsp
signatures <- "h.all.v2023.1.Hs.symbols.gmt"

# Single cell rank-based gene set enrichment analysis methods integration in the scgmt function.
# The current methods include:"AUCell","UCell","AddModuleScore","gsva","ssgsea","zscore","plage","VISION","JAS_likelihood","JAS_oddsratio" and "singscore".Fill in using the method parameter.
# The signatures parameter is filled in the gmt file path.

rds1 <- scgmt(rds=rds,method="UCell",signatures=signatures)


# The gene set scores of each signaling pathway in each cell are stored in the metadata of the Seurat object.
colnames(rds1@meta.data)
~~~

# Visualization
Visualization includes the following functions:   
scgmt_line_plot   
scgmt_merge_line_plot       
scgmt_heatmap_plot  
scgmt_scatter_plot  
scgmt_ridges_plot   
scgmt_density_plot  
~~~
p1 <- scgmt_line_plot(rds1,signatures="HALLMARK_HYPOXIA",group.by="orig.ident")
print(p1)

p2 <- scgmt_merge_line_plot(rds1,signatures=c("HALLMARK_HYPOXIA","HALLMARK_MITOTIC_SPINDLE"))
print(p2)

p3 <- scgmt_heatmap_plot(rds1,signatures=c("HALLMARK_TNFA_SIGNALING_VIA_NFKB","HALLMARK_HYPOXIA","HALLMARK_CHOLESTEROL_HOMEOSTASIS"))
print(p3)

p4 <- scgmt_scatter_plot(rds1,signature.x="HALLMARK_TNFA_SIGNALING_VIA_NFKB",signature.y="HALLMARK_HYPOXIA")
print(p4)

p5 <- scgmt_ridges_plot(rds1,signature="HALLMARK_HYPOXIA")
print(p5)

p6 <- scgmt_density_plot(rds1,signature="HALLMARK_HYPOXIA")
print(p6)
~~~

## Contact information
For any questions please contact fngseng12345@163.com or https://github.com/ShengFang-0829/scgmt/issues
