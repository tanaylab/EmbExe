
packages <- c(
    "Matrix",
    "metacell",
    "tgstat",
    "tgutil",
    "tglkmeans",
    "RColorBrewer",
    "ggrepel",
    "zoo",
    "umap",
    "tidyverse",
    "Seurat",
    "DoubletFinder",
    "slanter",
    "here",
    "devtools",
    "princurve",
    "qvalue",
    "data.table"
)

for (pkg in packages) {
    library(pkg, character.only = T)
}
