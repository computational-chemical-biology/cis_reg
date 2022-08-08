if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("Rgraphviz", "png", "KEGGgraph", "org.Hs.eg.db",
                       "gage", "pathview", "DESeq2", "tximport"))

install.packages(c("readr", "stringr", "ggplot2", "pheatmap", "venn", 'IRkernel'))
IRkernel::installspec()  # to register the kernel in the current R installation


