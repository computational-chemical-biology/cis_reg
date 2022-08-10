if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos='http://cran.us.r-project.org')

BiocManager::install(c("Rgraphviz", "png", "KEGGgraph", "org.Hs.eg.db",
                       "gage", "pathview", "DESeq2", "tximport", "universalmotif"))

install.packages(c("readr", "stringr", "ggplot2", "pheatmap", "venn", 
                   "IRkernel", "gplots"), repos='http://cran.us.r-project.org')
IRkernel::installspec()  # to register the kernel in the current R installation


