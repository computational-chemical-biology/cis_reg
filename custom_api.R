#* Project             RNA seq: differential expression analysis with DESeq2
#* Version information   V 10.0
#* Date                  16 out 2019
#* authors               Luisa C. Nora, Murilo H. A. Cassiano
#*                       
#* Copyright notice      This is the pipeline to get abundance files, import 
#*                       then using txinport and do the DEA with DESeq2
pkgs <- c("tximport", "readr", "stringr", "DESeq2", "ggplot2", 
	  "RColorBrewer", "pheatmap", "coda.base", "venn", "gage", 
	  "pathview", "rjson", "universalmotif", "gplots", "hrbrthemes")

suppressMessages(lapply(pkgs, require, character.only = TRUE)) 
#requires the DEG result matrix and the chart title
volcano_plot = function (mtx, title) {  #vulcano_plot(cane_filter, "sugar")
	pvaluetr = (-1*log(mtx$padj))
	pvalueT= pvaluetr[!is.na(pvaluetr)]
	fold_change = (mtx$log2FoldChange)
	fold_change = fold_change[!is.na(pvaluetr)]
	plot( (fold_change), (pvalueT),
      		type='n',
      		ylab='-1*log(p-value)',
      		xlab='Log2FoldChange',
	        xlim=c(-12, 12),
	        ylim=c(0, 700),
      		main=title)
	points(fold_change, pvalueT, col='black')
	points( fold_change[(pvalueT>3&fold_change>1)]  , pvalueT[(pvalueT>3&fold_change>1)]  ,col='red',pch=16)
	points( fold_change[(pvalueT>3&fold_change< -1)]  , pvalueT[(pvalueT>3&fold_change< -1)]  ,col='blue',pch=16)
	abline(h=3)
	abline(v=1)
	abline(v=-1)
}

#requires a dds object: deseq2 analized
sample_check = function (dds_object) {   #sample_check(dds)
	vsd <- vst(dds_object, blind=FALSE)
	sampleDists <- dist(t(assay(vsd)))
	sampleDistMatrix <- as.matrix(sampleDists)
	rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$samples, sep="-")
	colnames(sampleDistMatrix) <- NULL
	colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
	pheatmap(sampleDistMatrix,
        	 clustering_distance_rows=sampleDists,
         	 clustering_distance_cols=sampleDists,
         	 col=colors
	)
}

#requires the DEG result matrix
filter_DE_genes = function (res) {    #cane_filter = filter_DE_genes(res_cane)
	a = res[(res$padj <= 0.05 & !is.na(res$padj)), ]
	return(a[a$log2FoldChange >= 1 | a$log2FoldChange <= -1, ])
}

#read abundance files and make it ready to deseq analysis
cast_deseq2 = function (abundance_dir, sample_cond, tx_to_gene) {
	filess <- list.files(abundance_dir)
	files <- file.path(abundance_dir, filess)
	names(files) <- paste0(str_replace(filess, "_abundance.tsv", ""))
	samples <- read.table(file.path(sample_cond), sep = ",", header=TRUE, stringsAsFactors = TRUE)
	tx2gene <- read.csv(file.path(tx_to_gene), sep = ",", header = TRUE)
	txi <- tximport(files, type="kallisto", tx2gene=tx2gene)
	head(txi$counts)
	a = DESeqDataSetFromTximport(txi, colData = samples, design = ~ condition)
	keep <- rowSums(counts(a)) >= 10
	a <- a[keep,]
	head(a)
	return(a)
}

func_analy = function(base = "kog", data_filter = NULL, marg = c(4.25,20,0,1), legX = 0.5 , legY = 0.5, limi, colores = c("blue","red")){
  file_base = " "
  ann = NULL
  formated = NULL
  up_genes = row.names(data_filter[data_filter$log2FoldChange > 0,])
  down_genes = row.names(data_filter[data_filter$log2FoldChange < 0,])
  print("UP-regulated:"); print(length(up_genes));
  print("DOWN-regulated:"); print(length(down_genes))
  print("All-degs:"); print(length(down_genes)+length(up_genes))
  file_base = switch(base, "kog" = "aux_data/annotation/Rhoto_IFO0880_4_GeneCatalog_proteins_20170509_KOG.tab",
                     "kegg" = "aux_data/annotation/Rhoto_IFO0880_4_GeneCatalog_proteins_20170509_KEGG.tab", 
                     "go" = "aux_data/annotation/Rhoto_IFO0880_4_GeneCatalog_proteins_20170509_GO.tab")
  ann = read.csv(file = file_base, header = T, sep = "\t")
  if (base == "kog") {
    func_up = ann$kogClass[(ann$transcriptId %in% up_genes)]
    func_down = ann$kogClass[(ann$transcriptId %in% down_genes)]
    print("Func_UP-regulated:"); print(length(func_up));
    print("Func_DOWN-regulated:"); print(length(func_down))
    formated = as.data.frame(summary(func_up))
    formated = cbind(formated, as.data.frame(summary(func_down)))
    formated = t(formated)
    print(formated)
    print(colnames(formated))
    row.names(formated)= c("Up-regulated", "Down-regulated")
    formated = formated[ , order(colnames(formated), decreasing = T)]
    formated = formated[c(2,1),]
    par(mar = marg) 
    barplot(formated, names.arg = names(formated), horiz = T, las=1, col = colores, xlab = "Number of genes", beside = T, xlim = limi)
    legend(x = legX, y = legY, row.names(formated), bty = "n", fill = colores)
  }else{
    lkup = read.csv(file = "aux_data/annotation/Rtoruloides_lookup_table_transcriptId.csv", header = T, sep = ",")[,1:2]
    for (i in 1:length(ann$X.proteinId)) {
      idx = which(ann$X.proteinId[i] == lkup$proteinId)
      ann$transcriptId[i] = lkup$transcriptId[idx]
    }
    if(base == "kegg"){
      func_up = ann$pathway_class[(ann$transcriptId %in% up_genes)]
      func_down = ann$pathway_class[(ann$transcriptId %in% down_genes)]
      formated = as.data.frame(summary(func_up))
      formated = cbind(formated, as.data.frame(summary(func_down)))
      formated = t(formated)
      #formated = formated[ , order(names(formated))]
      row.names(formated)= c("Up-regulated", "Down-regulated")
      par(mar=marg) 
      barplot(formated, names.arg = names(formated), horiz = T, las=1, col = colores, xlab = "Number of genes", beside = T)
      legend(x = legX, y = legY, row.names(formated), bty = "n", fill = colores)
    }
    if(base == "go"){
      func_up = ann$goName[(ann$transcriptId %in% up_genes)]
      func_down = ann$goName[(ann$transcriptId %in% down_genes)]
      formated = as.data.frame(summary(func_up))
      formated = cbind(formated, as.data.frame(summary(func_down)))
      formated = t(formated)
      #formated = formated[ , order(names(formated))]
      row.names(formated)= c("Up-regulated", "Down-regulated")
      par(mar=marg) 
      barplot(formated[,1:25], names.arg = names(formated[,1:25]), horiz = T, las=1, col = colores, xlab = "Number of genes", beside = T)
      legend(x = legX, y = legY, row.names(formated), bty = "n", fill = colores)
      par(mar=c(4.25,25,0,10)) 
      barplot(formated[,26:50], names.arg = names(formated[,26:50]), horiz = T, las=1, col = colores, xlab = "Number of genes", beside = T)
      legend(x = legX, y = legY, row.names(formated), bty = "n", fill = colores)
      par(mar=c(4.25,25,0,10)) 
      barplot(formated[,51:75], names.arg = names(formated[,51:75]), horiz = T, las=1, col = colores, xlab = "Number of genes", beside = T)
      legend(x = legX, y = legY, row.names(formated), bty = "n", fill = colores)
      par(mar=c(4.25,25,0,10)) 
      barplot(formated[,76:100], names.arg = names(formated[,76:100]), horiz = T, las=1, col = colores, xlab = "Number of genes", beside = T)
      legend(x = legX, y = legY, row.names(formated), bty = "n", fill = colores)
    }
  }
}

gage_kegg_view = function(foldChange, path.dir = "/home/", name = "my_pathway", custon_id = NULL){
  resFC = foldChange["log2FoldChange"]
  indexes = match(row.names(foldChange), koala$V1)
  koala= (koala[indexes,])
  lilo = cbind(resFC, koala$V2)
  #row.names(lilo) = koala$V2
  print("biomol")
  #lilo = lilo[lilo$"koala$V2" != "",]
  lilo <- lilo[!is.na(lilo[,2]),]
  row.names(lilo) = lilo[,2] 
  lilo$"koala$V2" = NULL
  #testing the names' pattern, as recommended
  head(kegg.gs[[1]]); head(rownames(lilo))
  #running
  foldChange.kegg.p <- gage(as.matrix(lilo), gsets = kegg.gs)
  system(paste("mkdir ",path.dir,"/",name, sep=""))
  setwd(paste(path.dir,"/",name, sep=""))
  write.csv(foldChange.kegg.p$greater, file = "greater.csv")
  write.csv(foldChange.kegg.p$less, file = "less.csv")
  write.csv(foldChange.kegg.p$stats, file = "stats.csv")
  write.csv(foldChange.kegg.p$greater[(foldChange.kegg.p$greater[,"p.val"] <= 0.05 & !is.na(foldChange.kegg.p$greater[,"p.val"])),],file="greater_filter.csv")
  write.csv(foldChange.kegg.p$less[(foldChange.kegg.p$less[,"p.val"] <= 0.05 & !is.na(foldChange.kegg.p$less[,"p.val"])),],file = "less_filter.csv")
  #this is the custon anaysis
  if (!is.null(custon_id)) {
    system(paste("mkdir ",path.dir,"/",name,"/custom", sep = ""))
    setwd(paste(path.dir,"/",name,"/custom", sep = ""))
    pathview(gene.data=lilo, pathway.id=custon_id, species="ko", low = "blue" , mid = "white", high = "red", res = 400)
  }
  #here we attempt to draw all pathways with statistically relevance confidence (from gage results)
  great_p = substr(  row.names(  foldChange.kegg.p$greater[ (foldChange.kegg.p$greater[,"p.val"] <= 0.05 & !is.na(foldChange.kegg.p$greater[,"p.val"])) ,  ]  ) , 1, 7  )
  less_p =  substr(  row.names(  foldChange.kegg.p$less[ (foldChange.kegg.p$less[,"p.val"] <= 0.05 & !is.na(foldChange.kegg.p$less[,"p.val"])) ,  ]  )  , 1, 7  )
  # :) pathwayviewing
  system(paste("mkdir ",path.dir,"/",name,"/greater", sep = ""))
  setwd(paste(path.dir,"/",name,"/greater", sep = ""))
  pathview(gene.data=lilo,pathway.id=great_p[great_p != "ko04723"],species="ko",low = "blue" , mid = "white", high = "red")#,res = 400)
  system(paste("mkdir ",path.dir,"/",name,"/less", sep = ""))
  setwd(paste(path.dir,"/",name,"/less", sep = ""))
  pathview(gene.data=lilo, pathway.id=less_p, species="ko", low = "blue" , mid = "white", high = "red", res = 400)
}
