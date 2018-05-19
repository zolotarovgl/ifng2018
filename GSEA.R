#Here the GSEA and GAGE are performed using GO,KEGG and Rectome gene sets
rm(list=ls()) 
library(BiocInstaller)
#1.----Load packages----
list.of.packages <- c("ggplot2",
                      'clusterProfiler',
                      'ReactomePA',
                      'reactome.db',
                      'DOSE',
                      'bindrcpp',
                      'ggplotly',
                      'dplyr',
                      'org.Mm.eg.db')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
if(length(new.packages)) biocLite(new.packages)
lapply(list.of.packages, require, character.only = TRUE)
rm(list = c('list.of.packages','new.packages'))

#1.----Load DESeq2 results ----
setwd("~/interferon_gamma/DESeq2_results")
files = list.files(pattern = 'csv')
read.file = function(i){
  result = read.csv(files[i])
  name = unlist(strsplit(files[i],'\\.'))[1]
  comment(result) =readLines(list.files(pattern = paste0(unlist(strsplit(name,'_'))[1],'_info')))
  return(result)
}
desseq2 = lapply(seq_along(files),FUN = read.file)
names(desseq2) = sapply(files,function(x) unlist(strsplit(x,'_'))[1],USE.NAMES = F)

#2.----Define functions----
compute_group_go = function(sign.genes,level=1){
  list = list()
  ont = c('BP','MF','CC')
  for(i in 1:length(ont)){
    message(paste0('groupGO...',ont[i]))
    start_time <- Sys.time()
    ggo<- groupGO(gene     = sign.genes,
                  OrgDb    = org.Mm.eg.db,
                  ont      = ont[i],
                  level    = level,
                  readable = TRUE)
    list[ont[i]] = ggo
    end_time <- Sys.time()
    print(end_time - start_time)
  }
  return(list)
}

compute_enrich_go = function(sign.genes,level=1){
  list = list()
  ont = c('BP','MF','CC')
  for(i in 1:length(ont)){
    message(paste0('groupGO...',ont[i]))
    start_time <- Sys.time()
    ggo<- enrichGO(gene     = sign.genes,
                   OrgDb    = org.Mm.eg.db,
                   ont      = ont[i],
                   readable = TRUE)
    list[ont[i]] = ggo
    end_time <- Sys.time()
    print(end_time - start_time)
  }
  return(list)
}
compute_gsea_go = function(genes,ont = ont,pval = 0.05){
  list = list()
  ont = c('BP','MF','CC')
  for(i in 1:length(ont)){
    message(paste0('gseaGO...',ont[i]))
    start_time <- Sys.time()
    gse<- gseGO(geneList     = genes,
                OrgDb        = org.Mm.eg.db,
                ont          = ont[i],
                nPerm        = 1000,
                minGSSize    = 10,
                maxGSSize    = 2000,
                pvalueCutoff = pval,
                verbose      = T)
    list[ont[i]] = gse
    end_time <- Sys.time()
    print(end_time - start_time)
  }
  return(list)
}

#2.----Gene Enrichment and GSEA: clusterProfiler----
#----Cluster Profiler ----
#----GO enrichment----
enrich_go = list()
for(i in 1:length(names(desseq2))){
  print(names(desseq2)[i])
  #list of signifint genes and data
  result = desseq2[[i]]
  genes = result$log2FoldChange
  names(genes) = as.character(result$entrez)
  genes = sort(na.omit(genes))
  genes = sort(genes[!is.na(names(genes))],decreasing = T)
  sign.genes = result%>%
    filter(abs(log2FoldChange)>=1&padj<=0.01)
  sign.genes = as.character(sign.genes$entrez)
  enrich_go[[i]] = compute_enrich_go(sign.genes = sign.genes)
  names(enrich_go)[i] = names(desseq2)[i]
}
names(enrich_go) = unlist(strsplit(names(enrich_go),'_'))[seq(1,length(names(enrich_go))*2,2)]
#----GO gsea----
gsea_go = list()
for(i in 1:length(names(desseq2))){
  print(names(desseq2)[i])
  #list of signifint genes and data
  result = desseq2[[i]]
  genes = result$log2FoldChange
  names(genes) = as.character(result$entrez)
  genes = sort(na.omit(genes))
  genes = sort(genes[!is.na(names(genes))],decreasing = T)
  sign.genes = result%>%
    filter(abs(log2FoldChange)>=1&padj<=0.01)
  sign.genes = as.character(sign.genes$entrez)
  gsea_go[i] = compute_gsea_go(genes = genes, ont = ont)
  names(gsea_go)[i] = names(desseq2)[i]
}
names(gsea_go) = names(desseq2)
gsea_go = lapply(gsea_go,function(list) setReadable(list,OrgDb        = org.Mm.eg.db))
View(gsea_go$B16[,])
#----KEGG enrichment----
enrich_kegg = list()
for(i in 1:length(names(desseq2))){
  print(names(desseq2)[i])
  #list of signifint genes and data
  result = desseq2[[i]]
  genes = result$log2FoldChange
  names(genes) = as.character(result$entrez)
  genes = sort(na.omit(genes),decreasing = T)
  genes = sort(genes[!is.na(names(genes))],decreasing = T)
  sign.genes = result%>%
    filter(abs(log2FoldChange)>=1&padj<=0.01)
  sign.genes = as.character(sign.genes$entrez)
  enrich_kegg[i] = enrichKEGG(gene    = sign.genes,
                         organism     = 'mmu',
                         keyType = 'ncbi-geneid',
                         minGSSize    = 10,
                         pvalueCutoff = 0.05)
  names(enrich_kegg)[i] = names(desseq2)[i]
}
names(enrich_kegg) = unlist(strsplit(names(desseq2),'_'))[seq(1,length(names(desseq2))*2,2)]
#----KEGG gsea----
gsea_kegg = list()
for(i in 1:length(names(desseq2))){
  print(names(desseq2)[i])
  #list of signifint genes and ata
  result = desseq2[[i]]
  genes = result$log2FoldChange
  names(genes) = as.character(result$entrez)
  genes = sort(na.omit(genes),decreasing = T)
  genes = sort(genes[!is.na(names(genes))],decreasing = T)
  sign.genes = result%>%
    filter(abs(log2FoldChange)>=1&padj<=0.01)
  sign.genes = as.character(sign.genes$entrez)
  gsea_kegg[i] = gseKEGG(geneList     = genes,
                         organism     = 'mmu',
                         nPerm        = 3000,
                         keyType = 'ncbi-geneid',
                         minGSSize    = 10,
                         pvalueCutoff = 0.05,
                         verbose      = T)
  names(gsea_kegg)[i] = names(desseq2)[i]
}
names(gsea_kegg) = names(desseq2)

gsea_kegg = lapply(gsea_kegg,function(list) setReadable(list,OrgDb = org.Mm.eg.db))
setReadable(gsea_kegg$B16,OrgDb = org.Mm.eg.db)
#----Reactome enrichment----
enrich_reactome = list()
for(i in 1:length(names(desseq2))){
  print(names(desseq2)[i])
  #list of signifint genes and data
  result = desseq2[[i]]
  genes = result$log2FoldChange
  names(genes) = as.character(result$entrez)
  genes = sort(na.omit(genes),decreasing = T)
  genes = sort(genes[!is.na(names(genes))],decreasing = T)
  sign.genes = result%>%
    filter(abs(log2FoldChange)>=1&padj<=0.01)
  sign.genes = na.omit(sign.genes$entrez)
  enrich_reactome[i] = enrichPathway(sign.genes, pvalueCutoff=0.05,organism = 'mouse',
                                   pAdjustMethod="BH", readable=T)
  names(enrich_reactome)[i] = names(desseq2)[i]
}

#----Reactome gsea----

gsea_reactome = list()
for(i in 1:length(names(desseq2))){
  print(names(desseq2)[i])
  #list of signifint genes and data
  result = desseq2[[i]]
  genes = result$log2FoldChange
  names(genes) = as.character(result$entrez)
  genes = sort(na.omit(genes),decreasing = T)
  genes = sort(genes[!is.na(names(genes))],decreasing = T)
  sign.genes = result%>%
    filter(abs(log2FoldChange)>=1&padj<=0.01)
  sign.genes = as.character(sign.genes$entrez)
  gsea_reactome[i] = gsePathway(genes, nPerm=3000, organism = 'mouse',
                                  minGSSize=30, pvalueCutoff=0.05,
                                  pAdjustMethod="BH", verbose=T)
  names(gsea_reactome)[i] = names(desseq2)[i]
}
names(gsea_reactome) = names(desseq2)
gsea_reactome= lapply(gsea_reactome,function(list) setReadable(list,OrgDb = org.Mm.eg.db,keytype = 'auto'))
#3.----Export Results----
setwd("~/interferon_gamma/gsea")

#GO
lapply(seq_along(enrich_go),FUN = function(j) lapply(seq_along(enrich_go[[j]]),FUN = function(i) write.csv(enrich_go[[j]][i],paste0(names(enrich_go)[j],'_',names(enrich_go[[j]][i]),'_enrich_go.csv'),row.names = F)))
#KEGG
lapply(seq_along(enrich_kegg),FUN = function(i) write.csv(enrich_kegg[i],paste0(names(enrich_kegg[i]),'enrich_kegg.csv'),row.names = F))
lapply(seq_along(gsea_kegg),FUN = function(i) write.csv(gsea_kegg[i],paste0(names(gsea_kegg[i]),'gsea_kegg.csv'),row.names = F))
#Reactome
lapply(seq_along(enrich_reactome),FUN = function(i) write.csv(enrich_reactome[i],paste0(names(enrich_reactome[i]),'enrich_reactome.csv'),row.names = F))
lapply(seq_along(gsea_reactome),FUN = function(i) write.csv(gsea_reactome[i],paste0(names(gsea_reactome[i]),'gsea_reactome.csv'),row.names = F))



