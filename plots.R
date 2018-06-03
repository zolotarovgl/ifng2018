library('pheatmap')
library(AnnotationDbi)
install.packages('pheatmap')


vsd <- vst(dds.ifng, blind=FALSE)
rld <- rlog(dds.ifng, blind=FALSE)


anno <- as.data.frame(metadat[,-1])
rownames(anno) = metadat$sample_name
pheatmap(assay(vsd)[1:5,], cluster_rows=T, show_rownames=FALSE,
         cluster_cols=T,annotation_col = anno)
desseq2$IFNG%>%arrange(abs(log2FoldChange))%>%head

#heatmap for significant genes:
topgenes = desseq2$IFNG%>%filter(padj<0.05)%>% arrange(padj,-abs(log2FoldChange))%>%head(100)%>%select(ensembl)%>%pull()






pheatmap(assay(rld)[topgenes,], cluster_rows=T, show_rownames=FALSE,
         cluster_cols=T,annotation_col = anno,annotation_row = )
assay(dds.ifng)

#log2Fold changes for significant genes for each cell line
b16  = desseq2$B16%>%filter(padj<0.05)%>%select(symbol,ensembl,log2FoldChange)
rownames(b16) = b16$ensembl
tc1a9  = desseq2$TC1A9%>%filter(padj<0.05)%>%select(symbol,ensembl,log2FoldChange)
rownames(tc1a9) = tc1a9$ensembl

changes <- merge(b16, tc1a9,by='row.names', all=TRUE)
rownames(changes) = changes$Row.names
colnames(changes)[c(4,7)] = c('lfcb16','lfctc')
changes = changes[,-1]
data = data.frame(b16 = changes$lfcb16,
                  tc1a9 = changes$lfctc,
                  row.names = rownames(changes))


#heatmap for significant genes:
up = desseq2$IFNG%>%filter(padj<0.05)%>% arrange(-log2FoldChange)%>%head(10)%>%select(ensembl)%>%pull()
down = desseq2$IFNG%>%filter(padj<0.05)%>% arrange(-log2FoldChange)%>%tail(10)%>%select(ensembl)%>%pull()
topgenes = c(as.character(up),as.character(down))


library(reshape)

pvalb16 = function(gene){
  return(desseq2$B16%>%filter(ensembl == gene)%>%select(padj)%>%pull)
}
pvalb16('ENSMUSG00000001741')
data$gene = rownames(data)
ggdata = melt(data,id = 'gene',measure.vars = c('b16','tc1a9'))
colnames(ggdata)[2] = 'line'
ggdata$symbol <- mapIds(org.Mm.eg.db,
                     keys=ggdata$gene,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")

plot = ggdata[ggdata$gene%in%topgenes,]
plot$symbol = factor(plot$symbol,levels = unique(plot$symbol))
plot$pval = rep(0.5,nrow(plot))
ggplot(plot, aes(x = line, y= symbol)) +
  geom_tile(aes(fill = value), color = "black") +
  geom_text(aes(label = pval),cex=3)+
  ylab('Gene')+
  xlab('Cell line')+
  scale_fill_gradient2(low = "darkblue",midpoint=0,mid ='white', high = "red",name = 'log2 Change')+
  theme(panel.border = element_blank(),
        panel.background = element_blank()
        )






b16 = desseq2$B16%>%filter(padj<0.05)%>%select(symbol,ensembl,padj,log2FoldChange)
b16 = desseq2$B16[desseq2$B16$padj<0.05,c(2,6,7,9)]
tc = desseq2$TC1A9%>%filter(padj<0.05)%>%select(symbol,ensembl,padj,log2FoldChange)
rownames(tc)=NULL
b16$line = 'b16'
tc$line = 'tc1a9'
desseq2$B16$line = 'B16'
desseq2$TC1A9$line = 'TC1A9'
data = rbind(desseq2$B16,desseq2$TC1A9)
data = data%>%select(symbol,ensembl,padj,log2FoldChange,line)
data$padj = signif(data$padj, digits = 2)
plot = data[data$ensembl%in%topgenes,]
plot$symbol= factor(plot$symbol,levels = unique(plot$symbol))
ggplot(plot, aes(x = line, y= symbol)) +
  geom_tile(aes(fill = log2FoldChange), color = "black") +
  geom_text(aes(label = padj),cex=3)+
  ylab('Gene')+
  xlab('Cell line')+
  geom_text()
  ggtitle('Expression changes')+
  scale_fill_gradient2(low = "darkblue",midpoint=0,mid ='white', high = "red",name = 'log2 Change')+
  theme(panel.border = element_blank(),
        panel.background = element_blank()
  )
colnames(data)



View()



# gsea heatmap
View(gsea$IFNGgsea_go)
ids = c('regulation of cell migration','regulation of cell adhesion','I-kappaB kinase/NF-kappaB signaling')




gsea.b16 = gsea[grepl('B16',names(gsea))]
gsea.tc = gsea[grepl('TC1A9',names(gsea))]
gsea.ifng = gsea[grepl('IFNG',names(gsea))]

gsea.grouped = lapply(c('B16','TC1A9','IFNG'),function(name) return(gsea[grepl(name,names(gsea))]))
names(gsea.grouped) = c('B16','TC1A9','IFNG')

res = lapply(gsea.grouped,col1)
names(gsea.grouped$B16)
col1 = function(list){
  col2 = function(i){
    list[i]$type = names(list)[i]
  }
  lapply(seq_along(list),col2 )
}
list = gsea.grouped$B16


col = function(i){
  df = list[i]
  df$type = names(list)[[i]]
  return(df)
}
res=lapply(seq_along(list),col)



gsea$IFNGgsea_reactome
View(gsea$IFNGgsea_go)
View(gsea$IFNGgsea_kegg)


View(gsea$B16gsea_reactome)
View(gsea$B16gsea_go)
View(gsea$B16gsea_kegg)


View(gsea$TC1A9gsea_reactome)
View(gsea$TC1A9gsea_go)
View(gsea$TC1A9gsea_kegg)

go = gsea[grepl('go',names(gsea))]
list1= go$B16gsea_go
list2 = go$TC1A9gsea_go
sim = intersect(list1$ID,list2$ID)
diff= outersect(list1$ID,list2$ID)

list1$line = 'b16'
list2$line = 'tc'

data = rbind(list1,list2)
View(data[data$ID%in%diff,])
View(data[data$ID%in%sim,])
data%>%filter(Description%in%ids)%>%View()


sb16=desseq2$B16%>%filter(abs(log2FoldChange)>1&padj<0.01)
stc=desseq2$TC1A9%>%filter(abs(log2FoldChange)>1&padj<0.01)
length(intersect(sb16$ensembl,stc$ensembl))



length(intersect(desseq2$B16$ensembl,desseq2$TC1A9$ensembl))/(length(desseq2$B16$ensembl) + length(desseq2$TC1A9$ensembl))
b16 = desseq2$B16 %>% filter(abs(log2FoldChange)>1&padj<0.01)
tc = desseq2$TC1A9 %>% filter(abs(log2FoldChange)>1&padj<0.01)
length(tc$ensembl)-length(intersect(b16$ensembl,tc$ensembl))
length(intersect(b16$ensembl,tc$ensembl))/(length(tc$ensembl)+length(b16$ensembl))



gsea$IFNGgsea_reactome%>%filter(Description == 'Regulation of DNA replication')
View(gsea$IFNGgsea_reactome)
'M/G1 Transition'
gsea$IFNGgsea_reactome$Description[1]
  #PCA plot
rld <- rlog(dds)
data <- plotPCA(rld, intgroup=c("treatment", "cell_line"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))

ggplot(data, aes(PC1, PC2, color=treatment,shape = cell_line)) +
  geom_point(size=3) +
  ggtitle("Principal Components")+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_bw()

