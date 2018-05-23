#load packages
library(dplyr)
library(igraph)
library(STRINGdb)
library(BiocInstaller)




#------------STRING ANALYSIS---------
#1.----Load DESeq2 results ----
setwd("~/ifng2018/DESeq2_results")
files = list.files(pattern = 'csv')
read.file = function(i){
  result = read.csv(files[i])
  name = unlist(strsplit(files[i],'\\.'))[1]
  comment(result) =readLines(list.files(pattern = paste0(unlist(strsplit(name,'_'))[1],'_info')))
  return(result)
}
desseq2 = lapply(seq_along(files),FUN = read.file)
names(desseq2) = sapply(files,function(x) unlist(strsplit(x,'_'))[1],USE.NAMES = F)




#----STRING analysis in R-----
species = get_STRING_species(version="10", species_name=NULL) #selecting species id for Mus musculus
id = species[grepl('Mus musculus',species[,2]),]$species_id

#load string database for Mus musculus:
string_db <- STRINGdb$new( version="10", species=id,score_threshold=0, input_directory="" )
significant_genes = lapply(desseq2,FUN = function(result) result%>%filter(abs(log2FoldChange)>1&padj<=0.01)%>%arrange(padj)%>%dplyr::select(ensembl,entrez,symbol,log2FoldChange,padj) )


ifng = significant_genes$IFNG%>%dplyr::select(padj,log2FoldChange,symbol)
ifng_mapped <- string_db$map( ifng, "symbol", removeUnmappedRows = TRUE ) #map significant genes to STRING identifiers

b16 = significant_genes$B16%>%dplyr::select(padj,log2FoldChange,symbol)
b16_mapped <- string_db$map(b16, "symbol", removeUnmappedRows = TRUE ) #map significant genes to STRING identifiers


#1.network for ifng
#2.cluster network
#3.gene set enrichment analysis for clusters - which clusters are there?


interactome = string_db$get_graph()#convert STRING network into an igraph object


#get subgraph from graph
nodes = ifng_mapped$STRING_id[ifng_mapped$STRING_id%in%V(interactome)$name]
nodes = nodes[1:50]
ifng = induced.subgraph(graph = interactome,nodes)

nodes = b16_mapped$STRING_id[b16_mapped$STRING_id%in%V(interactome)$name]
nodes = nodes[1:50]
b16 = induced.subgraph(graph = interactome,nodes)
#delete vertices with no edges
b16 = delete.vertices(b16,degree(b16)==0)


ifng = delete.vertices(ifng,degree(ifng)==0)

# Set layout options
l <- layout.fruchterman.reingold(ifng)
l<-layout_with_dh(ifng)
l <- layout_in_circle(ifng)
l <- layout_in_circle(ifng)

# Plot graph and subgraph
plot.igraph(x=ifng,layout=l,vertex.label.cex = 1e-10,vertex.size = 5,edge.width=E(ifng)$coexpression/max(E(ifng)$coexpression)*2)
E(ifng)$width = E(ifng)$combined_score*0.001




l <- layout.fruchterman.reingold(b16)
closeness(b16, mode="all", weights=NA) 
deg()
plot.igraph(x=b16,layout=l,vertex.label.cex = 1e-10,vertex.size = 5,edge.width=E(ifng)$coexpression/max(E(ifng)$coexpression)*2)










#clusterisation
ceb <- cluster_edge_betweenness(ifng) 
dendPlot(ceb, mode="hclust")
plot(ceb, ifng,vertex.label.cex = 1e-10) 



#cenrality measurements
deg = degree(ifng, mode="in")

match_genes = function(gene_id){
  if(length(ifng_mapped[ifng_mapped$STRING_id==gene_id,]$symbol)>1){
    return(ifng_mapped[ifng_mapped$STRING_id==gene_id,]$symbol[1])
  }else{
    return(ifng_mapped[ifng_mapped$STRING_id==gene_id,]$symbol)
  }
}
symb = unlist(sapply(names(deg),match_genes))
gene_deg = data.frame(symbol = symb, id = names(symb) ,deg = deg )
gene_deg$logFC = ifng_mapped[match(gene_deg$symbol,ifng_mapped$symbol),3]
gene_deg$padj = ifng_mapped[match(gene_deg$symbol,ifng_mapped$symbol),2]
gene_deg = gene_deg%>%arrange(-deg)



#closeness 
closeness(ifng, mode="all", weights=NA) 
centr_clo(net, mode="all", normalized=T) 


#hubs and authorities
hs <- hub_score(ifng, weights=NA)$vector
as <- authority_score(ifng, weights=NA)$vector


match_genes = function(gene_id){
  if(length(ifng_mapped[ifng_mapped$STRING_id==gene_id,]$symbol)>1){
    return(ifng_mapped[ifng_mapped$STRING_id==gene_id,]$symbol[1])
  }else{
    return(ifng_mapped[ifng_mapped$STRING_id==gene_id,]$symbol)
  }
}
symb = unlist(sapply(names(hs),match_genes))
gene_hs = data.frame(symbol = symb, id = names(symb) ,hs = hs )
gene_hs =gene_hs%>%arrange(-hs)
gene_hs$logFC = ifng_mapped[match(gene_hs$symbol,ifng_mapped$symbol),3]
gene_hs$padj = ifng_mapped[match(gene_hs$symbol,ifng_mapped$symbol),2]
head(gene_hs)

#plot the neighborhood of selected proteins
V(ifng)$name = as.vector(gene_hs[match(V(ifng)$name,gene_hs$id),]$symbol)
node = gene_hs[465,]$symbol



#plot node neighborhood
nei = induced_subgraph(ifng, ego(ifng, 1, node)[[1]])
l <- layout_on_sphere(nei)
l <- layout_with_kk(nei)

E(nei)$width = E(nei)$coexpression/max(E(nei)$coexpression)*5
plot(nei,layout = l)
title(node)
dev.off()






#degree distribution for the network
deg <- degree(ifng, mode="all")
deg.dist <- degree_distribution(ifng, cumulative=T, mode="all")
plot( x=0:max(deg), y=1-deg.dist, pch=19, cex=1.2, col="orange", xlab="Degree", ylab="Cumulative Frequency")



deg <- degree(b16, mode="all")

deg.dist <- degree_distribution(b16, cumulative=T, mode="all")
plot( x=0:max(deg), y=1-deg.dist, pch=19, cex=1.2, col="orange", xlab="Degree", ylab="Cumulative Frequency")


hs <- hub_score(b16, weights=NA)$vector
match_genes = function(gene_id){
  if(length(b16_mapped[b16_mapped$STRING_id==gene_id,]$symbol)>1){
    return(b16_mapped[b16_mapped$STRING_id==gene_id,]$symbol[1])
  }else{
    return(b16_mapped[b16_mapped$STRING_id==gene_id,]$symbol)
  }
}

symb = unlist(sapply(names(hs),match_genes))
gene_hs = data.frame(symbol = symb, id = names(symb) ,hs = hs )
gene_hs =gene_hs%>%arrange(-hs)
gene_hs$logFC = b16_mapped[match(gene_hs$symbol,b16_mapped$symbol),3]
gene_hs$padj = b16_mapped[match(gene_hs$symbol,b16_mapped$symbol),2]
gene_hs%>%arrange(-logFC)














#Clusterisation

nodes = ifng_mapped$STRING_id[ifng_mapped$STRING_id%in%V(interactome)$name]
nodes[1:50]
ifng = induced.subgraph(graph = interactome,nodes)


ceb <- cluster_edge_betweenness(ifng)
plot(ceb, ifng)



#Stat1 appears as a core element with hub score = 1 
hist(deg)
deg.dist <- degree_distribution(ifng, cumulative=T, mode="all")
plot( x=0:max(deg), y=1-deg.dist, pch=19, cex=1.2, col="orange",xlab="Degree", ylab="Cumulative Frequency")

mean_distance(ifng, directed=T)




# filter by p-value and add a color column
# (i.e. green down-regulated gened and red for up-regulated genes)
example1_mapped_pval05 <- string_db$add_diff_exp_color( subset(example1_mapped, pvalue<0.05), logFcColStr="logFC" )
payload_id <- string_db$post_payload( example1_mapped_pval05$STRING_id,colors=example1_mapped_pval05$color )
string_db$plot_network( hits, payload_id=payload_id )




nodes <- read.csv("Dataset1-Media-Example-NODES.csv", header=T, as.is=T)
links <- read.csv("Dataset1-Media-Example-EDGES.csv", header=T, as.is=T)

links <- aggregate(links[,3], links[,-3], sum)

links <- links[order(links$from, links$to),]

colnames(links)[4] <- "weight"

rownames(links) <- NULL
net <- graph_from_data_frame(d=links, vertices=nodes, directed=T) 
net <- simplify(net, remove.multiple = F, remove.loops = T) 

net
plot(net, edge.arrow.size=.4,vertex.label=NA)

net









#----Network analysis----
library('igraph')

#compute the centrality for each node
node_degree = net%>%
  group_by(node1)%>%
  summarize(n=n())%>%
  arrange(-n)


#transform out data into an igraph object:

