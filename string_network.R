#load packages
library(dplyr)





#------------STRING ANALYSIS---------
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

#create a list with significant genes
getwd()
setwd("/Users/grygoriyzolotarov/interferon_gamma/")
significant = function(result){
  out = result%>%filter(abs(log2FoldChange)>1&padj<=0.01)%>%arrange(padj)%>%dplyr::select(ensembl)%>%pull()
  out = out[!is.na(out)]
  return(out)
}
significant_genes = lapply(desseq2,significant)




#write significant files
write = function(i){
  writeLines(significant_genes[i],sep='\n',"outfile.txt")
}
lapply(seq_along(significant_genes),FUN = function(i) writeLines(paste0(significant_genes[[i]],'\n'), paste0(names(significant_genes)[i],'sign_genes.txt')))


#read STRING output and convert to an igraph object
setwd("/Users/grygoriyzolotarov/interferon_gamma/STRING")
net = read.table('string_interactions.tsv',header = T)




#----STRING analysis in R-----
library(BiocInstaller)
library(STRINGdb)
data(diff_exp_example1)
head(diff_exp_example1)
example1_mapped <- string_db$map( diff_exp_example1, "gene", removeUnmappedRows = TRUE )

species = get_STRING_species(version="10", species_name=NULL)
id = species[grepl('Mus musculus',species[,2]),]$species_id

#load string database for Mus musculus:
string_db <- STRINGdb$new( version="10", species=id,score_threshold=0, input_directory="" )


STRINGdb$help("add_diff_exp_color")
STRINGdb$help("get_pubmed_interaction")
STRINGdb$help("get_png")
data(diff_exp_example1)


int = string_db$get_interactions(example1_mapped[1:100,]$STRING_id)
int = string_db$get_pubmed_interaction(example1_mapped[1:100,]$STRING_id)
graph = string_db$get_graph()


example1_mapped <- string_db$map( diff_exp_example1, "gene", removeUnmappedRows = TRUE )
hits <- example1_mapped$STRING_id[1:200]
string_db$plot_network( hits )
string_db$plot_network(example1_mapped$STRING_id[1:100])


significant_genes = lapply(desseq2,FUN = function(result) result%>%filter(abs(log2FoldChange)>1&padj<=0.01)%>%arrange(padj)%>%dplyr::select(ensembl,entrez,symbol,log2FoldChange,padj) )


ifng = significant_genes$IFNG%>%dplyr::select(padj,log2FoldChange,symbol)
ifng_mapped <- string_db$map( ifng, "symbol", removeUnmappedRows = TRUE )

# plot first 4 clusters
par(mfrow=c(2,2))
for(i in seq(1:4)){
  string_db$plot_network(clustersList_ifng[[i]])
}


library('devtools')


#1.network for ifng
#2.cluster network
#3.gene set enrichment analysis for clusters - which clusters are there?


library(igraph)
interactome = string_db$get_graph()



#get subgraph from graph
nodes = ifng_mapped$STRING_id[ifng_mapped$STRING_id%in%V(interactome)$name]
ifng = induced.subgraph(graph = interactome,nodes)


# Set layout options
l <- layout.fruchterman.reingold(ifng)
l<-layout_with_fr(ifng)
l <- layout_in_circle(ifng)
l <- layout_on_sphere(ifng)

# Plot graph and subgraph
plot.igraph(x=ifng,layout=l,vertex.label.cex = 1e-10,vertex.size = 5,edge.width=0.4)
E(ifng)$combined_score
E(ifng)$width = E(ifng)$combined_score*0.001
plot(ifng)

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

?closeness

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


#plot the neighborhood of selected proteins
V(ifng)$name = gene_hs[match(V(ifng)$name,gene_hs$id),]$symbol
ifng <- set.vertex.attribute(ifng, "label", value=gene_hs[match(V(ifng)$name,gene_hs$id),]$symbol)
n1 = neighbors(ifng, V(ifng)$name[1], 1)
plot(ifng[n],vertex.label.cex = 1e-10) 
V(ifng)$name[n1]









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

