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
biocLite('STRINGdb')
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
"get_bioc_graph"
STRINGdb$help("get_png")
int = string_db$get_interactions(example1_mapped[1:100,]$STRING_id)
int = string_db$get_pubmed_interaction(example1_mapped[1:100,]$STRING_id)
string_network = string_db$get_graph()
string_db$get_png(example1_mapped[1:100,]$STRING_id)

head(diff_exp_example1)
example1_mapped <- string_db$map( diff_exp_example1, "gene", removeUnmappedRows = TRUE )
hits <- example1_mapped$STRING_id[1:200]
string_db$plot_network( hits )
string_db$plot_network(example1_mapped$STRING_id[1:100])


significant_genes = lapply(desseq2,FUN = function(result) result%>%filter(abs(log2FoldChange)>1&padj<=0.01)%>%arrange(padj)%>%dplyr::select(ensembl,entrez,symbol,log2FoldChange,padj) )

b16 = significant_genes$B16%>%dplyr::select(padj,log2FoldChange,symbol)
b16_mapped <- string_db$map( b16, "symbol", removeUnmappedRows = TRUE )
hits <- b16_mapped$STRING_id[1:50]
clustersList_b16 <- string_db$get_clusters(b16_mapped$STRING_id[1:600])
# plot first 4 clusters
par(mfrow=c(2,2))
for(i in seq(1:4)){
  string_db$plot_network(clustersList_b16[[i]])
}

ifng = significant_genes$IFNG%>%dplyr::select(padj,log2FoldChange,symbol)
ifng_mapped <- string_db$map( ifng, "symbol", removeUnmappedRows = TRUE )
hits <- ifng_mapped$STRING_id[1:50]
clustersList_ifng <- string_db$get_clusters(ifng_mapped$STRING_id[1:600])
# plot first 4 clusters
par(mfrow=c(2,2))
for(i in seq(1:4)){
  string_db$plot_network(clustersList_ifng[[i]])
}

string_db$get_graph()
library('devtools')

#dowloading Mitools package:
install_github("vitkl/MItools")
library(MItools)
??fullInteractome
#1.network for ifng
#2.cluster network
#3.gene set enrichment analysis for clusters - which clusters are there?





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

