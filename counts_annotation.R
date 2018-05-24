#Install and load packages
library('rstudioapi') 
library('dplyr')
#1.----Loading packages -----
rm(list=ls())
list.of.packages <- c('DESeq2',
                      'BiocInstaller',
                      'rstudioapi',
                      "AnnotationDbi",
                      "org.Mm.eg.db")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
if(length(new.packages)) biocLite(new.packages)
lapply(list.of.packages, require, character.only = TRUE)


#2.----Data import -----
## read counts
current_path <- getActiveDocumentContext()$path 
setwd(dirname(current_path ))
setwd('data')
files.list<- list.files(pattern = 'Sample') #importing list of files

readfile = function(x){
  df<-read.csv(x,sep='\t',header = F)
  df<-df[,-c(2,3)]
  df<-df[-c(1,2,3,4),]
  colnames(df)<-c('gene',substr(x,1,8))
  return(df)
}
# let's define a function to read owr countdata for each sample and them apply it
list = lapply(files.list,FUN = readfile)
data = Reduce(function(x, y) merge(x, y, by = "gene", all = TRUE), list)
rownames(data)=data$gene
data = data[,-1]

##reading the metadata
metadat = read.csv("metadat.csv",sep = ';')
if(!all(metadat$sample_name == colnames(data))){
  message('Sample names and column names in metadat.csv are not matching!')
}else{message('Sample names and column names in metadat.csv are matching :)')}


#3.----Analysis----
## Preparing input for DESeq2 analysis
#----Differential expression analysis IFNG ----
#Create two groups ignoring the cell line information
dds.ifng <- DESeqDataSetFromMatrix(countData = data,
                              colData = metadat[,-2],
                              design = ~treatment)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds.ifng <- DESeq(dds) #running DESEq2 analysis
IFNG.result <- results(dds.ifng, contrast=c("treatment",'ifng',"control"))
IFNG.result = na.omit(IFNG.result) #omittnig NAs



dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = metadat,
                              design = ~cell_line+treatment+cell_line:treatment)
#filtering low-abundant transcripts that have < than 10 transcripts total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,] 
nrow(dds) # 20321 genes occured more 10 and more times for all samples
#----Differential expression analysis for cell lines----
#estimate the log2fc of genes for each cell line separately
dds <- DESeq(dds)
#-----B16 effect
#main effect
B16.result <- results(dds, contrast=c("treatment",'ifng',"control"))
#DESeq2 uses b16 cell line as a base level, so there is no need in adding an interaction coefficient
B16.result = na.omit(B16.result) #omittnig NAs

#----TC1A9----
# the condition effect for tc1a9 line
# this is, by definition, the main effect *plus* the interaction term
# I've checked the results provided by comparison using a dummy 'group' variable, the amount of significant genes is the same for both lines.
TC1A9.result <- results(dds, list( c("treatment_ifng_vs_control","cell_lineTC1A9.treatmentifng") ))
TC1A9.result = na.omit(TC1A9.result)
#4.----Annotation----
results = list(B16=B16.result,TC1A9=TC1A9.result,IFNG = IFNG.result)
for(i in 1:length(results)){
  print(names(results)[i])
  res = results[[i]]
  res$symbol <- mapIds(org.Mm.eg.db,
                       keys=row.names(res),
                       column="SYMBOL",
                       keytype="ENSEMBL",
                       multiVals="first")
  res$entrez <- mapIds(org.Mm.eg.db,
                       keys=row.names(res),
                       column="ENTREZID",
                       keytype="ENSEMBL",
                       multiVals="first")
  res$ensembl = rownames(res)
  results[[i]]=res
  
}

annotate_genes = function(res){
  res$symbol <- mapIds(org.Mm.eg.db,
                       keys=row.names(res),
                       column="SYMBOL",
                       keytype="ENSEMBL",
                       multiVals="first")
  res$entrez <- mapIds(org.Mm.eg.db,
                       keys=row.names(res),
                       column="ENTREZID",
                       keytype="ENSEMBL",
                       multiVals="first")
  res$ensembl = rownames(res)
  return(res)
}

results = lapply(results,annotate_genes)

#5.----Export----
setwd("../")
if (file.exists('DESeq2_results')){
  setwd('DESeq2_results')
  message('DESeq2_results directory already exists, settig wd...')
}else{
  dir.create('DESeq2_results')
  message('creating DESeq2_results directory, settig wd...')
  setwd('DESeq2_results')
}
for(i in 1:length(results)){
  names(results)[i]
  write.csv(results[[i]],paste0(names(results)[i],'_result.csv'),row.names = F)
  file = file(paste0(names(results)[i],"_info.txt")  )  
  writeLines(c(paste0(names(results)[i],' information:'),
               results[[i]]@elementMetadata$description[1:6]),
             file)
  close(file)
}
#exporting results info
writeLines(c(paste0('Created:'),capture.output(Sys.time()),capture.output(sessionInfo())), "sessionInfo.txt")
