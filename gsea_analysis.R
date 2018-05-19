getwd()
files = list.files(pattern = 'gsea')

#Find common and different pathways for both cell lines
#Load gsea results:
read_file = function(file){
  data = read.csv(file)
  colnames(data) = sapply(colnames(data),function(x) unlist(strsplit(x,'\\.'))[2],USE.NAMES = F)
  return(data)
}
gsea = lapply(files,read_file)
names(gsea)=sapply(files,function(x) unlist(strsplit(x,'\\.'))[1],USE.NAMES = F)


#Compute similarities and dissimilarities:
outersect <- function(x, y) {
  sort(c(setdiff(x, y),
         setdiff(y, x)))
}
#kegg
kegg = gsea[grepl('kegg',names(gsea))]
list1= kegg$B16gsea_kegg
list2 = kegg$TC1A9gsea_kegg
sim = intersect(list1$ID,list2$ID)
diff= outersect(list1$ID,list2$ID)

list1$line = 'b16'
list2$line = 'tc'

data = rbind(list1,list2)
View(data[data$ID%in%diff,])
View(data[data$ID%in%sim,])

reactome = gsea[grepl('reactome',names(gsea))]
list1= reactome$B16gsea_reactome
list2 = reactome$TC1A9gsea_reactome
sim = intersect(list1$ID,list2$ID)
diff= outersect(list1$ID,list2$ID)

list1$line = 'b16'
list2$line = 'tc'
data = rbind(list1,list2)
View(data[data$ID%in%diff,])
View(data[data$ID%in%sim,])



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


desseq2$B16%>%filter(symbol=='Ccl5')

