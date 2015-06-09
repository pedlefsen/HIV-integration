require(biomaRt)
require(dplyr)
require(stringr)
setwd("J:/MacLabUsers/Claire/Projects/HIV-integration")

#load data spreadsheet
load("MalderelliData.formatted.likeSCRIData.Rda")

#give the data frame an easier name
MalderelliData<-MalderelliData.formatted.likeSCRIData

Genes<-unique(MalderelliData$Gene)

######################## different method start here #################
library(org.Hs.eg.db)
library(topGO)

#make a named list of all GO ids in org.Hs.eg.db and the symbols that belong to it
Go2Gene <-annFUN.org("BP", mapping = "org.Hs.eg.db",
                             ID = "symbol")
#make into a df
Go2Genedf<-stack(Go2Gene)

#rename cols
colnames(Go2Genedf)<-c("symbol","GO_id")

#extract the just the annotated symbols that I also have in my symbol list
annotated<-Go2Genedf[Go2Genedf$symbol %in% Genes,]

#how many unique GO ids are there?
length(unique(annotated$GO_id))
#2596

########### COMPARE WITH  biomaRt METHOD ###########################

load("BPMalderelliAnnotations_08May15.Rda")

length(unique(BPMalderelliAnnotations$go_id))
#2926




