
##################Getting new go terms for Malderelli data ################

require(biomaRt)
require(dplyr)
require(stringr)
setwd("J:/MacLabUsers/Claire/Projects/HIV-integration")

#load data spreadsheet
load("MalderelliData.formatted.likeSCRIData.Rda")

#give the data frame an easier name
MalderelliData<-MalderelliData.formatted.likeSCRIData

Genes<-unique(as.character(MalderelliData$Gene))

#set up mart

load("ensemblMart11May15.Rda")
ensembl<-useMart("ensembl",dataset="hsapiens_gene_ensembl")




#get entrez id, go id, strand, go name and go category
MalderelliAnnotations<-getBM(attributes=c("hgnc_symbol",
                                "go_id","namespace_1003"),
                   mart = ensembl, values=Genes,
                   filters="hgnc_symbol")

load("MalderelliAnnotations8June15.Rda")

#BP only, then take out the BP column

BPMalderelliAnnotations<-MalderelliAnnotations%>%
  filter(namespace_1003=="biological_process")%>%
  select(-namespace_1003)

#eliminating any rows with  missing goid or symbol(none)
BPMalderelliAnnotations <- BPMalderelliAnnotations[
  nchar(BPMalderelliAnnotations$hgnc_symbol)>0 &
    nchar(BPMalderelliAnnotations$go_id)>0,]

#make sure entries are unique
UniqueBPMalderelliAnnotations<-unique(BPMalderelliAnnotations)
#same # of observations as the non-unique version


#convert the annotation df (with repeated symbols)to a list of go ids
# and the symbols that are in them

BPMalderelliGoAndSymbolList<-unstack(BPMalderelliAnnotations)



## For each gene symbol that I put into biomaRt (Genes),
# is that symbol present or absent from each named element in GOAndSymbolList?

LogicalBPMalderelliGOAndSymbolList<-lapply(BPMalderelliGoAndSymbolList,
          FUN=function(i)factor(as.integer(MalderelliData$Gene %in% i)))


#make that list into a data frame with GO ids as the columns

dfBPMalderelliGOAndSymbol<-data.frame(LogicalBPMalderelliGOAndSymbolList)



#combine with other ID columns from original data set
CLMalderelliData.formatted.likeSCRIData<-cbind(MalderelliData[,1:8],dfBPMalderelliGOAndSymbol)


#look back in annotation list to check that things are aligned 

head(BPMalderelliGoAndSymbolList)

#so in rows RAB1A and ULK2, column GO:00000045 should !=0 since
#I know from BPMalderelliGOAndSymbolList that both symbols are in that
#go id
CLMalderelliData.formatted.likeSCRIData[CLMalderelliData.formatted.likeSCRIData$"GO.0000045"!=0,c("Gene","GO.0000045")]

CLMalderelliData.formatted.likeSCRIData[CLMalderelliData.formatted.likeSCRIData$Gene=="RAB1A","GO.0000045"]






################################   SAVE ####################################


save(CLMalderelliData.formatted.likeSCRIData,
     file = "CLMalderelliData.formatted.likeSCRIData.Rda")




######################## Filtering for Hallmark terms ##############

hallmark<-read.csv("hallmarksCancer_GO_CL.csv")
hallmark$go_id<-as.character(hallmark$go_id)

c<-unique(BPMalderelliAnnotations$go_id %in% hallmark$go_id)

#filter out just the hallmark go_ids that overlapped with those that I
#already found

hallmarkFiltered<-dfBPMalderelliGOAndSymbol[,c]


CLMalderelliHallmarkFiltered<-cbind(MalderelliData[,1:8],hallmarkFiltered)
