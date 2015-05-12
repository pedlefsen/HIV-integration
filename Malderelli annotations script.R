
##################Getting new go terms for Malderelli data ################

require(biomaRt)

setwd("J:/MacLabUsers/Claire/Projects/HIV-integration")

#load data spreadsheet
load("MalderelliData.formatted.likeSCRIData.Rda")

#give the data frame an easier name
MalderelliData<-MalderelliData.formatted.likeSCRIData


#set up mart

ensembl<-useMart("ensembl",dataset="hsapiens_gene_ensembl")


save(ensembl, file = "ensemblMart11May15")


#get entrez id, go id, strand, go name and go category
MalderelliAnnotations<-getBM(attributes=c("hgnc_symbol",
                                "go_id","namespace_1003"),
                   mart = ensembl, values=MalderelliData$Gene,
                   filters="hgnc_symbol")

#BP only, then take out the BP column

BPMalderelliAnnotations<-MalderelliAnnotations%>%
  filter(namespace_1003=="biological_process")%>%
  select(-namespace_1003)

#eliminating any rows with  missing goid or symbol(none)
BPMalderelliAnnotations <- BPMalderelliAnnotations[
  nchar(BPMalderelliAnnotations$hgnc_symbol)>0 &
    nchar(BPMalderelliAnnotations$go_id)>0,]


#convert the annotation df (with repeated symbols)to a list of go ids
# and the symbols that are in them

BPMalderelliGoAndSymbolList<-unstack(BPMalderelliAnnotations)

#what is the length of each element (how many genes in each go id)
x<-sapply(BPMalderelliGoAndSymbolList,FUN=length)
#how many have <1 gene in them?
sum(x<1) # none

LogicalBPMalderelliGOAndSymbolList<- lapply(BPMalderelliGoAndSymbolList,
                                  FUN=function(i)factor(as.integer(MalderelliData$Gene %in% i)))



#http://stackoverflow.com/questions/18747800/fast-way-of-converting-large-list-to-dataframe
#make a data frame out of the list with one column per go id and symbols as
#rownames. values are 0's and 1's

n<-length(LogicalBPMalderelliGOAndSymbolList[[1]])

dfBPMalderelliGOAndSymbol<-structure(LogicalBPMalderelliGOAndSymbolList,
                                     row.names = c(NA, -n), class = "data.frame")

#change colons to periods
colnames(dfBPMalderelliGOAndSymbol)<-str_replace_all(colnames(dfBPMalderelliGOAndSymbol),":","\\.")


CLMalderelliData.formatted.likeSCRIData<-cbind(MalderelliData[,1:8],dfBPMalderelliGOAndSymbol)

###################   SAVE ####################################


save(CLMalderelliData.formatted.likeSCRIData,
     file = "CLMalderelliData.formatted.likeSCRIData.Rda")












#look back in annotation list to check that things are aligned 

head(BPMalderelliGoAndSymbolList)

#so in rows RAB1A and ULK2, column GO:00000045 should !=0 since
#I know from BPMalderelliGOAndSymbolList that both symbols are in that
#go id
CLMalderelliData[CLMalderelliData$"GO:0000045"!=0,c("Gene","GO:0000045")]

CLMalderelliData[CLMalderelliData$Gene=="RAB1A","GO:0000045"]



