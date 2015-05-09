
#################### Adding GO terms to Malderelli data ################
require(biomaRt)

setwd("J:/MacLabUsers/Claire/Projects/HIV-integration")

#load data spreadsheet
load("MalderelliData.formatted.likeSCRIData.Rda")

#give the data frame an easier name
MalderelliData<-MalderelliData.formatted.likeSCRIData

GeneSymbols<-as.character(MalderelliData$Gene) #gene symbols



#set up mart
ensembl<-useMart("ensembl",dataset="hsapiens_gene_ensembl")
save(ensembl, file = "ensemblMart08May15")

#get entrez id, go id, strand, go name and go category
MalderelliAnnotations<-getBM(attributes=c("hgnc_symbol","entrezgene",
                                "go_id","strand","namespace_1003", "name_1006"),
                   mart = ensembl, values=GeneSymbols,
                   filters="hgnc_symbol")

#BP only, then take out the BP column

BPMalderelliAnnotations<-MalderelliAnnotations%>%
  filter(namespace_1003=="biological_process")%>%
  select(-namespace_1003)

#eliminating any rows with  missing goid or symbol(none)
BPMalderelliAnnotations <- BPMalderelliAnnotations[
  nchar(BPMalderelliAnnotations$hgnc_symbol)>0 &
    nchar(BPMalderelliAnnotations$go_id)>0,]


save(BPMalderelliAnnotations, file="BPMalderelliAnnotations_08May15.Rda")

Short_BPMalderelliAnnotations = unique(BPMalderelliAnnotations[,c(1,3)])

BPMalderelliGoAndSymbolList<-unstack(Short_BPMalderelliAnnotations)

save(BPMalderelliGoAndSymbolList, file = "BPMalderelliGoAndSymbolList_08May15.Rda")


