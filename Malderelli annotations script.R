
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

#get entrez id, go id, strand, go name and go category
MalderelliAnnotations<-getBM(attributes=c("hgnc_symbol","entrezgene",
                                "go_id","strand","namespace_1003", "name_1006"),
                   mart = ensembl, values=GeneSymbols,
                   filters="hgnc_symbol")



save(MalderelliAnnotations, file="MalderelliAnnotations_7May15.Rda")


