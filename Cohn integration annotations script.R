################## Annotations of Cohn Integration site genes ##########
# Started 7May15
require(biomaRt)
require(stringr)
require(dplyr)
setwd("J:/MacLabUsers/Claire/Projects/HIV-integration")
#read in spreadsheet of integration info

CohnIntegrations<-read.csv("Cohn et al Integration list.csv")

#select just the symbol.isoform column
symbols.isoforms<-select(CohnIntegrations, Symbol.Isoform)
#make into characters
symbols.isoforms$Symbol.Isoform<-as.character(symbols.isoforms$Symbol.Isoform)

#symbol and isoform are separated by |
#need to escape the |
symbols.isoforms<-str_split_fixed(symbols.isoforms$Symbol.Isoform,"\\|",2)

#just want the first column (the symbols)
symbols<-symbols.isoforms[,1]

#create a mapping from hgnc symbol to go id, strand and entrezid
#based on hsapiens_gene_ensembl information in ensembl mart

#set up mart
ensembl<-useMart("ensembl",dataset="hsapiens_gene_ensembl")

#get entrez id, go id, strand 
CohnIntegrationAnnotations<-getBM(attributes=c("hgnc_symbol","entrezgene",
                                "go_id","strand","namespace_1003","name_1006"),
                   mart = ensembl, values=symbols,
                   filters = "hgnc_symbol")


save(CohnIntegrationAnnotations, file="CohnIntegrationAnnotations_7May15.Rda")

write.csv(CohnIntegrationAnnotations,"CohnIntegrationAnnotations.csv")

