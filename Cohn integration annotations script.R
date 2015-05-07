################## Annotations of Cohn Integration site genes ##########

require(biomaRt)
require(stringr)
require(dplyr)

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
                                "go_id","strand","namespace_1003"),
                   mart = ensembl, values=symbols,
                   filters = "hgnc_symbol")


save(CohnIntegrationAnnotations, file="CohnIntegrationAnnotations_7May15.Rda")

write.csv(CohnIntegrationAnnotations,"CohnIntegrationAnnotations.csv")





#get the gene universe: this includes hgnc_symbol,entrezgene,
#go_id,strand, namespace_1003 (go term name)

load("geneUniverse7May15.Rda")

# limit the universe to UNIQUE symbols to use
#for names in the overlap vector
UniqueAnnotations<-unique(annotations)

#0's and 1's vector (overlaps):
#which symbols in UniqueAnnotations overlap with those in symbols? (logical)

overlaps <-factor(as.integer(UniqueAnnotations$hgnc_symbol %in% symbols))

CohnIntegrationAnnotations<-UniqueAnnotations[overlaps,]

#assign the names from allHumanGeneSymbols to the 1's and 0's in overlaps
names(overlaps)<-allHumanGeneSymbols













#set up data for GO2Genes argument in topGO:

#make a list mapping GO ids and their corresponding symbols
GOtoGene<- split(annotations_bp$hgnc_symbol,annotations_bp$go_id)


#remove duplicates
GOtoGene<- lapply(GOtoGene, unique)  
save(GOtoGene, file= "CohnIntegrationGOtoGene.Rda")
###################################################################

#Go2Genes argument needs to be a mapping list (that I just made above)
topGOdata <- new("topGOdata",description = "Simple session",
                 ontology = "BP",
                 allGenes = overlaps,
                 nodeSize = 5,
                 annot = annFUN.GO2genes,
                 GO2genes = GOtoGene)
