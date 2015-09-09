################## Annotations of Cohn Integration site genes ##########
# Started 7May15
require(biomaRt)
require(stringr)
require(dplyr)
require(reshape2)
require(HGNChelper)
setwd("J:/MacLabUsers/Claire/Projects/HIV-integration")

#read in Mald data to get terms for comparison
load("MalderelliData.formatted.likeSCRIData.Rda")

#here are the terms that I want to check against the annotations
#that I get for the Malderelli list
termsToCheck<-colnames(MalderelliData.formatted.likeSCRIData[,9:48])


#read in spreadsheet of Cohn integration info
CohnData<-read.csv("Cohn et al Integration list.csv")

#select just the symbol.isoform column
symbols.isoforms<-select(CohnData, Symbol.Isoform)
#make into characters
symbols.isoforms$Symbol.Isoform<-as.character(symbols.isoforms$Symbol.Isoform)

#symbol and isoform are separated by |
#need to escape the |
symbols.isoforms<-str_split_fixed(symbols.isoforms$Symbol.Isoform,"\\|",2)
colnames(symbols.isoforms)<-c("symbol","isoform")


#just want the first column (the symbols)
CohnDataCL<-cbind(CohnData,symbols.isoforms)
CohnDataCL<-select(CohnDataCL,-Symbol.Isoform)

hgncCheck<-checkGeneSymbols(CohnDataCL$symbol,
                    unmapped.as.na = FALSE)

CohnDataCL<-cbind(hgncCheck,CohnDataCL)

#keep the correct names
CohnDataCL<-CohnDataCL[,-c(1:2,11)]
#change the col name
names(CohnDataCL)[1]<-"symbol"

#create a mapping from hgnc symbol to go id, strand and entrezid
#based on hsapiens_gene_ensembl information in ensembl mart

load("ensemblMart11May15.Rda")
#set up mart
#ensembl<-useMart("ensembl",dataset="hsapiens_gene_ensembl")

#get go id, strand 

load("BPCohnAnnotations_08May15.Rda")
# CohnAnnotations<-getBM(attributes=c("hgnc_symbol",
#                                 "go_id","namespace_1003"),
#                    mart = ensembl, values=CohnDataCL$symbol,
#                    filters = "hgnc_symbol")
# 
# BPCohnAnnotations<-CohnAnnotations%>%
#   filter(namespace_1003=="biological_process")%>%
#   select(-namespace_1003)
# 
# #eliminating any rows with  missing goid or symbol(none)
# BPCohnAnnotations <- BPCohnAnnotations[
#   nchar(BPCohnAnnotations$hgnc_symbol)>0 &
#     nchar(BPCohnAnnotations$go_id)>0,]


#convert the annotation df (with repeated symbols)to a list of gene symbols symbols
# and the go ids annotated to them.

#the "form" argument says which variable should be the sub elements
#of the list and which should be the names of the elements.

SymbolTermList<-unstack(BPCohnAnnotations, form = go_id~hgnc_symbol)



#Look at each gene symbol in the list and put a 1 if that element
#contains a go term from termsToCheck and 0 if it doesn't.There
#should be either a 0 or 1 for each of the 40 terms to check
#under each symbol

TermCheckList<-lapply(SymbolTermList,
                      FUN=function(i)factor(as.integer(termsToCheck %in% i)))


#turn the list into a data frame (makes symbols the columns)
#then transpose the rows and columns using t()
df<-t(data.frame(TermCheckList))


#make the go ids the col names
colnames(df)<-termsToCheck

#it got turned into a matrix so change back to df
df<-as.data.frame(df)

#make the gene symbols a real column instead of row names
#and put it first
df<-df%>%
  mutate("Gene"=rownames(df))%>%
  select(Gene,1:40)

#I could merge in the original info from the MalderelliData.formatted...
#file but I don't trust merge to do it all correctly so I'm
#just saving this


save(df,file="CohnData BP only checked terms 9Sept15.Rda")



