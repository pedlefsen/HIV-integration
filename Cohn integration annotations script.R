################## Annotations of Cohn Integration site genes ##########
# Started 7May15
require(biomaRt)
require(stringr)
require(dplyr)
require(reshape2)
setwd("J:/MacLabUsers/Claire/Projects/HIV-integration")
#read in spreadsheet of integration info

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




#create a mapping from hgnc symbol to go id, strand and entrezid
#based on hsapiens_gene_ensembl information in ensembl mart

#set up mart
ensembl<-useMart("ensembl",dataset="hsapiens_gene_ensembl")

#get go id, strand 
CohnAnnotations<-getBM(attributes=c("hgnc_symbol",
                                "go_id","namespace_1003"),
                   mart = ensembl, values=CohnDataCL$symbol,
                   filters = "hgnc_symbol")

BPCohnAnnotations<-CohnAnnotations%>%
  filter(namespace_1003=="biological_process")%>%
  select(-namespace_1003)

#eliminating any rows with  missing goid or symbol(none)
BPCohnAnnotations <- BPCohnAnnotations[
  nchar(BPCohnAnnotations$hgnc_symbol)>0 &
    nchar(BPCohnAnnotations$go_id)>0,]


save(BPCohnAnnotations, file="BPCohnAnnotations_08May15.Rda")

#convert the annotation df (with repeated symbols)to a list of go ids
# and the symbols that are in them

BPCohnGoAndSymbolList<-unstack(BPCohnAnnotations)


save(BPCohnGoAndSymbolList, file = "BPCohnGoAndSymbolList_08May15.Rda")


#make a logical vector indicating which symbols are in go ids.For each go id,
#need a 0 (not in that group) or 1 (in that group) for each symbol in the list


x<- lapply(BPCohnGoAndSymbolList,
           FUN=function(i)factor(as.integer(CohnDataCL$symbol %in% i)))

save(x,file = "BPGOAndSymbolListLogical.Rda")

#make a data frame out of the list with one column per go id and symbols as
#rownames. values are 0's and 1's
#http://stackoverflow.com/questions/18747800/fast-way-of-converting-large-list-to-dataframe
nCohn<-length(x[[1]])
dfCohn<-structure(x, row.names = c(NA, -n), class = "data.frame")


CohnData.formatted.likeSCRIData<-cbind(CohnDataCL,df)

colnames(CohnData.formatted.likeSCRIData)[1:3]<-c("Pt.","Read","Strand")
#change read to char vector
CohnData.formatted.likeSCRIData$Read<-as.character(CohnData.formatted.likeSCRIData$Read)

#combine read and site to match mald format
CohnData.formatted.likeSCRIData$Read<-paste(CohnData.formatted.likeSCRIData$Read,
                                            CohnData.formatted.likeSCRIData$site,
                                            sep="+")

#remove site column
CohnData.formatted.likeSCRIData<-select(CohnData.formatted.likeSCRIData,-site)


#reorder columns to match mald

#make cols in correct order and rename as needed
colsIwant<-CohnData.formatted.likeSCRIData%>%
  select(Pt., Read, symbol,Strand)

names(colsIwant)[3]<-"Gene"

#remove offending cols
CohnData.formatted.likeSCRIData<-CohnData.formatted.likeSCRIData%>%
  select(-c(Pt.,Read,Strand,symbol))

#re add corrected cols
CohnData.formatted.likeSCRIData<-cbind(colsIwant,CohnData.formatted.likeSCRIData)

save(CohnData.formatted.likeSCRIData, file="CohnData.formatted.likeSCRIData.Rda")


#This should be ILR and KPNA1
CohnData.formatted.likeSCRIData[CohnData.formatted.likeSCRIData$"GO:0000018"!=0,
                                c("Gene","GO:0000018")]

#ok
