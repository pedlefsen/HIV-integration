
##################Getting new go terms for Malderelli data ################

require(biomaRt)
require(dplyr)
require(stringr)
require(HGNChelper)
require(reshape2)

setwd("J:/MacLabUsers/Claire/Projects/HIV-integration")

#load data spreadsheet
load("MalderelliData.formatted.likeSCRIData.Rda")

#here are the terms that I want to check against the annotations
#that I get for the Malderelli list
termsToCheck<-colnames(MalderelliData.formatted.likeSCRIData[,9:48])

#using HGNChelper to check symbols and provide corrections
hgncCheck<-checkGeneSymbols(MalderelliData.formatted.likeSCRIData$Gene,
                    unmapped.as.na = FALSE)

#if hgnc help can't find a good symbol, it can leave as NA or leave
#the input if unmapped.as.na=FALSE



MalderelliData.formatted.likeSCRIData<-cbind(hgncCheck,MalderelliData.formatted.likeSCRIData)

#check what the wrong ones look like
y<-filter(hgncCheck,Approved==FALSE)
w<-unique(y)

#keep the correct names
MalderelliData.formatted.likeSCRIData<-MalderelliData.formatted.likeSCRIData[,-c(1:2,6)]
#change the col name
names(MalderelliData.formatted.likeSCRIData)[1]<-"Gene"
#fix this weird one
MalderelliData.formatted.likeSCRIData$Gene<-str_replace(MalderelliData.formatted.likeSCRIData$Gene,
                                                        "MARC2 /// MARCH2","MARCH2")
#new gene list for mart has all correct symbols
Genes<-unique(MalderelliData.formatted.likeSCRIData$Gene)




#set up mart

load("ensemblMart11May15.Rda")
# ensembl<-useMart("ensembl",dataset="hsapiens_gene_ensembl")
# 
# #get entrez id, go id, strand, go name and go category
# MalderelliAnnotations<-getBM(attributes=c("hgnc_symbol",
#                                 "go_id","namespace_1003"),
#                    mart = ensembl, values=Genes,
#                    filters="hgnc_symbol")

load("MalderelliAnnotations3Sept15.Rda")


#BP only, then take out the "BP" column

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

#replace the colon in "GO:xxxxxxx" with a period so it matches
#the formatting in the SCRI data

UniqueBPMalderelliAnnotations$go_id<-str_replace(UniqueBPMalderelliAnnotations$go_id,
                                                 ":","\\.")

#convert the annotation df (with repeated symbols)to a list of gene symbols symbols
# and the go ids annotated to them.

#the "form" argument says which variable should be the sub elements
#of the list and which should be the names of the elements.

SymbolTermList<-unstack(UniqueBPMalderelliAnnotations, form = go_id~hgnc_symbol)



#Look at each gene symbol in the list and put a 1 if that element
#contains a go term from termsToCheck and 0 if it doesn't.

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


save(df,file="MalderelliData BP only checked terms 9Sept15.Rda")


symbolCheck<-cbind(unique(hgncCheck))
colnames(symbolCheck)<-c("origNames","Approved","biomaRtInput")
bioMartOut<-unique(MalderelliAnnotations$hgnc_symbol)
biomaRtBPfiltered<-unique(BPMalderelliAnnotations$hgnc_symbol)

save(symbolCheck,file="symbolCheck.Rda")
save(bioMartOut,file="bioMartOutput.Rda")
save(biomaRtBPfiltered, file="bioMartBPfilteredOutput.Rda")
