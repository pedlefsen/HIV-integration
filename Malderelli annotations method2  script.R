
##################### NOTES ###################################
# 
#This is what this script does:
#Start with MalderelliData.formatted.likeSCRIdata.Rda from Paul/Charles
# Use the HGNChelper package to check for any incorrect/excel mogrified symbols and correct them.
# Use the or.Hs.eg.db database (Genome wide annotation for Human, primarily based on mapping using Entrez Gene identifiers) whole genome annotations.
# 
# Use topGo "classic" algorithm for Fisher's test to see which, if any, of the Go term annotations of 
# the Malderelli gene list were in different proportions to what is expected from the whole genome 
# annotations from org.Hs.eg.db
# 
# Out of all the Malderelli symbols that were annotated, make a list showing which belong to the 
# significant annotations (0 if no, 1 if yes). Convert to a data frame and add other columns to format 
# the df the same way as SCRI data.

require(dplyr)
require(stringr)
library(HGNChelper)
setwd("J:/MacLabUsers/Claire/Projects/HIV-integration")

#load data spreadsheet
load("MalderelliData.formatted.likeSCRIData.Rda")


#using HGNChelper to check symbols and provide corrections
r<-checkGeneSymbols(MalderelliData.formatted.likeSCRIData$Gene,
                    unmapped.as.na = FALSE)
r<-r[,c(1,3)]
names(r)[1]<-"Gene"
r<-unique(r)

#I tried to merge without making the values in r unique and it 
#gave me ~25k entries in the result...


#merge in the correct "suggested symbols" with my Malderelli symbol list
MalderelliData.formatted.likeSCRIData<-merge(r,MalderelliData.formatted.likeSCRIData,
                                             by ="Gene")

MalderelliData.formatted.likeSCRIData<-MalderelliData.formatted.likeSCRIData[,2:262]

names(MalderelliData.formatted.likeSCRIData)[1]<-"Gene"

MalderelliData.formatted.likeSCRIData$Gene<-str_replace(MalderelliData.formatted.likeSCRIData$Gene,
                                                        "MARC2 /// MARCH2","MARCH2")
Genes<-unique(MalderelliData.formatted.likeSCRIData$Gene)

######################## different method start here #################
library(org.Hs.eg.db)
library(topGO)
library(stringr)
library(reshape2)

#Make a list of all the symbols in org.Hs.eg.db and under each
#symbol, put the GO ids that annotate that symbol. The symbols are
#the "names" of the list items (the GO ids)
Gene2GO <-inverseList(annFUN.org("BP", mapping = "org.Hs.eg.db",
                                 ID = "symbol"))

#Get the org.Hs.eg.db symbols from the list made above
geneNames<-names(Gene2GO)

#Genes = unique gene symbols in Malderelli data set

#geneList = a numeric representation of which symbols are both in my Malderelli set ("Genes")
#AND in the list I got from org.Hs.eg.db. 0 = not in my list,
# 1 = in my list and the org.Hs.eg.db list

geneList<-factor(as.integer(geneNames %in% Genes))

#add symbols to the numbers so we know which genes are being
#represented
names(geneList)<-geneNames

#allGenes = gene universe = a named (with symbols) list of
#0's and 1's showing which genes in the or.Hs.eg.db list were
# found in Malderelli
GOdata<-new("topGOdata", ontology = "BP", allGenes = geneList,
            annot = annFUN.gene2GO, gene2GO = Gene2GO, nodeSize=5)


#test to see if GO terms showed up in the Malderelli data in the 
#same proportions at they do in the org.Hs.eg.db list??

resultClassic = runTest(GOdata,algorithm = 'classic',
                        statistic =  'fisher')


print(resultClassic)
# 872 significant genes,535 terms p<0.01
#535 terms had showed up in proportions different from in
#org.Hs.eg. db??

resultelim = runTest(GOdata,algorithm = 'elim',
                     statistic =  'fisher')

print(resultelim)# 872significant genes, 122 terms p <0.01

results.table = GenTable(GOdata, Rclassic = resultClassic,
                         Relim = resultelim,
                         topNodes = length(resultClassic@score), 
                         ranksOf = "Relim", orderBy = "Rclassic")

#selecting just rows with GO terms that were significant in the Rclassic test
results.table.bh <- results.table[as.numeric(results.table$Rclassic) < 0.01,]

#just the terms
sigTerms<-results.table.bh$GO.ID


# Get the gene symbols that belong to the significant annotations
#as defined by the classic test. I did this instead
#of the for loop in Kavita's code (it is built into topGo)
genesInGOterms<-genesInTerm(GOdata,sigTerms)

# Extracting all the gene symbols that got annotated by TopGo
#(Called signif in topGO, since I used a pre-defined list,
#all genes that got annotated are considered signif)

signifGenes<-sigGenes(GOdata)

#Make a list showing numerically: of all the gene symbols that
#got annotated, which of them  belong to significant annotations. 

overlapList<-lapply(genesInGOterms,
                    FUN=function(i)(factor(as.integer(signifGenes %in% i))))

#make the list into a dataframe
df<-data.frame(overlapList)

#add the genes list as a column in the df
df<-cbind(signifGenes,df)
names(df)[1]<-"Gene"

#merge in the id columns from Malderelli data by Gene.
#Mald data has genes repeated a lot of times since there are different
#reads for the same gene
df2<-merge(MalderelliData.formatted.likeSCRIData [,1:8], df, by = "Gene")



################################## TEST ###############################
#Check to make sure that merging didn't mess things up

#melt the pre merge df and arrange by gene and go term
melteddf<-melt(df, id.vars = "Gene", variable.name="go_term")
melteddf<-melteddf%>%
  arrange(Gene, go_term)

#filter for just one gene and the go terms it is in
testdf<-melteddf %>%
  filter(Gene == "GNB1", value==1)


x<-unique(testdf$go_term)#GNB1 is a part of 78 go terms

# repeat filtering for the merged df
melteddf2<-melt(df2, id.vars = colnames(df2[1:8]),
                variable.name= "go_term")
melteddf2<-melteddf2%>%
  arrange(Gene, go_term)

testdf2 <- melteddf2 %>%
  filter(Gene == "GNB1", value==1)# there are 234 obs but...

y <- unique(testdf2$go_term) #only 78 unique ones

identical(x,y)# and merging didn't change terms that GNB1 belongs to
####################################################################


#topGO was able to annotate 845 genes but not all were significant
#so there are still some genes that show up in the list, but do
#not belong to any GO term.(I guess? see below)

totalAnnotations<-melteddf2 %>%
  dplyr::select(Gene, go_term, value) %>%
  group_by(Gene) %>%
  summarise(total=sum(as.integer((value))))%>%
  arrange(total)


m<-as.character(totalAnnotations[1:20,"Gene"])

n<-checkGeneSymbols(m,unmapped.as.na = FALSE)


#just want genes that were annotated to at least one term
GenesToKeep<- totalAnnotations %>%
  filter(total != 0)%>%
  dplyr::select(Gene)

#filter df2 for  just the genes to keep found above
df3<-df2[df2$Gene%in% GenesToKeep$Gene,]

#put the columns in the correct order
CLMalderelliData.formatted.likeSCRIData<-df3 %>%
  dplyr::select(Pt., Read, Gene, 4:541)

############################# SAVE ###############################

save(CLMalderelliData.formatted.likeSCRIData,
     file = "CLMalderelliData.formatted.likeSCRIData.Rda")







