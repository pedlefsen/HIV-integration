
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

#make a named list of all GO ids in org.Hs.eg.db and the symbols that belong to it
Gene2GO <-inverseList(annFUN.org("BP", mapping = "org.Hs.eg.db",
                                 ID = "symbol"))

geneNames<-names(Gene2GO)

geneList<-factor(as.integer(geneNames %in% Genes))

names(geneList)<-geneNames


GOdata<-new("topGOdata", ontology = "BP", allGenes = geneList,
            annot = annFUN.gene2GO, gene2GO = Gene2GO, nodeSize=5)

resultClassic = runTest(GOdata,algorithm = 'classic',
                        statistic =  'fisher')


print(resultClassic)
# 872 significant genes,535 terms p<0.01

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


#named list of go terms and the genes that go in them. I did this instead
#of the for loop in Kavita's code
genesInGOterms<-genesInTerm(GOdata,sigTerms)

#annotated genes (Called signif in topGO, since I used a pre-defined list,
#all genes that got annotated are considered signif)
signifGenes<-sigGenes(GOdata)

#get logical
overlapList<-lapply(genesInGOterms,
                    FUN=function(i)(factor(as.integer(signifGenes %in% i))))


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







