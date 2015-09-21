##################Getting new go terms for Malderelli data ################
## Writen by Claire Levy
## Modified by Paul E on September 9, 2015

require(biomaRt)
require(dplyr)
require(stringr)
require(HGNChelper)
require(reshape2)

#load data spreadsheet
load( "MalderelliData.formatted.likeSCRIData.Rda" );

######
#using HGNChelper to check gene symbols and provide corrections
######
hgncCheck <- checkGeneSymbols( MalderelliData.formatted.likeSCRIData$Gene,
                              unmapped.as.na = TRUE );
#if hgnc help can't find a good symbol, it can leave as NA or leave
#the input if unmapped.as.na=FALSE

## These didn't map:
# > hgncCheck[ is.na( hgncCheck[ , "Suggested.Symbol" ] ), ]
#                  x Approved Suggested.Symbol
# 209      LOC283050    FALSE             <NA>
# 210      LOC283050    FALSE             <NA>
# 316      LOC643733    FALSE             <NA>
# 324   LOC100526771    FALSE             <NA>
# 325   LOC100526771    FALSE             <NA>
# 578        HDGFRP3    FALSE             <NA>
# 949   LOC100294362    FALSE             <NA>
# 1099  LOC100506012    FALSE             <NA>
# 1257  LOC100505746    FALSE             <NA>
# 1620  LOC100294145    FALSE             <NA>
# 1697     LOC729852    FALSE             <NA>
# 1698     LOC729852    FALSE             <NA>
# 1774 GIMAP1-GIMAP5    FALSE             <NA>
# 1775 GIMAP1-GIMAP5    FALSE             <NA>
# 1777        SGK223    FALSE             <NA>
# 1778        SGK223    FALSE             <NA>

#keep the correct names
MalderelliData.formatted.likeSCRIData.fixed <- MalderelliData.formatted.likeSCRIData;
MalderelliData.formatted.likeSCRIData.fixed$Gene <- hgncCheck$Suggested.Symbol;

# fix this weird one
MalderelliData.formatted.likeSCRIData.fixed$Gene[ grep( "///", MalderelliData.formatted.likeSCRIData.fixed$Gene ) ] <- "MARCH2";

#new gene list for mart has all correct symbols
unique.genes <- unique( MalderelliData.formatted.likeSCRIData.fixed$Gene );

####
# Mapping GO terms
####
#here are the terms that I want to check against the annotations
#that I get for the Malderelli list
terms.we.care.about <-
    grep( "^GO", colnames( MalderelliData.formatted.likeSCRIData ), value = T );

#set up mart

load( "ensemblMart11May15.Rda" );
ensembl<-useMart("ensembl",dataset="hsapiens_gene_ensembl")
# 
# #get entrez id, go id, strand, go name and go category
# unique.genes.GO.annotations <-
#     getBM( attributes=c(
#                "hgnc_symbol",
#                "go_id",
#                "namespace_1003"),
#           mart = ensembl,
#           values = unique.genes,
#           filters="hgnc_symbol" );
# save( unique.genes.GO.annotations, file = "MalderelliAnnotations9Sept15.Rda" );
load( file = "MalderelliAnnotations9Sept15.Rda" );

# BP only, and take out the extra column
unique.genes.GO.annotations.BP <-
    unique.genes.GO.annotations[
      unique.genes.GO.annotations$namespace == "biological_process", 1:2
    ];

.go.terms.columns <- grep( "^GO", colnames( MalderelliData.formatted.likeSCRIData.fixed.further ) );
MalderelliData.formatted.likeSCRIData.fixed.further <-
    apply( MalderelliData.formatted.likeSCRIData.fixed, 1, function( .row ) {
        .row[ .go.terms.columns ] <- 0;
        .go.terms.for.row <- unique.genes.GO.annotations.BP[ unique.genes.GO.annotations.BP[ , "hgnc_symbol" ] == .row[ "Gene" ], "go_id" ];
        # Replace ":" with "." for this.
        .go.terms.for.row <- gsub( ":", ".", .go.terms.for.row );
        if( length( .go.terms.for.row ) == 0 ) {
            ## This means we don't know that they are or are not part
            ## of a term, but since these terms are collections of
            ## what we know and omit many things we do not know, we
            ## let 0 code for "don't know" and leave things as they
            ## are.
        } else {
            # Further filter it.
            .go.terms.for.row.among.those.we.care.about <-
                setdiff( terms.we.care.about, .go.terms.for.row );
            if( length( .go.terms.for.row.among.those.we.care.about ) > 0 ) {
                stopifnot( all( .go.terms.for.row.among.those.we.care.about %in% names( .row ) ) );
                .row[ .go.terms.for.row.among.those.we.care.about ] <- 1;
            }
        }
        return( .row );
    } );
MalderelliData.formatted.likeSCRIData.fixed <-
    t( MalderelliData.formatted.likeSCRIData.fixed.further );
save( MalderelliData.formatted.likeSCRIData.fixed, file = "MalderelliData.formatted.likeSCRIData.fixed.Rda" );
