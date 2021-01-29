#!/usr/bin/env Rscript

library( tidyverse )
library( gelnet )     # https://github.com/ArtemSokolov/gelnet

# Maps ENSEMBL IDs to HUGO
# Use srcType = "ensembl_gene_id" for Ensembl IDs
# Use srcType = "entrezgene" for Entrez IDs
genes2hugo <- function( v, srcType = "ensembl_gene_id" )
{
    ## Retrieve the EMSEMBL -> HUGO mapping
    ensembl <- biomaRt::useMart( "ENSEMBL_MART_ENSEMBL", host="www.ensembl.org",
                                dataset="hsapiens_gene_ensembl" )
    ID <- biomaRt::getBM( attributes=c(srcType, "hgnc_symbol"),
                         filters=srcType, values=v, mart=ensembl )

    ## Make sure there was at least one mapping
    if( nrow(ID) < 1 ) top( "No IDs mapped successfully" )
    
    ## Drop empty duds
    j <- which( ID[,2] == "" )
    if( length(j) > 0 ) ID <- ID[-j,]
    stopifnot( all( ID[,1] %in% v ) )
    ID
}

cat( "Loading PCBC data...\n" )

## Load RNAseq data
X <- read.delim("/data/PCBC/rnaseq_norm.tsv") %>%
    mutate( tracking_id = str_split( tracking_id, "\\.", simplify=TRUE )[,1] ) %>%
    column_to_rownames( "tracking_id" ) %>%
    as.matrix()

## Load metadata
## Fill missing labels by hand
Y <- read_csv("/data/PCBC/meta.csv", col_types=cols()) %>%
    select( UID, Class=Diffname_short ) %>%
    mutate( across(UID, ~gsub("-", ".", .x)) ) %>%
    add_row(UID = c("SC11.014BEB.133.5.6.11", "SC12.039ECTO.420.436.92.16"),
            Class = c("EB", "ECTO"))

## Isolate all labels for which we have data
y <- with( Y, set_names(Class, UID) )[colnames(X)]

cat( "Mapping Ensembl IDs to HGNC...\n" )
V <- genes2hugo( rownames(X) )
X <- X[V[,1],]
rownames(X) <- V[,2]

## Mean-center the data
m <- apply( X, 1, mean )
X <- X - m

## Identify stem cell samples
stopifnot( identical(colnames(X), names(y)) )
X.tr <- X[,which( y == "SC" )]

mdef <- gelnet( t(X.tr) ) + model_oclr() + rglz_L2(1)
mdl <- gelnet_train( mdef )

## Compose the signature and write it to a file
stemsig <- tibble(Gene = rownames(X), Weight = mdl$w)
stemsig %>% write_csv( "/data/pcbc-stemsig.csv" )

cat( "Loading PanCan33 data...\n" )
PCraw <- read_tsv( "/data/PanCan33/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv",
                  col_types=cols())

cat( "Cleaing up IDs...\n" )
PC <- PCraw %>% filter( !grepl("\\?", gene_id) ) %>%
    mutate( across(gene_id, ~str_split(.x, "\\|", simplify=TRUE)[,1]) ) %>%
    filter( !duplicated(gene_id) ) %>%
    rename( Gene = gene_id )

## Score each sample via Spearman correlation against the signature
##  to produce stemness indices
cat( "Scoring PanCan33 data...\n" )
mRNAsi <- stemsig %>% inner_join( PC, by="Gene" ) %>%
    summarize( across(c(-Gene, -Weight), cor, Weight,
                      method="sp", use = "complete.obs") ) %>%
    gather( Sample, mRNAsi ) %>%
    mutate( across(mRNAsi, ~(.x - min(.x))/(max(.x) - min(.x))) )

## Write out stemness indices to file
write_csv( mRNAsi, "/data/pancan33-mRNAsi.csv" )
