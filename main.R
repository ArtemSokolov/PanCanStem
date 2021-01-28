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

## Load RNAseq data
X <- read_tsv("./data/PCBC/rnaseq_norm.tsv", col_types=cols()) %>%
    mutate( tracking_id = str_split( tracking_id, "\\.", simplify=TRUE )[,1] ) %>%
    column_to_rownames( "tracking_id" ) %>%
    as.matrix()

## Load metadata
## Fill missing labels by hand
Y <- read_csv("./data/PCBC/meta.csv", col_types=cols()) %>%
    select( UID, Class=Diffname_short ) %>%
    mutate( across(UID, ~gsub("-", ".", .x)) ) %>%
    add_row(UID = c("SC11.014BEB.133.5.6.11", "SC12.039ECTO.420.436.92.16"),
            Class = c("EB", "ECTO"))

## Isolate all labels for which we have data
y <- with( Y, set_names(Class, UID) )[colnames(X)]

## Map Ensembl IDs to HUGO
V <- genes2hugo( rownames(X) )
X <- X[V[,1],]
rownames(X) <- V[,2]

## Mean-center the data
m <- apply( X, 1, mean )
X <- X - m

## Identify stem cell samples
stopifnot( identical(colnames(X), names(y)) )
X.tr <- X[,which( y == "SC" )]

## Train a one-class model
mdef <- gelnet( t(X.tr) ) + model_oclr() + rglz_L2(1)
mdl <- gelnet_train( mdef )

## Compose the signature and write it to a file
stemsig <- tibble(Gene = rownames(X), Weight = mdl$w)
stemsig %>% write_csv( "pcbc-stemsig.csv" )

