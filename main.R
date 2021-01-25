library( tidyverse )

## library( gelnet )     # https://github.com/ArtemSokolov/gelnet

## Load RNAseq data
X <- read.delim("./data/PCBC/rnaseq_norm.tsv") %>%
    tibble::column_to_rownames( "tracking_id" ) %>%
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

