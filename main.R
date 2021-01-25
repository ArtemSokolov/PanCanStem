library( tidyverse )

## The script also requires the following:
##
## library( synapser )   # https://github.com/Sage-Bionetworks/synapser
## library( synExtra )   # https://github.com/ArtemSokolov/synExtra
## library( gelnet )     # https://github.com/ArtemSokolov/gelnet

## Define the Synapse interface
synapser::synLogin()
syn <- synExtra::synDownloader( "./data", ifcollision="overwrite.local" )

## Load RNAseq data
X <- syn( "syn2701943" ) %>% read.delim() %>%
    tibble::column_to_rownames( "tracking_id" ) %>%
    as.matrix()

## Retrieve metadata
## Fill missing labels by hand
synMeta <- synapser::synTableQuery( "SELECT UID, Diffname_short FROM syn3156503" )
Y <- read_csv(synMeta$filepath, col_types=cols()) %>%
    select( UID, Class=Diffname_short ) %>%
    mutate( across(UID, ~gsub("-", ".", .x)) ) %>%
    add_row(UID = c("SC11.014BEB.133.5.6.11", "SC12.039ECTO.420.436.92.16"),
            Class = c("EB", "ECTO"))

## Isolate all labels for which we have data
y <- with( Y, set_names(Class, UID) )[colnames(X)]
