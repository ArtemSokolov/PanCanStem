#!/usr/bin/env Rscript

library( synapser )

argv <- commandArgs( trailingOnly=TRUE )
if( length(argv) < 2 ) stop("Please provide Synapse credentials")

## All data will be downloaded from Synapse
synLogin(argv[1], argv[2])

## Download PCBC RNAseq data
synGet( "syn2701943", downloadLocation="/data/PCBC",
       ifcollision="overwrite.local" )

## Download PCBC annotations
qry <- synTableQuery( "SELECT UID, Diffname_short FROM syn3156503" )
file.copy( qry$filepath, "./data/PCBC/meta.csv" )

## Download PAnCan33 data
synGet( "syn4976369", downloadLocation="/data/PanCan33",
       ifcollision="overwrite.local" )
