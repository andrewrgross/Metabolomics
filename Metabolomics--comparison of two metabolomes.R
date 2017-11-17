### Comparison of two metabolomes -- Andrew R Gross -- 2017/11/10
### A script for reading in two metabolomes and associated data and comparing shared features

########################################################################
### Header

########################################################################
### Import data

setwd('Z:/Data/Andrew/E068-Metabolomics/')
met.neurons <- read.csv('metabolome--iPSC_neurons.csv')
met.plasma <-  read.csv('metabolome--patient_plasma.csv')

########################################################################
### Determine distribution of each metabolome

summary(met.neurons)
