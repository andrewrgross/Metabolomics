### Metabolomics compound converter -- Andrew R Gross -- 2017/12/01
### A script for converting compound names

########################################################################
### Header

########################################################################
### Functions

convert.IDs <- function(dataframe) {
  for(row.num in 1:nrow(dataframe)) {       # Find matching name for each row
    old.ID <- as.character(dataframe[row.num,][,1])       # Define the old name
    old.ID <- tolower(paste0('^',old.ID,'$'))
    new.ID <- as.character(name.map$FINAL.NAME[grep(old.ID, name.map$ORIGINAL.NAME)]) # Define the new name
    dataframe[row.num,][,1] <- new.ID
  }
  return(dataframe)
}

convert.to.numeric <-function(dataframe) {
  for(col.num in 2:ncol(dataframe)) {
    column <- as.numeric(dataframe[,col.num])
    dataframe[col.num] <- column
  }
  return(dataframe)
}

sum.duplicate.rows <- function(dataframe) {
  full.ids <- tolower(dataframe[,1])
  unique.ids <- unique(full.ids)
  unique.ids <- unique.ids[!is.na(unique.ids)]
  new.matrix <- matrix(1, length(unique.ids), ncol(dataframe)-1)
  new.df <- data.frame(new.matrix)
  names(new.df) <- names(dataframe[-1])
  row.names(new.df) <- unique.ids
  
  for(row.num in 1:nrow(new.df)) {
    metabolite <- unique.ids[row.num]
    metabolite <- paste0('^',metabolite,'$')
    positions <- grep(metabolite, tolower(dataframe[,1]))
    temp.df <- dataframe[-1][positions,]
    if(nrow(temp.df) >1) {
      #print(metabolite)
      #print(temp.df)
      #print(new.row)
    }
    new.row <- apply(temp.df,2,sum)
    new.df[row.num,] <- new.row
  }
  return(new.df)
}

rank.by.mean <- function(dataframe) {
  row.means <- apply(dataframe,1, mean)
  dataframe <- cbind(dataframe,row.means)
  dataframe <- dataframe[order(dataframe$row.means,decreasing= TRUE),]
  dataframe <- dataframe[ncol(dataframe)]
  return(dataframe)
}

########################################################################
### Import data

setwd('Z:/Data/Andrew/E068-Metabolomics/')
name.map <- read.csv('name_map_MAIN.csv')[1:2]
name.map$ORIGINAL.NAME <- tolower(name.map$ORIGINAL.NAME)
met.neurons <- read.csv('metabolome--iPSC_neurons.csv', stringsAsFactors = FALSE)
met.plasma <-  read.csv('metabolome--patient_plasma.csv', stringsAsFactors = FALSE)

fc_list_neurons <- read.csv('Temp files/neuron_fc_list.csv')
fc_list_plasma <- read.csv('Temp files/plasma_fc_list.csv')

### Convert NA to 0

met.neurons[is.na(met.neurons)] <- 0

########################################################################
### Sort

met.neurons.ctr <- met.neurons[c(1,4,6,7,9,11,15)][-1,]
met.neurons.als <- met.neurons[c(1,3,5,8,10,12,13,14)][-1,]

met.plasma.ctr <- met.plasma[c(1,12:21)][-1,]
met.plasma.als <- met.plasma[1:11][-1,]

########################################################################
### Convert IDs

met.neurons <- convert.IDs(met.neurons)
met.neurons.ctr <- convert.IDs(met.neurons.ctr)
met.neurons.als <- convert.IDs(met.neurons.als)

met.plasma <- convert.IDs(met.plasma)
met.plasma.ctr <- convert.IDs(met.plasma.ctr)
met.plasma.als <- convert.IDs(met.plasma.als)

### Convert to numeric

met.neurons <- convert.to.numeric(met.neurons)
met.neurons.ctr <- convert.to.numeric(met.neurons.ctr)
met.neurons.als <- convert.to.numeric(met.neurons.als)

met.plasma <- convert.to.numeric(met.plasma)
met.plasma.ctr <- convert.to.numeric(met.plasma.ctr)
met.plasma.als <- convert.to.numeric(met.plasma.als)

### Drop NA rows

met.neurons <- met.neurons[!is.na(met.neurons[,1]),]
met.neurons.ctr <- met.neurons.ctr[!is.na(met.neurons.ctr[,1]),]
met.neurons.als <- met.neurons.als[!is.na(met.neurons.als[,1]),]

met.plasma <- met.plasma[!is.na(met.plasma[,1]),]
met.plasma.ctr <- met.plasma.ctr[!is.na(met.plasma.ctr[,1]),]
met.plasma.als <- met.plasma.als[!is.na(met.plasma.als[,1]),]

########################################################################
### Sum duplicate rows

met.neurons.summed <- sum.duplicate.rows(met.neurons)
met.neurons.ctr.summed <- sum.duplicate.rows(met.neurons.ctr)
met.neurons.als.summed <- sum.duplicate.rows(met.neurons.als)

met.plasma.summed <- sum.duplicate.rows(met.plasma)
met.plasma.ctr.summed <- sum.duplicate.rows(met.plasma.ctr)
met.plasma.als.summed <- sum.duplicate.rows(met.plasma.als)

########################################################################
### Rank by mean

met.neurons.ranked <- rank.by.mean(met.neurons.summed)
met.neurons.ctr.ranked <- rank.by.mean(met.neurons.ctr.summed)
met.neurons.als.ranked <- rank.by.mean(met.neurons.als.summed)

met.plasma.ranked <- rank.by.mean(met.plasma.summed)
met.plasma.ctr.ranked <- rank.by.mean(met.plasma.ctr.summed)
met.plasma.als.ranked <- rank.by.mean(met.plasma.als.summed)

########################################################################
### Filter 

#met.neurons.all.ranked <- cbind(met.neurons, s)

########################################################################
### write

write.csv(met.neurons.ranked,'met.neurons.ranked.csv')

write.csv(met.plasma.ranked,'met.plasma.ranked.csv')

########################################################################
### scratchwork
met.neurons[,1]
row.names(met.neurons[1,])
new.df <- met.plasma
for(row.num in 1:nrow(new.df)) {       # Find matching name for each row
  old.ID <- as.character(new.df[row.num,][,1])       # Define the old name
  old.ID <- tolower(paste0('^',old.ID,'$'))
  new.ID <- as.character(name.map$FINAL.NAME[grep(old.ID, name.map$ORIGINAL.NAME)]) # Define the new name
  new.df[row.num,][,1] <- new.ID
}

row.num <- 2
met.neurons[row.num,][,1]

(old.ID <- as.character(met.neurons[row.num,][,1]))

grep(tolower(old.ID), name.map$ORIGINAL.NAME)
as.character(name.map$FINAL.NAME[grep(tolower(old.ID), name.map$ORIGINAL.NAME)])


grep(tolower(test), name.map$ORIGINAL.NAME)

test <- paste0('^',old.ID,'$')

df1 <- head(numeric.dataframe)
df1[] <- lapply(df1, function(x) {
  if(is.factor(x)) as.numeric(as.character(x)) else x
})
sapply(df1, class)
