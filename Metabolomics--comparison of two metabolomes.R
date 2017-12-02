### Comparison of two metabolomes -- Andrew R Gross -- 2017/11/10
### A script for reading in two metabolomes and associated data and comparing shared features

########################################################################
### Header

########################################################################
### Functions

log.transform.and.rm.na <- function(data.frame) {
  data.frame <- log(data.frame + 1)
  data.frame[is.na(data.frame)] <- 1
  return(data.frame)
}

rm.empty.rows <- function(data.frame) {
  maxes <- apply(data.frame,1, max)
  data.frame <- data.frame[maxes>1,]
  return(data.frame)
}

split.groups <- function(data.frame) {     # An intial metabolome is split into control and variable data frames
  labels.vector <- data.frame[1,]
  labels.vector.no.na <- labels.vector[!is.na(labels.vector)]
  groups <- levels(as.factor(labels.vector.no.na))
  
  data.frame <- data.frame[-1,]       # Remove first row and convert character strings to numbers
  indx <- sapply(data.frame, is.factor)
  data.frame[indx] <- lapply(data.frame[indx], function(x) as.numeric(as.character(x)))
  
  data.frame <- log.transform.and.rm.na(data.frame)
  data.frame <- rm.empty.rows(data.frame)
  
  group1.df <- data.frame[which(labels.vector==groups[1])]
  group2.df <- data.frame[which(labels.vector==groups[2])]
  output.list <- list(group1.df,group2.df)
  names(output.list) <- groups
  return(output.list)
}   
  
t.test.between.groups <- function(input.list, control.name = 'CTR'){
  groups <- names(input.list)
  ctrl.df <- input.list[[which(groups==control.name)]]
  var.df <- input.list[[which(groups!=control.name)]]
  
  output.df <- ctrl.df[1:2]
  ### Calculate means
  output.df[1] <- apply(ctrl.df,1,mean,na.rm = TRUE)
  output.df[2] <- apply(var.df, 1,mean,na.rm = TRUE)
  names(output.df) <- c('CTR', 'ALT')
  ### Calculate p.val
  p.val.vector <- c()
  for(row.num in 1:nrow(ctrl.df)) {
    p.val <- try(p.val <- t.test(ctrl.df[row.num,], var.df[row.num,])[[3]],silent = FALSE)
    if(is.numeric(p.val) == FALSE ) { p.val = 1}
    p.val.vector <- c(p.val.vector, p.val)
  }
  output.df <- cbind(output.df, 'p.val' = p.val.vector)
  
  output.df <- output.df[order(output.df$p.val), ]
}



########################################################################
### Import data

setwd('Z:/Data/Andrew/E068-Metabolomics/')
met.neurons <- read.csv('metabolome--iPSC_neurons.csv',row.names = 1)
met.plasma <-  read.csv('metabolome--patient_plasma.csv', row.names = 1)

fc_list_neurons <- read.csv('Temp files/neuron_fc_list.csv')
fc_list_plasma <- read.csv('Temp files/plasma_fc_list.csv')

########################################################################
### Determine distribution of each metabolome

summary(met.neurons)
summary(met.plasma)

output.neurons <- split.groups(met.neurons)
input.list <- output.neurons 
comparison.neurons <- t.test.between.groups(output.neurons)

head(comparison.neurons,20)

output.plasma <- split.groups(met.plasma)

comparison.plasma <- t.test.between.groups(output.plasma)

### Filter

test <- comparison.neurons[comparison.neurons$p.val < 0.05,]

comparison.neurons2 <- comparison.neurons[comparison.neurons$p.val < 0.1,]

comparison.plasma2 <- comparison.plasma[comparison.plasma$p.val < 0.1,]

intersect(row.names(comparison.plasma2), row.names(comparison.neurons2))

