### Metabolomic analysis -- Andrew R Gross -- 2017-03-16
### A script for examining metabolomics data

### Header

### Functions

### Input

df.metabolomics <- read.csv("Z:/Data/Andrew/E068-Metabolomics/metabolome--iPSC_neurons.csv",stringsAsFactors = FALSE, na.strings = 'N/A', skip = 2)

df.metabolomics <- df.metabolomics[1:15]

### Format
### Format names
labels <- c('labels', 'Blank_1',	'30iALS_1',	'03iCTR_1',	'134iALS_1',	'14iCTR_1',	'03iCTR_2',	'52iALS_1',	'25iCTR_1',	'52iALS_2',	'25iCTR_2',	'Blank_2',	'134iALS_2',	'30iALS_2',	'14iCTR_2')
names(df.metabolomics) <- labels

### Change order
df.metabolomics2 <- df.metabolomics[c(1,2,12,4,7,6,15,9,11,3,14,8,10,5,13)]

### Replace NA with 0
df.metabolomics2[is.na(df.metabolomics2)] = 0

### Optional log transformation
#df.metabolomics2[2:15] <- log(df.metabolomics2[2:15])        # Didn't seem to help

### Generate avgs of each sample
df.metabolomics.avg <- df.metabolomics2[1]
df.metabolomics.avg$blank <- apply(df.metabolomics2[2:3],1,mean)
df.metabolomics.avg$'03iCTR' <- apply(df.metabolomics2[4:5],1,mean)
df.metabolomics.avg$'14iCTR' <- apply(df.metabolomics2[6:7],1,mean)
df.metabolomics.avg$'25iCTR' <- apply(df.metabolomics2[8:9],1,mean)
df.metabolomics.avg$'30iALS' <- apply(df.metabolomics2[10:11],1,mean)
df.metabolomics.avg$'52iALS' <- apply(df.metabolomics2[12:13],1,mean)
df.metabolomics.avg$'134iALS' <- apply(df.metabolomics2[14:15],1,mean)

### Calculate the variance of each sample
df.metabolomics.var <- df.metabolomics[1]
df.metabolomics.var$blank <- apply(df.metabolomics2[2:3],1,diff)
df.metabolomics.var$'03iCTR' <- apply(df.metabolomics2[4:5],1,diff)
df.metabolomics.var$'14iCTR' <- apply(df.metabolomics2[6:7],1,diff)
df.metabolomics.var$'25iCTR' <- apply(df.metabolomics2[8:9],1,diff)
df.metabolomics.var$'30iALS' <- apply(df.metabolomics2[10:11],1,diff)
df.metabolomics.var$'52iALS' <- apply(df.metabolomics2[12:13],1,diff)
df.metabolomics.var$'134iALS' <- apply(df.metabolomics2[14:15],1,diff)

### Calculate the ratio of the variance of each sample
df.metabolomics.dev <- df.metabolomics.var
df.metabolomics.dev[2:8] <- abs(round(df.metabolomics.var[2:8]/(df.metabolomics.avg[2:8]+0.1)*100))

head(df.metabolomics.dev)
summary(df.metabolomics.dev)




