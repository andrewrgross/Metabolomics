### Getting the shape of metabolomics data

### Header
library(ggplot2)

### Functions

generate.stats <- function(dataframe) {
  stats <- apply(dataframe,2,quantile)
  stats <- t(as.matrix(stats))
  stats <- as.data.frame(stats)
  stats$line <- row.names(stats)
  names(stats) <- c('ymin', 'lower', 'middle', 'upper', 'ymax', 'line')
  return(stats)
}

plot.boxplot <- function(stats) {
  plot <- ggplot(data = stats, aes(x = line)) +
    geom_boxplot(aes(ymin = ymin, lower = lower, middle = middle, upper = upper, ymax = ymax), stat = 'identity') +
    coord_flip()
  return(plot)
}

plot.density.within.range <- function(array,lower.bound,upper.bound) {
  data <- as.matrix(array)
  data <- data[data > lower.bound]
  data <- data[data < upper.bound]
  plot(density(data))
}

plot.density.range.adjust <- function(array,lower.bound,upper.bound, adjust = 1) {
  data <- as.matrix(array)
  data <- data[data > lower.bound]
  data <- data[data < upper.bound]
  plot(density(data,adjust))
}

### Data Input
metabo.df <- read.csv('Z:/Data/Andrew/Metabolomics/031517QvdardovMets1_18 analyzed_TEMP_no-blank.csv', stringsAsFactors = FALSE)
samples.df <- metabo.df

### Data filtering
head(metabo.df)
metabo.df <- metabo.df[-1,]
metabo.df <- metabo.df[-1]
metabo.df[is.na(metabo.df)] <- 1
metabo.df[1:ncol(metabo.df)] <- as.numeric(as.matrix(metabo.df[1:ncol(metabo.df)]))
summary(metabo.df)

### Format to plot
stats <- apply(metabo.df,2,quantile)
stats <- t(as.matrix(stats))
stats <- as.data.frame(stats)
stats$line <- row.names(stats)
names(stats) <- c('ymin', 'lower', 'middle', 'upper', 'ymax', 'line')

### Plot
# Boxplots
ggplot(data = stats, aes(x = line)) +
  geom_boxplot(aes(ymin = ymin, lower = lower, middle = middle, upper = upper, ymax = ymax), stat = 'identity') +
  coord_flip()

# Density plot of concentration values
plot(density(t(as.matrix(metabo.df))))

# Filtered density plot of concentrations
metabo.filtered <- metabo.df[metabo.df>1.4e7]
metabo.filtered2 <- metabo.filtered[metabo.filtered <1.9e7]
plot(density(metabo.filtered2))

# Filtered density plot of concentrations
metabo.filtered <- metabo.df[metabo.df>0]
metabo.filtered2 <- metabo.filtered[metabo.filtered <1.9e5]
plot(density(metabo.filtered2))


### Transforming
log.metabo.df <- log(metabo.df,2)
summary(log.metabo.df)

### Format to plot
stats <- generate.stats(log.metabo.df)

### Plot
# Boxplots
(g <- plot.boxplot(stats))

# Density plot of concentration values
plot(density(as.matrix(log.metabo.df)))

# Filtered density plot of concentrations
plot.density.within.range(log.metabo.df,12,30)
plot.density.within.range(log.metabo.df,16.7,17.4)



### Transforming -- Cube root
cr.metabo.df <- (metabo.df)^(1/3)
summary(cr.metabo.df)

### Format to plot
stats <- generate.stats(cr.metabo.df)

### Plot
# Boxplots
(g <- plot.boxplot(stats))

# Density plot of concentration values
plot(density(as.matrix(cr.metabo.df)))

# Filtered density plot of concentrations
plot.density.within.range(log.metabo.df,95,110)



### Examine means of both groups

means.df <- data.frame('ctr.mean' = apply(metabo.df[c(1,3,4,6,8,11)],1,mean),
                       'als.mean' = apply(metabo.df[c(2,5,7,9,10)],1,mean),
                       'ALL.mean' = apply(metabo.df,1,mean))

means2.df <- means.df[rev(order(means.df$ALL.mean)),]

length(means.df$ALL.mean[means.df$ALL.mean<2])

empty.rows <- row.names(means.df[means.df$ALL.mean <2,])

bad.metabo.df <- metabo.df[empty.rows,]

good.rows <- row.names(means.df[means.df$ALL.mean >500,])
good.metabo.df <- metabo.df[row.names(means.df[means2.df$ALL.mean >10,]),]
summary(means2.df[good.rows,])



### Examining pruned set
metabo.df <- good.metabo.df
### Format to plot
stats <- generate.stats(metabo.df)

### Plot
# Boxplots
ggplot(data = stats, aes(x = line)) +
  geom_boxplot(aes(ymin = ymin, lower = lower, middle = middle, upper = upper, ymax = ymax), stat = 'identity') +
  coord_flip()

# Density plot of concentration values
plot(density(t(as.matrix(metabo.df))))

# Filtered density plot of concentrations
metabo.filtered <- metabo.df[metabo.df>1.4e7]
metabo.filtered2 <- metabo.filtered[metabo.filtered <1.9e7]
plot(density(metabo.filtered2))

# Filtered density plot of concentrations
metabo.filtered <- metabo.df[metabo.df>0]
metabo.filtered2 <- metabo.filtered[metabo.filtered <1.9e5]
plot(density(metabo.filtered2))


### Transforming
log.metabo.df <- log(metabo.df,2)
summary(log.metabo.df)

### Format to plot
stats <- apply(log.metabo.df,2,quantile)
stats <- t(as.matrix(stats))
stats <- as.data.frame(stats)
stats$line <- row.names(stats)
names(stats) <- c('ymin', 'lower', 'middle', 'upper', 'ymax', 'line')

### Plot
# Boxplots
ggplot(data = stats, aes(x = line)) +
  geom_boxplot(aes(ymin = ymin, lower = lower, middle = middle, upper = upper, ymax = ymax), stat = 'identity') +
  coord_flip()

# Density plot of concentration values
plot(density(as.matrix(log.metabo.df)))

# Filtered density plot of concentrations
metabo.filtered <- log.metabo.df[log.metabo.df>12]
metabo.filtered2 <- metabo.filtered[metabo.filtered <30]
plot(density(metabo.filtered2))


metabo.filtered <- log.metabo.df[log.metabo.df>16.7]
metabo.filtered2 <- metabo.filtered[metabo.filtered <17.4]
plot(density(metabo.filtered2))




### Transforming -- Cube root
cr.metabo.df <- (metabo.df)^(1/3)
summary(cr.metabo.df)

### Format to plot
stats <- apply(cr.metabo.df,2,quantile)
stats <- t(as.matrix(stats))
stats <- as.data.frame(stats)
stats$line <- row.names(stats)
names(stats) <- c('ymin', 'lower', 'middle', 'upper', 'ymax', 'line')

### Plot
# Boxplots
ggplot(data = stats, aes(x = line)) +
  geom_boxplot(aes(ymin = ymin, lower = lower, middle = middle, upper = upper, ymax = ymax), stat = 'identity') +
  coord_flip()

# Density plot of concentration values
plot(density(as.matrix(cr.metabo.df)))

# Filtered density plot of concentrations
metabo.filtered <- cr.metabo.df[cr.metabo.df>95]
metabo.filtered2 <- metabo.filtered[metabo.filtered <110]
plot(density(metabo.filtered2))


metabo.filtered <- log.metabo.df[log.metabo.df>16.7]
metabo.filtered2 <- metabo.filtered[metabo.filtered <17.4]
plot(density(metabo.filtered2))




####################################
### Generate one df with no zero rows

non.zero.rows <- apply(metabo.df,1,min) > 4
non.zero.metabo.df <- metabo.df[non.zero.rows,]
summary(non.zero.metabo.df)
nrow(non.zero.metabo.df)

stats <- generate.stats(non.zero.metabo.df)
plot.boxplot(stats)

plot(density(as.matrix(non.zero.metabo.df)))

plot.density.within.range(non.zero.metabo.df,0,1e7)

### Transform after filtering
log.non.zero.metabo.df <- log(non.zero.metabo.df,2)
stats <- generate.stats(log.non.zero.metabo.df)
plot.boxplot(stats)

plot(density(as.matrix(log.non.zero.metabo.df)))

plot.density.within.range(log.non.zero.metabo.df,0,1e7)

hist(log.non.zero.metabo.df)

### Cube root transform
cr.non.zero.metabo.df <- (non.zero.metabo.df)^(1/3)

summary(cr.non.zero.metabo.df)

stats <- generate.stats(cr.non.zero.metabo.df)
plot.boxplot(stats)

plot(density(as.matrix(cr.non.zero.metabo.df)))

plot.density.within.range(cr.non.zero.metabo.df,0,1e3)

hist(log.non.zero.metabo.df)

### Normalize by median

norm.metabo.df <- log.non.zero.metabo.df
summary(log.non.zero.metabo.df)

subtract.median <- function(dataframe) {
  for(column in 1:ncol(dataframe)) {
    dataframe[column] <- dataframe[,column] - median(dataframe[,column])
  }
  return(dataframe)
}

norm.metabo.df <- subtract.median(norm.metabo.df)
stats <- generate.stats(norm.metabo.df)
plot.boxplot(stats)

plot(density(as.matrix(norm.metabo.df)))

hist(norm.metabo.df)

plot(density(norm.metabo.df[,4]))


#####################################################
### Looking for bimodality in the original data

### Data Input
metabo.df <- read.csv('Z:/Data/Andrew/Metabolomics/031517QvdardovMets1_18 analyzed_TEMP_no-blank.csv', stringsAsFactors = FALSE)
df.fullwood <- read.table("Z:/Data/Andrew/reference_data/geo/GSE69360/GSE69360_RNAseq.counts.txt", sep = '\t', header = TRUE, row.names = 1)[6:24]
df.yu <- read.csv('Z:/Data/Andrew/reference_data/geo/yu_GSE9440.csv', row.names = 1)
df.han <- read.csv('Z:/Data/Andrew/reference_data/geo/han_GSE35108.csv', row.names = 1)

remove.rows.by.median <- function(dataframe,cutoff){
  dataframe.medians <- apply(dataframe,1,median)
  rows.to.keep <- dataframe.medians > cutoff
  return(dataframe[rows.to.keep,])
}

summary(df.fullwood)
fullwood.filtered <- as.matrix(remove.rows.by.median(df.fullwood, 1))
summary(fullwood.filtered)

plot(density(log(as.matrix(df.fullwood)),adjust = 0.1))
plot(density(log2(fullwood.filtered),adjust = 0.1))

summary(df.yu)
yu.filtered <- as.matrix(remove.rows.by.median(df.yu, 0.75))
plot(density(log2(10*yu.filtered),adjust=0.1))

summary(df.han)
plot(density(as.matrix(df.han),adjust = 0.01))

summary(metabo.df)
plot(density(log(as.matrix(metabo.df),2), adjust = 1.0))






# Import, format, filter rows in the lowest 10 percentiles
fullwood.medians <- apply(df.fullwood,1,median)
fullwood.mins <- apply(df.fullwood,1,min)
q <- quantile(fullwood.medians)
m.fullwood <- as.matrix(df.fullwood)

fullwood.max <- apply(df.fullwood,1,max)
non.zero.rows <- fullwood.max > 0
df.fullwood <- df.fullwood[non.zero.rows,]

### Yu
m.yu <- as.matrix(df.yu)

## Han
m.han <- as.matrix(df.han)

TPMdata <- read.csv("z://Data/RNAseq HT neurons and tissue/Andrews_files/20160317_Sareen_rerun-29299281.tpm.csv", row.names=1)
summary(TPMdata)
tpm.median <- apply(TPMdata, 1, median)
tpm.max <- apply(TPMdata, 1, max)
non.zero.rows <- tpm.median > 1
df.tpm <- TPMdata[non.zero.rows,]
summary(df.tpm)
m.tpm <- as.matrix(df.tpm)




### Data filtering
head(metabo.df)
metabo.df <- metabo.df[-1,]
metabo.df <- metabo.df[-1]
metabo.df <- na.omit(metabo.df)
metabo.df[1:ncol(metabo.df)] <- as.numeric(as.matrix(metabo.df[1:ncol(metabo.df)]))
summary(metabo.df)


plot(density(metabo.df[,1]))
plot(density(metabo.df[,1],adjust = 0.1))
plot.density.within.range(metabo.df[,1],0,1e8)
plot.density.within.range(metabo.df[,1],0,1e6)
plot.density.range.adjust(metabo.df[,1],0,1e6,100)
plot.density.range.adjust(metabo.df,0,1e6,100)
plot.density.within.range(metabo.df[,1],0,1e4)

plot.density.within.range(log(metabo.df[,1]),0,1e8)
plot.density.range.adjust(log(metabo.df[,1]),0,1e8,0.01)
plot.density.range.adjust(log(metabo.df),0,1e8,0.01)
plot.density.range.adjust(log(metabo.df),0,1e8,0.1)
### Generate an actual guasian curve and see if this does the same thing


plot(density(as.matrix(df.fullwood)))
plot.density.within.range(as.matrix(df.fullwood),0,1e3)
plot(density(log(as.matrix(df.fullwood)),adjust = 0.1),1)
plot.density.within.range(log(df.fullwood),0,1e3)

plot(density(m.yu))
plot.density.within.range(m.yu,0,10)
plot(density(log(log(log(m.yu+1)+1)+1),adjust = 5))

plot(density(m.han, adjust = 0.1))

plot(density(m.tpm))
plot(density(log(m.tpm+50),adjust = 0.1))

###################################
### Online example
library(e1071)
data(airquality)
ozone <- airquality$Ozone
ozone <- ozone[!is.na(ozone)]
hist(ozone, col = "tomato")
plot(density(ozone,adjust = 0.1))
skewness(ozone)
skew.score <- function(c, x) (skewness(log(x + c)))^2
cval <- seq(0, 20, l = 101)
skew <- cval * 0
for (i in 1:length(cval)) skew[i] <- skewness(log(cval[i] + ozone))
plot(cval, skew, type = "l", ylab = expression(b[3](c)), xlab = expression(c))
abline(h = 0, lty = 3)
best.c <- optimise(skew.score, c(0, 20), x = ozone)$minimum
best.c
ozone.transformed <- log(ozone + best.c)
hist(ozone.transformed, col = "azure")
plot(density(ozone.transformed,adjust = 1))

qqnorm(ozone)
qqnorm(log(ozone))
qqline(log(ozone))
qqnorm(ozone.transformed)
qqline(ozone.transformed)

### Now with my data
hist(df.fullwood[1])
skewness(df.fullwood[,1])

cval <- seq(0, 20, l = 101)
skew <- cval * 0
for (i in 1:length(cval)) skew[i] <- skewness(log(cval[i] + df.fullwood[,1]))
plot(cval, skew, type = "l", ylab = expression(b[3](c)), xlab = expression(c))
abline(h = 0, lty = 3)
best.c <- optimise(skew.score, c(0, 20), x = df.fullwood[,1])$minimum
best.c

df.f.transformed <- log(df.fullwood[,1] + best.c)
hist(df.f.transformed, col = "azure")
plot(density(df.f.transformed,adjust = 1))

qqnorm(ozone)
qqnorm(log(ozone))
qqline(log(ozone))
qqnorm(df.f.transformed)
qqline(ozone.transformed)
