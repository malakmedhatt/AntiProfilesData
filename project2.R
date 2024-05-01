# downloading and accessing the data
install.packages("BiocManager")
library(BiocManager)
BiocManager::install("antiProfilesData")
data <- antiProfilesData::apColonData

#1.1
class(data)
pdata = pData(data)
edata = exprs(data)
fdata = fData(data)

#1.2
pdata_Classes <- sapply(pdata,class)
pdata_Classes

#1.3
ColumnNames <- colnames(pdata)
ColumnNames
RNames <- rownames(pdata)
RNames

#1.4
summary(edata)
summary(pdata)

#1.5
table(pdata$ExperimentID, useNA="ifany")
table(pdata$Tissue, useNA="ifany")
table(pdata$SubType, useNA="ifany")
table(pdata$ClinicalGroup, useNA="ifany")
table(pdata$Status, useNA="ifany")
sapply(pdata[,c("ExperimentID","Tissue","SubType","ClinicalGroup","Status")],table,useNA="ifany" )

#1.6
corr <- cor(edata[,1:10])
covv <- cov(edata[,1:10])
corr
covv
heatmap(corr)

#1.7
#The results show that they have a positive linear relation (as one increase the other increase too)
x <- edata[,'GSM95478']
y <- edata[,'GSM95473']
cor(x, y)
plot(x, y, xlab='GSM95478', ylab='GSM95473', col = c('black', 'blue'))
fit1 <- lm(x ~ y)
abline(fit1, col="red")
equation = paste("y =", round(fit1$coefficients[1],3), "+", round(fit1$coefficients[2],3), "x")
text(4, 20, equation)


#1.8
library(preprocessCore)
norm_edata = normalize.quantiles(as.matrix(edata))

#1.9
normalized <- t(t(norm_edata) - colMeans(norm_edata))
pca <- prcomp(normalized)
eigenVec <- pca$rotation[, 1:2]  #rotation is the V matrix that contains the eigenvectors

plot(eigenVec[, 1], eigenVec[, 2],col= ifelse(pdata$Status == 1, "red", "green"))

#1.10
dist1 = dist(t(edata[,1:10]))
hclust1 = hclust(dist1)
plot(hclust1,hang = -1)

#1.11
tedata <- t(edata)
kmeansData <- kmeans(tedata, centers = 3)
table(kmeansData$cluster)

#2.1
zodiacs <- as.factor(c(rep("aries" ,29), rep("taurus",24), rep("gemini",22), rep("cancer",19), rep("leo",21), rep("virgo",18), rep("libra",19), rep("scorpio",20), rep("sagittarius",23), rep("capricorn",18), rep("Aquarius",20), rep("pisces",23)))
freq = table(zodiacs) 
p <- rep(1/length(freq),length(freq)) #assuming they are evenly distributed that means they all have the same frequency (the sum should be 1)
chisq.test(freq, p=p) #the chi test needs the observed and the expected to be the same length, so we repeated the expected 12 times,

#2.2
#The null hypothesis (H0) is that the zodiac signs are evenly distributed across visual artists (the frequencies are equal)
#in this case the p-val = 0.763 > 0.05 which means we are accepting the null hypothesis that they are evenly distributed
#and rejecting The alternative hypothesis (H1) that states that there is a significant difference between the zodiac signs among the visual artists

#2.3
#Assuming that H0 is that there is no significant difference between the expected and observed frequencies (the observed frequencies represent the true distribution of zodiacs)
p2 <- freq/sum(freq)
chisq.test(freq, p=p2)
freq1=c(29,24,22,19,21,18,19,20,23,18,20,23)
p3=freq1/sum(freq1)
chisq.test(freq, p=p3)
freq


