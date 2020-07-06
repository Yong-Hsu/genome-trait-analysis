library(gtools)
#load sorted files
dataFiles <- sapply(mixedsort(Sys.glob("data/gene_info/gene_*.dat")), read.table )
locus = read.table("snp.dat",header = T)

#create list to store data
genelm = list(list(NULL))
length(genelm) = 300

#loop to store data
for (i in 1:300){
  a = cbind2( c(rep(0,500),rep(1,500)) , locus[ , match(dataFiles[[i]] , names(locus)) ] )

  genelm[[i]] <- lm(x~.,data = a)  #multiple llinear regression
}

#create null vector to store pvalue in F-test of regression
pvalue<-vector(mode="numeric",length=300)
#calculate pvalue
for (i in 1:300){
  pvalue[i] = pf(summary(genelm[[i]])$fstatistic[1L], summary(genelm[[i]])$fstatistic[2L], summary(genelm[[i]])$fstatistic[3L], lower.tail = FALSE)
}
#minor procession
p = -log(pvalue,10)

#draw pvalue barplot
pdata = data.frame(gene = rep(1:300),pvalue = t(p))
library(ggplot2)
ggplot(data=pdata, aes(x=gene, y=p)) +
  geom_bar(stat = "identity", width = 0.07 ,color = "steelblue") +
  theme_minimal() +
  theme(axis.text = element_text(size = 13),axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))+
  geom_hline(aes(yintercept=2), colour="#990000", linetype="dashed")

#calculate r-squared value
adjrsquared = vector(mode="numeric",length=300)
for (i in 1:300){
  adjrsquared[i] = summary(genelm[[i]])$adj.r.squared
}

library(ggplot2)

pdata = data.frame(gene = rep(1:300),adjrsquared = t(adjrsquared))
shift_trans = function(d = 0) {
  scales::trans_new("shift", transform = function(x) x - d, inverse = function(x) x + d)
}
ggplot(data=pdata, aes(x=gene, y=adjrsquared)) +
  geom_bar( stat = "identity", width = 0.07 ,color = "steelblue") +
  theme_minimal() +
  theme(axis.text = element_text(size = 13),axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))+
  scale_y_continuous(trans = shift_trans(0)) +  
  geom_hline(aes(yintercept=0.0175), colour="#990000", linetype="dashed")
  
