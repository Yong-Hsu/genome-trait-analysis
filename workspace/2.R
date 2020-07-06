#import dataset from the last question
locus = read.table("snp.dat",header = T)

#create probability list to store data
p = matrix(nrow = 1,ncol = 9445)

#Chi-square test
for (i in 1:9445){
  data1 = locus[1:500,i]#healthy
  data2 = locus[501:1000,i]#disease
  
  table = matrix( c( length(which((data2[]==0))) , length(which((data2[]==1))) , length(which((data2[]==2))) , length(which((data1[]==0))) , length(which((data1[]==1))) , length(which((data1[]==2))) ) , nrow = 2 , ncol = 3 ,byrow = T)
  
  p[i] = -log(chisq.test(table)$p.value,10)
}

#deal with DD/DI/II special points
p[1,6070] = 0
p[1,8473] = 0
rm(data1,data2,i,table)

#find the snp
colnames(locus)[which(p[]>-log(0.05/9445,10))]

#draw plots
library(ggplot2)
pdata = data.frame(snp = rep(1:9445),pvalue = t(p))

ggplot(data=pdata, aes(x=snp, y=pvalue)) +
  geom_bar(stat = "identity", width = 0.05 ,color = "steelblue") +
  theme_minimal() +
  theme(axis.text = element_text(size = 20),axis.title.x=element_text(size=20),axis.title.y=element_text(size=20))+
  geom_hline(aes(yintercept=1.3013), colour="#990000", linetype="dashed")+
  geom_hline(aes(yintercept=5.275724), colour="#990000", linetype="dashed")+
  annotate("text", x = 1000, y = 5.5, label = "bonferroni correction")

