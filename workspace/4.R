#import dataset
locus = read.table("data/locus012.dat",header = T)
multi_pheno = read.table("data/multi_phenos.txt")

#calculate correlation coefficient between different phenotypes
corr = cor(multi_pheno)

#calculate with chi_suqare test ,same in the second answer
p = matrix( nrow = 10 , ncol = 9445 )

for (i in 1:10){
  d1 = which(multi_pheno[,i]==0) #healthy 
  d2 = which(multi_pheno[,i]==1) #disease
  for (j in 1:9445){
    data1 = locus[d1,j]#healthy
    data2 = locus[d2,j]#disease
    
    table = matrix( c( length(which((data2[]==0))) , length(which((data2[]==1))) , length(which((data2[]==2))) , length(which((data1[]==0))) , length(which((data1[]==1))) , length(which((data1[]==2))) ) , nrow = 2 , ncol = 3 ,byrow = T)
    
    p[i,j] = -log(chisq.test(table)$p.value,10)
  }
}
#release memory
rm(d1,d2,data1,data2,i,j,table)

#draw plots
pdata = data.frame(gene = rep(1:9445),t(p))

ggplot(pdata, aes(gene , y = pvalue, color = phenotype)) + 
  geom_point(aes(y = X1, col = "pheno1")) +
  geom_point(aes(y = X2, col = "pheno2")) +
  geom_point(aes(y = X3, col = "pheno3")) +
  geom_point(aes(y = X4, col = "pheno4")) +
  geom_point(aes(y = X5, col = "pheno5")) +
  geom_point(aes(y = X6, col = "pheno6")) +
  geom_point(aes(y = X7, col = "pheno7")) +
  geom_point(aes(y = X8, col = "pheno8")) +
  geom_point(aes(y = X9, col = "pheno9")) +
  geom_point(aes(y = X10, col = "pheno10")) +
  geom_hline(aes(yintercept=3), colour="#990000", linetype="dashed") +
  theme_minimal()
  #scale_color_brewer(palette="Dark2")
  
#step 3 
#find snps that have at least one pvalue > 3
count = apply( p , 2 , function(col) length( which(col[]>3) ))
which(count[]>0)
length(which(count[]>0))

#canonical analysis
library(CCP)
library(GGally)
library(CCA)

locus_list = locus[,which(count[]>0)] 

cancorr = cc( locus_list , multi_pheno )
rho <- cancorr$cor

#cca p_value
n <- dim(locus_list)[1]
p <- length(locus_list)
q <- length(multi_pheno)
p.asym(rho, n, p, q, tstat = "Hotelling")

#draw plots
ccadata = data.frame(locus = rep(1:55), fit_coef = c(cancorr$xcoef[,1],cancorr$xcoef[,2]) , attribute = c(rep('1',55),rep('2',55)))

ggplot(ccadata, aes(locus , y = fit_coef, color = attribute , group =attribute)) + geom_point()  + geom_line() + theme_minimal()

#find which locus
which(cancorr$xcoef[,2]>1)
