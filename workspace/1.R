#import dataset
getype = read.table("data/genotype.dat",header = T)

#use a list to record type 
type = matrix( nrow = 1 , ncol = 9445 )
for (n in 1:9445) {
  c = levels(getype[,n])
  type[n] = c[2]
}
rm(n,c)

for (n in 1:9445) { 
  if(type[n]== "AT" || type[n]=="TA")
  { 
    getype[,n] = gsub("AA",0,getype[,n])
    getype[,n] = gsub("AT",1,getype[,n])
    getype[,n] = gsub("TA",1,getype[,n])
    getype[,n] = gsub("TT",2,getype[,n]) 
  }
  if(type[n]== "AG" || type[n]=="GA")
  {
    getype[,n] = gsub("AA",0,getype[,n])
    getype[,n] = gsub("AG",1,getype[,n])
    getype[,n] = gsub("GA",1,getype[,n])
    getype[,n] = gsub("GG",2,getype[,n]) 
  }
  if(type[n]== "AC" || type[n]=="CA")
  {
    getype[,n] = gsub("AA",0,getype[,n])
    getype[,n] = gsub("AC",1,getype[,n])
    getype[,n] = gsub("CA",1,getype[,n])
    getype[,n] = gsub("CC",2,getype[,n]) 
  }
  if(type[n]== "CG" || type[n]=="GC")
  {
    getype[,n] = gsub("CC",0,getype[,n])
    getype[,n] = gsub("CG",1,getype[,n])
    getype[,n] = gsub("GC",1,getype[,n])
    getype[,n] = gsub("GG",2,getype[,n]) 
  }
  if(type[n]== "CT" || type[n]=="TC")
  {
    getype[,n] = gsub("CC",0,getype[,n])
    getype[,n] = gsub("CT",1,getype[,n])
    getype[,n] = gsub("TC",1,getype[,n])
    getype[,n] = gsub("TT",2,getype[,n]) 
  }
  if(type[n]== "GT" || type[n]=="TG")
  {
    getype[,n] = gsub("GG",0,getype[,n])
    getype[,n] = gsub("GT",1,getype[,n])
    getype[,n] = gsub("TG",1,getype[,n])
    getype[,n] = gsub("TT",2,getype[,n]) 
  }
  if(type[n]== "DI" || type[n]=="ID")
  {
    getype[,n] = gsub("DD",0,getype[,n])
    getype[,n] = gsub("DI",1,getype[,n])
    getype[,n] = gsub("ID",1,getype[,n])
    getype[,n] = gsub("II",2,getype[,n]) 
  }
}
rm(n)

locus = apply(getype,2,as.numeric)
write.table(locus,"snp.dat",col.names = TRUE)