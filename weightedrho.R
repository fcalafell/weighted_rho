#This script computes a mutation-weighted version of rho for Y-STR haplotypes
setwd("c:/docum180612/Francesc/Y/DF27/ABC Ferrara")

#input file: haplotyes as rows, no individual or STR labels, no missing values
haps<-read.csv("J2a.csv", header=FALSE)

#mutation rates per STR and generation
mutationrates<-read.csv("mutationrates.csv", header=FALSE)

#male generation length
genlength<-30

nind<-nrow(haps)
nloci<-ncol(haps)

yearspermut<-array(0:0,nloci)
yearspermut<-genlength/mutationrates
totyearspermut<-genlength/sum(t(mutationrates))
meanmut<-mean(t(mutationrates))


jj<-c(1:nind)
ii<-c(1:nloci)
wdistl<-array(0:0,dim=c(nind,nloci))
rwdistl<-array(0:0,dim=c(nind,nloci))

hapdif<-array(0:0,dim=c(nind,nloci))
rhoh<-array(0:0,dim=nind)
rhoh2<-array(0:0,dim=nind)
rrhoh<-array(0:0,dim=nind)
rrhoh2<-array(0:0,dim=nind)
medianhap<-array(0:0,dim=nloci)

for (i in ii){
   medianhap[i]<-median(haps[,i])
 }


freqhap<-array(0:0,dim=nind)
nhapdif<-1
for (i in ii){
  hapdif[1,i]<-haps[1,i]
}
freqhap[1]<-1

if (nind>1){
for (z in c(2:nind)){
  trobat<-0
  zz<-1
  while (zz<=nhapdif){
    nmatch<-0
    for (i in ii){
      if (haps[z,i]==hapdif[zz,i]){nmatch<-nmatch+1}
          }
      if (nmatch==nloci) {freqhap[zz]<-freqhap[zz]+1
  trobat<-1
      break}
  else {zz<-zz+1}
  }
  if (trobat==0){
    nhapdif<-nhapdif+1
    for (i in ii){
      hapdif[nhapdif,i]<-haps[z,i]}
  freqhap[nhapdif]<-1}
  
}
}

for (i in ii){
  for (k in c(1:nhapdif)){ 
    wdistl[k,i]<-abs((hapdif[k,i]-medianhap[i]))*(meanmut/t(mutationrates[i]))
    rwdistl[k,i]<-abs((hapdif[k,i]-medianhap[i]))
  }
}
for (k in c(1:nhapdif)){
  rhoh[k]<-sum(wdistl[k,])*freqhap[k]
  rhoh2[k]<-sum(wdistl[k,])*freqhap[k]*freqhap[k]
  rrhoh[k]<-sum(rwdistl[k,])*freqhap[k]
  rrhoh2[k]<-sum(rwdistl[k,])*freqhap[k]*freqhap[k]
}
wrho<-sum(rhoh)/nind
wage<-wrho*totyearspermut
sdwrho<-sqrt(sum(rhoh2))/nind
sdwage<-sdwrho*totyearspermut

rho<-sum(rrhoh)/nind
age<-rho*totyearspermut
sdrho<-sqrt(sum(rrhoh2))/nind
sdage<-sdrho*totyearspermut


print (paste0("N=",round(nind, digits=0)))
print (paste0("Nhap=",round(nhapdif, digits=0)))

print (paste0("rho=",round(rho, digits=4)))
print (paste0("sdrho=",round(sdrho, digits=4)))       
print (paste0("age=",round(age, digits=0)))       
print (paste0("sdage=",round(sdage, digits=0)))

print (paste0("weighted rho=",round(wrho, digits=4)))
print (paste0("weighted sdrho=",round(sdwrho, digits=4)))       
print (paste0("weighted age=",round(wage, digits=0)))       
print (paste0("weighted sdage=",round(sdwage, digits=0)))




