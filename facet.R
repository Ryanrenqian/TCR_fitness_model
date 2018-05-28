options(echo=TRUE)
Args <- commandArgs(trailingOnly = TRUE)
print(Args)
library(facets)
datafile=Args[1]
output=Args[2]
rcmat = readSnpMatrix(datafile)
set.seed(1234)
xx = preProcSample(rcmat)
oo=procSample(xx,cval=150)
fit=emcncf(oo)
jpeg(file=paste(output,".jpeg",sep = ''))
print(paste(output,".jpeg",sep=''))
plotSample(x=oo,emfit=fit)
dev.off()
# chromosome    start   end     copy_number     minor_cn        major_cn        cellular_prevalence
#cnv=matrix()
cncf=fit$cncf
chromosome=cncf$chrom
start=cncf$start
end=cncf$end
copy_number=cncf$tcn.em
minor_cn=cncf$lcn.em
major_cn=cncf$tcn.em-cncf$lcn.em
cellular_prevalence=cncf$cf.em
cnv=data.frame(chromosome,start,end,copy_number,minor_cn,major_cn,cellular_prevalence)
write.table(cnv,file=paste(output,'.txt',sep=''),quote = F,row.names=F)
#print(paste(fit$purity,fit$ploidy,sep=' '))
cat(c(fit$purity,fit$ploidy),file=paste(output,'_fit.txt',sep=''),fill = T)

