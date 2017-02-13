#MCC
inp <- scan("5/NCBIblast_W7_541010_top_26.mcc",list(mcc=0.0))
NCBIblastWU<-inp$mcc
inp <- scan("5/WUblast_W7_541010_top_86.mcc",list(mcc=0.0))
WUblastWU<-inp$mcc
inp <- scan("5/fasta_541010_z11_top_ktup6_85.mcc",list(mcc=0.0))
fastaWU<-inp$mcc
inp <- scan("5/paralign_541010_b1000_top_91.mcc",list(mcc=0.0))
paralignWU<-inp$mcc
inp <- scan("5/ssearch_541010_z11_top_91.mcc",list(mcc=0.0))
ssearchWU<-inp$mcc

inp <- scan("5/NCBIblast_W7_541010_top_26.mcc",list(mcc=0.0))
NCBIblast<-inp$mcc
inp <- scan("5/WUblast_W3_top_108.mcc",list(mcc=0.0))
WUblast<-inp$mcc
inp <- scan("5/fasta_z11_top_89.mcc",list(mcc=0.0))
fasta<-inp$mcc
inp <- scan("5/paralign_541010_b1000_top_91.mcc",list(mcc=0.0))
paralign<-inp$mcc
inp <- scan("5/ssearch_z11_top_97.mcc",list(mcc=0.0))
ssearch<-inp$mcc

inp <- scan("5/HMMer2_hmmfs_-3.mcc",list(mcc=0.0))
HMMer<-inp$mcc
inp <- scan("5/SAM35_sw2_0306_-6.28.mcc",list(mcc=0.0))
SAM<-inp$mcc

inp <- scan("5/erpin_PCW0.1.mcc",list(mcc=0.0))
erpin<-inp$mcc
inp <- scan("5/erpin_PSR.mcc",list(mcc=0.0))
erpinPSR<-inp$mcc
inp <- scan("5/infernal7-local_11.74.mcc",list(mcc=0.0))
infernal<-inp$mcc
inp <- scan("5/ravenna2f-ML.mcc",list(mcc=0.0))
ravenna<-inp$mcc
inp <- scan("5/rsearch.mcc",list(mcc=0.0))
rsearch<-inp$mcc
inp <- scan("5/RSmatch12_15.mcc",list(mcc=0.0))
RSmatch12<-inp$mcc

#RANKS
inp <- scan("5/NCBIblast_W7_541010_top_26.ranks",list(mcc=0.0))
NCBIblastWUranks<-inp$mcc
inp <- scan("5/WUblast_W7_541010_top_86.ranks",list(mcc=0.0))
WUblastWUranks<-inp$mcc
inp <- scan("5/fasta_541010_z11_top_ktup6_85.ranks",list(mcc=0.0))
fastaWUranks<-inp$mcc
inp <- scan("5/paralign_541010_b1000_top_91.ranks",list(mcc=0.0))
paralignWUranks<-inp$mcc
inp <- scan("5/ssearch_541010_z11_top_91.ranks",list(mcc=0.0))
ssearchWUranks<-inp$mcc

inp <- scan("5/NCBIblast_W7_541010_top_26.ranks",list(mcc=0.0))
NCBIblastranks<-inp$mcc
inp <- scan("5/WUblast_W3_top_108.ranks",list(mcc=0.0))
WUblastranks<-inp$mcc
inp <- scan("5/fasta_z11_top_89.ranks",list(mcc=0.0))
fastaranks<-inp$mcc
inp <- scan("5/paralign_541010_b1000_top_91.ranks",list(mcc=0.0))
paralignranks<-inp$mcc
inp <- scan("5/ssearch_z11_top_97.ranks",list(mcc=0.0))
ssearchranks<-inp$mcc

inp <- scan("5/HMMer2_hmmfs_-3.ranks",list(mcc=0.0))
HMMerranks<-inp$mcc
inp <- scan("5/SAM35_sw2_0306_-6.28.ranks",list(mcc=0.0))
SAMranks<-inp$mcc

inp <- scan("5/erpin_PCW0.1.ranks",list(mcc=0.0))
erpinranks<-inp$mcc
inp <- scan("5/erpin_PSR.ranks",list(mcc=0.0))
erpinPSRranks<-inp$mcc
inp <- scan("5/infernal7-local_11.74.ranks",list(mcc=0.0))
infernalranks<-inp$mcc
inp <- scan("5/ravenna2f-ML.ranks",list(mcc=0.0))
ravennaranks<-inp$mcc
inp <- scan("5/rsearch.ranks",list(mcc=0.0))
rsearchranks<-inp$mcc
inp <- scan("5/RSmatch12_15.ranks",list(mcc=0.0))
RSmatch12ranks<-inp$mcc

#SENS
inp <- scan("5/NCBIblast_W7_541010_top_26.tps",list(id1="",id2="",id3="",mcc=0.0))
NCBIblastWUsens<-inp$mcc
inp <- scan("5/WUblast_W7_541010_top_86.tps",list(id1="",id2="",id3="",mcc=0.0))
WUblastWUsens<-inp$mcc
inp <- scan("5/fasta_541010_z11_top_ktup6_85.tps",list(id1="",id2="",id3="",mcc=0.0))
fastaWUsens<-inp$mcc
inp <- scan("5/paralign_541010_b1000_top_91.tps",list(id1="",id2="",id3="",mcc=0.0))
paralignWUsens<-inp$mcc
inp <- scan("5/ssearch_541010_z11_top_91.tps",list(id1="",id2="",id3="",mcc=0.0))
ssearchWUsens<-inp$mcc

inp <- scan("5/NCBIblast_W7_541010_top_26.tps",list(id1="",id2="",id3="",mcc=0.0))
NCBIblastsens<-inp$mcc
inp <- scan("5/WUblast_W3_top_108.tps",list(id1="",id2="",id3="",mcc=0.0))
WUblastsens<-inp$mcc
inp <- scan("5/fasta_z11_top_89.tps",list(id1="",id2="",id3="",mcc=0.0))
fastasens<-inp$mcc
inp <- scan("5/paralign_541010_b1000_top_91.tps",list(id1="",id2="",id3="",mcc=0.0))
paralignsens<-inp$mcc
inp <- scan("5/ssearch_z11_top_97.tps",list(id1="",id2="",id3="",mcc=0.0))
ssearchsens<-inp$mcc

inp <- scan("5/HMMer2_hmmfs_-3.tps",list(id1="",id2="",id3="",mcc=0.0))
HMMersens<-inp$mcc
inp <- scan("5/SAM35_sw2_0306_-6.28.tps",list(id1="",id2="",id3="",mcc=0.0))
SAMsens<-inp$mcc

inp <- scan("5/erpin_PCW0.1.tps",list(id1="",id2="",id3="",mcc=0.0))
erpinsens<-inp$mcc
inp <- scan("5/erpin_PSR.tps",list(id1="",id2="",id3="",mcc=0.0))
erpinPSRsens<-inp$mcc
inp <- scan("5/infernal7-local_11.74.tps",list(id1="",id2="",id3="",mcc=0.0))
infernalsens<-inp$mcc
inp <- scan("5/ravenna2f-ML.tps",list(id1="",id2="",id3="",mcc=0.0))
ravennasens<-inp$mcc
inp <- scan("5/rsearch.tps",list(id1="",id2="",id3="",mcc=0.0,id4=""))
rsearchsens<-inp$mcc
inp <- scan("5/RSmatch12_15.tps",list(id1="",id2="",id3="",mcc=0.0))
RSmatch12sens<-inp$mcc

#SPEC
inp <- scan("5/NCBIblast_W7_541010_top_26.fps",list(id1="",id2="",id3="",mcc=0.0))
NCBIblastWUspec<-1-inp$mcc
inp <- scan("5/WUblast_W7_541010_top_86.fps",list(id1="",id2="",id3="",mcc=0.0))
WUblastWUspec<-1-inp$mcc
inp <- scan("5/fasta_541010_z11_top_ktup6_85.fps",list(id1="",id2="",id3="",mcc=0.0))
fastaWUspec<-1-inp$mcc
inp <- scan("5/paralign_541010_b1000_top_91.fps",list(id1="",id2="",id3="",mcc=0.0))
paralignWUspec<-1-inp$mcc
inp <- scan("5/ssearch_541010_z11_top_91.fps",list(id1="",id2="",id3="",mcc=0.0))
ssearchWUspec<-1-inp$mcc

inp <- scan("5/NCBIblast_W7_541010_top_26.fps",list(id1="",id2="",id3="",mcc=0.0))
NCBIblastspec<-1-inp$mcc
inp <- scan("5/WUblast_W3_top_108.fps",list(id1="",id2="",id3="",mcc=0.0))
WUblastspec<-1-inp$mcc
inp <- scan("5/fasta_z11_top_89.fps",list(id1="",id2="",id3="",mcc=0.0))
fastaspec<-1-inp$mcc
inp <- scan("5/paralign_541010_b1000_top_91.fps",list(id1="",id2="",id3="",mcc=0.0))
paralignspec<-1-inp$mcc
inp <- scan("5/ssearch_z11_top_97.fps",list(id1="",id2="",id3="",mcc=0.0))
ssearchspec<-1-inp$mcc


inp <- scan("5/HMMer2_hmmfs_-3.fps",list(id1="",id2="",id3="",mcc=0.0))
HMMerspec<-1-inp$mcc
inp <- scan("5/SAM35_sw2_0306_-6.28.fps",list(id1="",id2="",id3="",mcc=0.0))
SAMspec<-1-inp$mcc

inp <- scan("5/erpin_PCW0.1.fps",list(id1="",id2="",id3="",mcc=0.0))
erpinspec<-1-inp$mcc
inp <- scan("5/erpin_PSR.fps",list(id1="",id2="",id3="",mcc=0.0))
erpinPSRspec<-1-inp$mcc
inp <- scan("5/infernal7-local_11.74.fps",list(id1="",id2="",id3="",mcc=0.0))
infernalspec<-1-inp$mcc
inp <- scan("5/ravenna2f-ML.fps",list(id1="",id2="",id3="",mcc=0.0))
ravennaspec<-1-inp$mcc
inp <- scan("5/rsearch.fps",list(id1="",id2="",id3="",mcc=0.0,id4=""))
rsearchspec<-1-inp$mcc
inp <- scan("5/RSmatch12_15.fps",list(id1="",id2="",id3="",mcc=0.0))
RSmatch12spec<-1-inp$mcc

## x11()
## op <- par(mfrow=c(1,1),las=2,cex.axis=1.0)
## boxplot(NCBIblast,WUblast,fasta,paralign,ssearch,HMMer,SAM,erpin,infernal,ravenna,rsearch,RSmatch12,names=c("NCBI-BLAST","WU-BLAST","FASTA","ParAlign","SSEARCH","HMMer","SAM","ERPIN","Infernal","RaveNnA","RSEARCH","RSmatch"),varwidth=TRUE,notch=TRUE,ylab="MCC",ylim=c(0,1),col=c(6,7,5,"sienna","cyan4",3,4,1,2,"cornflowerblue",8,"orange"),main="5 sequences")
## dev.copy2eps(file="boxplot5.eps")





###############

mn.MCC <- c(median(infernal),median(rsearch),median(ravenna),median(SAM),median(HMMer),median(WUblast),median(ssearch),median(fasta),median(NCBIblast),median(paralign),median(RSmatch12),median(erpin))

mn.SENS <- c(median(infernalsens),median(rsearchsens),median(ravennasens),median(SAMsens),median(HMMersens),median(WUblastsens),median(ssearchsens),median(fastasens),median(NCBIblastsens),median(paralignsens),median(RSmatch12sens),median(erpinsens))

mn.SPEC <- c(median(infernalspec),median(rsearchspec),median(ravennaspec),median(SAMspec),median(HMMerspec),median(WUblastspec),median(ssearchspec),median(fastaspec),median(NCBIblastspec),median(paralignspec),median(RSmatch12spec),median(erpinspec))

x11()
op <- par(mfrow=c(1,1),las=2,cex.axis=1.0)
boxplot(infernalranks,rsearchranks,ravennaranks,SAMranks,HMMerranks,WUblastranks,ssearchranks,fastaranks,NCBIblastranks,paralignranks,RSmatch12ranks,erpinranks,names=c("Infernal","RSEARCH","RaveNnA","SAM","HMMer","WU-BLAST","SSEARCH","FASTA","NCBI-BLAST","ParAlign","RSmatch","ERPIN"),varwidth=TRUE,notch=TRUE,outline=FALSE,ylab="rank(MCC)",ylim=c(35,1),col=c(8,8,8,"gray30","gray30",0,0,0,0,0,8,8),main="5 sequences")

X<-seq(1+0.45,12-0.45, length.out=12 )

op <- par(new=T,las=0)
plot(X,mn.MCC,axes=F,ylim=c(0,1),xlim=c(1,12),xlab="",ylab="",pch="*",col="red", cex = 2)
lines(X,mn.MCC,col="red", lwd = 2)
axis(side=3, labels = FALSE, tick = FALSE)
axis(side=4)
mtext("MCC",side=4,line=1,col="red")

par(new=T)
plot(X,mn.SENS,axes=F,ylim=c(0,1),xlim=c(1,12),xlab="",ylab="",pch="*",col="blue", cex = 2)
lines(X,mn.SENS,col="blue", lwd = 2)

par(new=T)
plot(X,mn.SPEC,axes=F,ylim=c(0,1),xlim=c(1,12),xlab="",ylab="",pch="*",col="green", cex = 2)
lines(X,mn.SPEC,col="green", lwd = 2)

dev.copy2eps(file="boxplot5_ranks_nocol.eps")

x11()
op <- par(mfrow=c(1,1),las=2,cex.axis=1.0)
boxplot(infernal,rsearch,ravenna,SAM,HMMer,WUblast,ssearch,fasta,NCBIblast,paralign,RSmatch12,erpin,names=c("Infernal","RSEARCH","RaveNnA","SAM","HMMer","WU-BLAST","SSEARCH","FASTA","NCBI-BLAST","ParAlign","RSmatch","ERPIN"),varwidth=TRUE,notch=TRUE,ylab="MCC",ylim=c(0,1),col=c(8,8,8,"gray30","gray30",0,0,0,0,0,8,8),main="5 sequences",outline=FALSE)
dev.copy2eps(file="boxplot5_nocol.eps")

x11()
op <- par(las=2,cex.axis=1.0)
boxplot(infernalsens,rsearchsens,ravennasens,SAMsens,HMMersens,WUblastsens,ssearchsens,fastasens,NCBIblastsens,paralignsens,RSmatch12sens,erpinsens,names=c("Infernal","RSEARCH","RaveNnA","SAM","HMMer","WU-BLAST","SSEARCH","FASTA","NCBI-BLAST","ParAlign","RSmatch","ERPIN"),varwidth=TRUE,notch=TRUE,ylab="Sensitivity",ylim=c(0,1),col=c(8,8,8,"gray30","gray30",0,0,0,0,0,8,8),main="5 sequences",outline=FALSE)
dev.copy2eps(file="boxplot5_sens_nocol.eps")

x11()
op <- par(las=2,cex.axis=1.0)
boxplot(infernalspec,rsearchspec,ravennaspec,SAMspec,HMMerspec,WUblastspec,ssearchspec,fastaspec,NCBIblastspec,paralignspec,RSmatch12spec,erpinspec,names=c("Infernal","RSEARCH","RaveNnA","SAM","HMMer","WU-BLAST","SSEARCH","FASTA","NCBI-BLAST","ParAlign","RSmatch","ERPIN"),varwidth=TRUE,notch=TRUE,ylab="Specificity",ylim=c(0.992,1),col=c(8,8,8,"gray30","gray30",0,0,0,0,0,8,8),main="5 sequences",outline=FALSE)
dev.copy2eps(file="boxplot5_spec_nocol.eps")

#####20 sequence datasets######
#MCC
inp <- scan("20/NCBIblast_W7_541010_top_26.mcc",list(mcc=0.0))
NCBIblastWU20<-inp$mcc
inp <- scan("20/WUblast_W7_541010_top_86.mcc",list(mcc=0.0))
WUblastWU20<-inp$mcc
inp <- scan("20/fasta_541010_z11_top_ktup6_85.mcc",list(mcc=0.0))
fastaWU20<-inp$mcc
inp <- scan("20/paralign_541010_b1000_top_91.mcc",list(mcc=0.0))
paralignWU20<-inp$mcc
inp <- scan("20/ssearch_541010_z11_top_91.mcc",list(mcc=0.0))
ssearchWU20<-inp$mcc

inp <- scan("20/NCBIblast_W7_541010_top_26.mcc",list(mcc=0.0))
NCBIblast20<-inp$mcc
inp <- scan("20/WUblast_W3_top_108.mcc",list(mcc=0.0))
WUblast20<-inp$mcc
inp <- scan("20/fasta_z11_top_89.mcc",list(mcc=0.0))
fasta20<-inp$mcc
inp <- scan("20/paralign_541010_b1000_top_91.mcc",list(mcc=0.0))
paralign20<-inp$mcc
inp <- scan("20/ssearch_z11_top_97.mcc",list(mcc=0.0))
ssearch20<-inp$mcc

inp <- scan("20/HMMer2_hmmfs_-3.mcc",list(mcc=0.0))
HMMer20<-inp$mcc
inp <- scan("20/SAM35_sw2_0306_-6.28.mcc",list(mcc=0.0))
SAM20<-inp$mcc

inp <- scan("20/erpin_PCW0.1_20.mcc",list(mcc=0.0))
erpin20<-inp$mcc
inp <- scan("20/infernal7-local_11.74.mcc",list(mcc=0.0))
infernal20<-inp$mcc
inp <- scan("20/ravenna2f-ML.mcc",list(mcc=0.0))
ravenna20<-inp$mcc
inp <- scan("20/rsearchD100s10n1000.mcc",list(mcc=0.0))
rsearch20<-inp$mcc
inp <- scan("20/RSmatch12_15.mcc",list(mcc=0.0))
RSmatch1220<-inp$mcc

#RANKs
inp <- scan("20/NCBIblast_W7_541010_top_26.ranks",list(mcc=0.0))
NCBIblastWU20ranks<-inp$mcc
inp <- scan("20/WUblast_W7_541010_top_86.ranks",list(mcc=0.0))
WUblastWU20ranks<-inp$mcc
inp <- scan("20/fasta_541010_z11_top_ktup6_85.ranks",list(mcc=0.0))
fastaWU20ranks<-inp$mcc
inp <- scan("20/paralign_541010_b1000_top_91.ranks",list(mcc=0.0))
paralignWU20ranks<-inp$mcc
inp <- scan("20/ssearch_541010_z11_top_91.ranks",list(mcc=0.0))
ssearchWU20ranks<-inp$mcc

inp <- scan("20/NCBIblast_W7_541010_top_26.ranks",list(mcc=0.0))
NCBIblast20ranks<-inp$mcc
inp <- scan("20/WUblast_W3_top_108.ranks",list(mcc=0.0))
WUblast20ranks<-inp$mcc
inp <- scan("20/fasta_z11_top_89.ranks",list(mcc=0.0))
fasta20ranks<-inp$mcc
inp <- scan("20/paralign_541010_b1000_top_91.ranks",list(mcc=0.0))
paralign20ranks<-inp$mcc
inp <- scan("20/ssearch_z11_top_97.ranks",list(mcc=0.0))
ssearch20ranks<-inp$mcc

inp <- scan("20/HMMer2_hmmfs_-3.ranks",list(mcc=0.0))
HMMer20ranks<-inp$mcc
inp <- scan("20/SAM35_sw2_0306_-6.28.ranks",list(mcc=0.0))
SAM20ranks<-inp$mcc

inp <- scan("20/erpin_PCW0.1_20.ranks",list(mcc=0.0))
erpin20ranks<-inp$mcc
inp <- scan("20/infernal7-local_11.74.ranks",list(mcc=0.0))
infernal20ranks<-inp$mcc
inp <- scan("20/ravenna2f-ML.ranks",list(mcc=0.0))
ravenna20ranks<-inp$mcc
inp <- scan("20/rsearchD100s10n1000.ranks",list(mcc=0.0))
rsearch20ranks<-inp$mcc
inp <- scan("20/RSmatch12_15.ranks",list(mcc=0.0))
RSmatch1220ranks<-inp$mcc

#SENS
inp <- scan("20/NCBIblast_W7_541010_top_26.tps",list(id1="",id2="",id3="",mcc=0.0))
NCBIblastWU20sens<-inp$mcc
inp <- scan("20/WUblast_W7_541010_top_86.tps",list(id1="",id2="",id3="",mcc=0.0))
WUblastWU20sens<-inp$mcc
inp <- scan("20/fasta_541010_z11_top_ktup6_85.tps",list(id1="",id2="",id3="",mcc=0.0))
fastaWU20sens<-inp$mcc
inp <- scan("20/paralign_541010_b1000_top_91.tps",list(id1="",id2="",id3="",mcc=0.0))
paralignWU20sens<-inp$mcc
inp <- scan("20/ssearch_541010_z11_top_91.tps",list(id1="",id2="",id3="",mcc=0.0))
ssearchWU20sens<-inp$mcc

inp <- scan("20/NCBIblast_W7_541010_top_26.tps",list(id1="",id2="",id3="",mcc=0.0))
NCBIblast20sens<-inp$mcc
inp <- scan("20/WUblast_W3_top_108.tps",list(id1="",id2="",id3="",mcc=0.0))
WUblast20sens<-inp$mcc
inp <- scan("20/fasta_z11_top_89.tps",list(id1="",id2="",id3="",mcc=0.0))
fasta20sens<-inp$mcc
inp <- scan("20/paralign_541010_b1000_top_91.tps",list(id1="",id2="",id3="",mcc=0.0))
paralign20sens<-inp$mcc
inp <- scan("20/ssearch_z11_top_97.tps",list(id1="",id2="",id3="",mcc=0.0))
ssearch20sens<-inp$mcc

inp <- scan("20/HMMer2_hmmfs_-3.tps",list(id1="",id2="",id3="",mcc=0.0))
HMMer20sens<-inp$mcc
inp <- scan("20/SAM35_sw2_0306_-6.28.tps",list(id1="",id2="",id3="",mcc=0.0))
SAM20sens<-inp$mcc

inp <- scan("20/erpin_PCW0.1_20.tps",list(id1="",id2="",id3="",mcc=0.0))
erpin20sens<-inp$mcc
inp <- scan("20/infernal7-local_11.74.tps",list(id1="",id2="",id3="",mcc=0.0))
infernal20sens<-inp$mcc
inp <- scan("20/ravenna2f-ML.tps",list(id1="",id2="",id3="",mcc=0.0))
ravenna20sens<-inp$mcc
inp <- scan("20/rsearchD100s10n1000.tps",list(id1="",id2="",id3="",mcc=0.0))
rsearch20sens<-inp$mcc
inp <- scan("20/RSmatch12_15.tps",list(id1="",id2="",id3="",mcc=0.0))
RSmatch1220sens<-inp$mcc

#SPEC
inp <- scan("20/NCBIblast_W7_541010_top_26.fps",list(id1="",id2="",id3="",mcc=0.0))
NCBIblastWU20spec<-1-inp$mcc
inp <- scan("20/WUblast_W7_541010_top_86.fps",list(id1="",id2="",id3="",mcc=0.0))
WUblastWU20spec<-1-inp$mcc
inp <- scan("20/fasta_541010_z11_top_ktup6_85.fps",list(id1="",id2="",id3="",mcc=0.0))
fastaWU20spec<-1-inp$mcc
inp <- scan("20/paralign_541010_b1000_top_91.fps",list(id1="",id2="",id3="",mcc=0.0))
paralignWU20spec<-1-inp$mcc
inp <- scan("20/ssearch_541010_z11_top_91.fps",list(id1="",id2="",id3="",mcc=0.0))
ssearchWU20spec<-1-inp$mcc

inp <- scan("20/NCBIblast_W7_541010_top_26.fps",list(id1="",id2="",id3="",mcc=0.0))
NCBIblast20spec<-1-inp$mcc
inp <- scan("20/WUblast_W3_top_108.fps",list(id1="",id2="",id3="",mcc=0.0))
WUblast20spec<-1-inp$mcc
inp <- scan("20/fasta_z11_top_89.fps",list(id1="",id2="",id3="",mcc=0.0))
fasta20spec<-1-inp$mcc
inp <- scan("20/paralign_541010_b1000_top_91.fps",list(id1="",id2="",id3="",mcc=0.0))
paralign20spec<-1-inp$mcc
inp <- scan("20/ssearch_z11_top_97.fps",list(id1="",id2="",id3="",mcc=0.0))
ssearch20spec<-1-inp$mcc

inp <- scan("20/HMMer2_hmmfs_-3.fps",list(id1="",id2="",id3="",mcc=0.0))
HMMer20spec<-1-inp$mcc
inp <- scan("20/SAM35_sw2_0306_-6.28.fps",list(id1="",id2="",id3="",mcc=0.0))
SAM20spec<-1-inp$mcc

inp <- scan("20/erpin_PCW0.1_20.fps",list(id1="",id2="",id3="",mcc=0.0))
erpin20spec<-1-inp$mcc
inp <- scan("20/infernal7-local_11.74.fps",list(id1="",id2="",id3="",mcc=0.0))
infernal20spec<-1-inp$mcc
inp <- scan("20/ravenna2f-ML.fps",list(id1="",id2="",id3="",mcc=0.0))
ravenna20spec<-1-inp$mcc
inp <- scan("20/rsearchD100s10n1000.fps",list(id1="",id2="",id3="",mcc=0.0))
rsearch20spec<-1-inp$mcc
inp <- scan("20/RSmatch12_15.fps",list(id1="",id2="",id3="",mcc=0.0))
RSmatch1220spec<-1-inp$mcc


mn.MCC <- c(median(infernal20),median(rsearch20),median(ravenna20),median(SAM20),median(HMMer20),median(WUblast20),median(ssearch20),median(fasta20),median(NCBIblast20),median(paralign20),median(RSmatch1220),median(erpin20))

mn.SENS <- c(median(infernal20sens),median(rsearch20sens),median(ravenna20sens),median(SAM20sens),median(HMMer20sens),median(WUblast20sens),median(ssearch20sens),median(fasta20sens),median(NCBIblast20sens),median(paralign20sens),median(RSmatch1220sens),median(erpin20sens))

mn.SPEC <- c(median(infernal20spec),median(rsearch20spec),median(ravenna20spec),median(SAM20spec),median(HMMer20spec),median(WUblast20spec),median(ssearch20spec),median(fasta20spec),median(NCBIblast20spec),median(paralign20spec),median(RSmatch1220spec),median(erpin20spec))


x11()
op <- par(mfrow=c(1,1),las=2,cex.axis=1.0)
boxplot(infernal20ranks,rsearch20ranks,ravenna20ranks,SAM20ranks,HMMer20ranks,WUblast20ranks,ssearch20ranks,fasta20ranks,NCBIblast20ranks,paralign20ranks,RSmatch1220ranks,erpin20ranks,names=c("Infernal","RSEARCH","RaveNnA","SAM","HMMer","WU-BLAST","SSEARCH","FASTA","NCBI-BLAST","ParAlign","RSmatch","ERPIN"),varwidth=TRUE,notch=TRUE,outline=FALSE,ylab="rank(MCC)",ylim=c(35,1),col=c(8,8,8,"gray30","gray30",0,0,0,0,0,8,8),main="20 sequences")

X<-seq(1+0.45,12-0.45, length.out=12 )

op <- par(new=T,las=0)
plot(X,mn.MCC,axes=F,ylim=c(0,1),xlim=c(1,12),xlab="",ylab="",pch="*",col="red", cex = 2)
lines(X,mn.MCC,col="red", lwd = 2)
axis(side=3, labels = FALSE, tick = FALSE)
axis(side=4)
mtext("MCC",side=4,line=1,col="red")

par(new=T)
plot(X,mn.SENS,axes=F,ylim=c(0,1),xlim=c(1,12),xlab="",ylab="",pch="*",col="blue", cex = 2)
lines(X,mn.SENS,col="blue", lwd = 2)

par(new=T)
plot(X,mn.SPEC,axes=F,ylim=c(0,1),xlim=c(1,12),xlab="",ylab="",pch="*",col="green", cex = 2)
lines(X,mn.SPEC,col="green", lwd = 2)

dev.copy2eps(file="boxplot20_ranks_nocol.eps")

x11()
op <- par(las=2,cex.axis=1.0)
boxplot(infernal20,rsearch20,ravenna20,SAM20,HMMer20,WUblast20,ssearch20,fasta20,NCBIblast20,paralign20,RSmatch1220,erpin20,names=c("Infernal","RSEARCH","RaveNnA","SAM","HMMer","WU-BLAST","SSEARCH","FASTA","NCBI-BLAST","ParAlign","RSmatch","ERPIN"),varwidth=TRUE,notch=TRUE,ylab="MCC",ylim=c(0,1),col=c(8,8,8,"gray30","gray30",0,0,0,0,0,8,8),main="20 sequences",outline=FALSE)
dev.copy2eps(file="boxplot20_nocol.eps")


x11()
op <- par(las=2,cex.axis=1.0)
boxplot(infernal20sens,rsearch20sens,ravenna20sens,SAM20sens,HMMer20sens,WUblast20sens,ssearch20sens,fasta20sens,NCBIblast20sens,paralign20sens,RSmatch1220sens,erpin20sens,names=c("Infernal","RSEARCH","RaveNnA","SAM","HMMer","WU-BLAST","SSEARCH","FASTA","NCBI-BLAST","ParAlign","RSmatch","ERPIN"),varwidth=TRUE,notch=TRUE,ylab="Sensitivity",ylim=c(0,1),col=c(8,8,8,"gray30","gray30",0,0,0,0,0,8,8),main="20 sequences",outline=FALSE)
dev.copy2eps(file="boxplot20_sens_nocol.eps")

x11()
op <- par(las=2,cex.axis=1.0)
boxplot(infernal20spec,rsearch20spec,ravenna20spec,SAM20spec,HMMer20spec,WUblast20spec,ssearch20spec,fasta20spec,NCBIblast20spec,paralign20spec,RSmatch1220spec,erpin20spec,names=c("Infernal","RSEARCH","RaveNnA","SAM","HMMer","WU-BLAST","SSEARCH","FASTA","NCBI-BLAST","ParAlign","RSmatch","ERPIN"),varwidth=TRUE,notch=TRUE,ylab="Specificity",ylim=c(0.975,1),col=c(8,8,8,"gray30","gray30",0,0,0,0,0,8,8),main="20 sequences",outline=FALSE)
dev.copy2eps(file="boxplot20_spec_nocol.eps")

######################################################################
##############################WU SCORING##############################
######################################################################

x11()
op <- par(mfrow=c(1,2),las=2,cex.axis=1.0)
boxplot(NCBIblastWUranks,WUblastWUranks,fastaWUranks,paralignWUranks,ssearchWUranks,names=c("NCBI-BLAST","WU-BLAST","FASTA","ParAlign","SSEARCH"),varwidth=TRUE,notch=TRUE,col=c(0,0,0,0,0),ylab="rank(MCC)",ylim=c(35,1),main="5 sequences",outline=FALSE)
boxplot(NCBIblastWU20ranks,WUblastWU20ranks,fastaWU20ranks,paralignWU20ranks,ssearchWU20ranks,names=c("NCBI-BLAST","WU-BLAST","FASTA","ParAlign","SSEARCH"),outline=FALSE,varwidth=TRUE,notch=TRUE,col=c(0,0,0,0,0),ylab="rank(MCC)",ylim=c(35,1),main="20 sequences")
dev.copy2eps(file="WU-scoring_ranks_nocol.eps")

mn.MCC <- c(median(ssearchWU),median(fastaWU),median(WUblastWU),median(NCBIblastWU),median(paralignWU))
mn.SENS <- c(median(ssearchWUsens),median(fastaWUsens),median(WUblastWUsens),median(NCBIblastWUsens),median(paralignWUsens))
mn.SPEC <- c(median(ssearchWUspec),median(fastaWUspec),median(WUblastWUspec),median(NCBIblastWUspec),median(paralignWUspec))

x11()
op <- par(mfrow=c(1,1),las=2,cex.axis=1.0)
boxplot(ssearchWUranks,fastaWUranks,WUblastWUranks,NCBIblastWUranks,paralignWUranks,names=c("SSEARCH","FASTA","WU-BLAST","NCBI-BLAST","ParAlign"),varwidth=TRUE,notch=TRUE,col=c(0,0,0,0,0),ylab="rank(MCC)",ylim=c(35,1),main="5 sequences",outline=FALSE)

X<-seq(1+0.45,5-0.45, length.out=5 )
op <- par(new=T,las=0)
plot(X,mn.MCC,axes=F,ylim=c(0.8,1),xlim=c(1,5),xlab="",ylab="",pch="*",col="red", cex = 2)
lines(X,mn.MCC,col="red", lwd = 2)
axis(side=3, labels = FALSE, tick = FALSE)
axis(side=4)
mtext("MCC",side=4,line=1,col="red")

par(new=T)
plot(X,mn.SENS,axes=F,ylim=c(0.8,1),xlim=c(1,5),xlab="",ylab="",pch="*",col="blue", cex = 2)
lines(X,mn.SENS,col="blue", lwd = 2)

par(new=T)
plot(X,mn.SPEC,axes=F,ylim=c(0.8,1),xlim=c(1,5),xlab="",ylab="",pch="*",col="green", cex = 2)
lines(X,mn.SPEC,col="green", lwd = 2)

dev.copy2eps(file="WU-scoring_ranks_nocol5.eps")

mn.MCC <- c(median(ssearchWU20),median(fastaWU20),median(WUblastWU20),median(NCBIblastWU20),median(paralignWU20))
mn.SENS <- c(median(ssearchWU20sens),median(fastaWU20sens),median(WUblastWU20sens),median(NCBIblastWU20sens),median(paralignWU20sens))
mn.SPEC <- c(median(ssearchWU20spec),median(fastaWU20spec),median(WUblastWU20spec),median(NCBIblastWU20spec),median(paralignWU20spec))

x11()
op <- par(mfrow=c(1,1),las=2,cex.axis=1.0)

boxplot(ssearchWU20ranks,fastaWU20ranks,WUblastWU20ranks,NCBIblastWU20ranks,paralignWU20ranks,names=c("SSEARCH","FASTA","WU-BLAST","NCBI-BLAST","ParAlign"),varwidth=TRUE,notch=TRUE,col=c(0,0,0,0,0),ylab="rank(MCC)",ylim=c(35,1),main="20 sequences",outline=FALSE)
#boxplot(NCBIblastWU20ranks,WUblastWU20ranks,fastaWU20ranks,paralignWU20ranks,ssearchWU20ranks,names=c("NCBI-BLAST","WU-BLAST","FASTA","ParAlign","SSEARCH"),outline=FALSE,varwidth=TRUE,notch=TRUE,col=c(0,0,0,0,0),ylab="rank(MCC)",ylim=c(35,1),main="20 sequences")

X<-seq(1+0.45,5-0.45, length.out=5 )
op <- par(new=T,las=0)
plot(X,mn.MCC,axes=F,ylim=c(0.92,1),xlim=c(1,5),xlab="",ylab="",pch="*",col="red", cex = 2)
lines(X,mn.MCC,col="red", lwd = 2)
axis(side=3, labels = FALSE, tick = FALSE)
axis(side=4)
mtext("MCC",side=4,line=1,col="red")

par(new=T)
plot(X,mn.SENS,axes=F,ylim=c(0.92,1),xlim=c(1,5),xlab="",ylab="",pch="*",col="blue", cex = 2)
lines(X,mn.SENS,col="blue", lwd = 2)

par(new=T)
plot(X,mn.SPEC,axes=F,ylim=c(0.92,1),xlim=c(1,5),xlab="",ylab="",pch="*",col="green", cex = 2)
lines(X,mn.SPEC,col="green", lwd = 2)

dev.copy2eps(file="WU-scoring_ranks_nocol20.eps")

######

x11()
op <- par(mfrow=c(1,1),las=2,cex.axis=1.0)
boxplot(ssearchWU,fastaWU,WUblastWU,NCBIblastWU,paralignWU,names=c("SSEARCH","FASTA","WU-BLAST","NCBI-BLAST","ParAlign"),varwidth=TRUE,notch=TRUE,col=c(0,0,0,0,0),ylab="MCC",ylim=c(0,1),main="5 sequences",outline=FALSE)
dev.copy2eps(file="WU-scoring_nocol.eps")

x11()
op <- par(mfrow=c(1,1),las=2,cex.axis=1.0)
boxplot(ssearchWU20,fastaWU20,WUblastWU20,NCBIblastWU20,paralignWU20,names=c("SSEARCH","FASTA","WU-BLAST","NCBI-BLAST","ParAlign"),varwidth=TRUE,notch=TRUE,col=c(0,0,0,0,0),ylab="MCC",ylim=c(0,1),main="20 sequences",outline=FALSE)
dev.copy2eps(file="WU-scoring20_nocol.eps")

x11()
op <- par(mfrow=c(1,1),las=2,cex.axis=1.0)
boxplot(ssearchWUsens,fastaWUsens,WUblastWUsens,NCBIblastWUsens,paralignWUsens,names=c("SSEARCH","FASTA","WU-BLAST","NCBI-BLAST","ParAlign"),varwidth=TRUE,notch=TRUE,col=c(0,0,0,0,0),ylab="Sensitivity",ylim=c(0,1),main="5 sequences",outline=FALSE)
dev.copy2eps(file="WU-scoring_sens_nocol.eps")

x11()
op <- par(mfrow=c(1,1),las=2,cex.axis=1.0)
boxplot(ssearchWU20sens,fastaWU20sens,WUblastWU20sens,NCBIblastWU20sens,paralignWU20sens,names=c("SSEARCH","FASTA","WU-BLAST","NCBI-BLAST","ParAlign"),varwidth=TRUE,notch=TRUE,col=c(0,0,0,0,0),ylab="Sensitivity",ylim=c(0,1),main="20 sequences",outline=FALSE)
dev.copy2eps(file="WU-scoring20_sens_nocol.eps")

x11()
op <- par(mfrow=c(1,1),las=2,cex.axis=1.0)
boxplot(ssearchWUspec,fastaWUspec,WUblastWUspec,NCBIblastWUspec,paralignWUspec,names=c("SSEARCH","FASTA","WU-BLAST","NCBI-BLAST","ParAlign"),varwidth=TRUE,notch=TRUE,col=c(0,0,0,0,0),ylab="Specificity",ylim=c(0.984,1),main="5 sequences",outline=FALSE)
dev.copy2eps(file="WU-scoring_spec_nocol.eps")

x11()
op <- par(mfrow=c(1,1),las=2,cex.axis=1.0)
boxplot(ssearchWU20spec,fastaWU20spec,WUblastWU20spec,NCBIblastWU20spec,paralignWU20spec,names=c("SSEARCH","FASTA","WU-BLAST","NCBI-BLAST","ParAlign"),varwidth=TRUE,notch=TRUE,col=c(0,0,0,0,0),ylab="Specificity",ylim=c(0.984,1),main="20 sequences",outline=FALSE)
dev.copy2eps(file="WU-scoring20_spec_nocol.eps")

###############################
##RNA scoring matrices:
#MCC
inp <- scan("5/fasta_z11_top_89.mcc",list(mcc=0.0))
fasta<-inp$mcc
inp <- scan("5/fasta_U_z11_152.mcc",list(mcc=0.0))
fastaU<-inp$mcc
inp <- scan("5/fasta_z11_R9540_1010_top_6.mcc",list(mcc=0.0))
fastaR<-inp$mcc
inp <- scan("5/fasta_z11_sdf_1010_top_392.mcc",list(mcc=0.0))
fastaF<-inp$mcc

inp <- scan("5/WUblast_W7_541010_top_86.mcc",list(mcc=0.0))
WUblastWU<-inp$mcc
inp <- scan("5/WUblast_W7_pupy_top_92.mcc",list(mcc=0.0))
WUblastPUPY<-inp$mcc

##############
#Rank
inp <- scan("5/fasta_z11_top_89.ranks",list(mcc=0.0))
fastaranks<-inp$mcc
inp <- scan("5/fasta_U_z11_152.ranks",list(mcc=0.0))
fastaUranks<-inp$mcc
inp <- scan("5/fasta_z11_R9540_1010_top_6.ranks",list(mcc=0.0))
fastaRranks<-inp$mcc
inp <- scan("5/fasta_z11_sdf_1010_top_392.ranks",list(mcc=0.0))
fastaFranks<-inp$mcc

inp <- scan("5/WUblast_W7_541010_top_86.ranks",list(mcc=0.0))
WUblastWUranks<-inp$mcc
inp <- scan("5/WUblast_W7_pupy_top_92.ranks",list(mcc=0.0))
WUblastPUPYranks<-inp$mcc

###############
#Sens
inp <- scan("5/fasta_z11_top_89.tps",list(id1="",id2="",id3="",mcc=0.0))
fastasens<-inp$mcc
inp <- scan("5/fasta_U_z11_152.tps",list(id1="",id2="",id3="",mcc=0.0))
fastaUsens<-inp$mcc
inp <- scan("5/fasta_z11_R9540_1010_top_6.tps",list(id1="",id2="",id3="",mcc=0.0))
fastaRsens<-inp$mcc
inp <- scan("5/fasta_z11_sdf_1010_top_392.tps",list(id1="",id2="",id3="",mcc=0.0))
fastaFsens<-inp$mcc

inp <- scan("5/WUblast_W7_541010_top_86.tps",list(id1="",id2="",id3="",mcc=0.0))
WUblastWUsens<-inp$mcc
inp <- scan("5/WUblast_W7_pupy_top_92.tps",list(id1="",id2="",id3="",mcc=0.0))
WUblastPUPYsens<-inp$mcc

##########
#Spec
inp <- scan("5/fasta_z11_top_89.fps",list(id1="",id2="",id3="",mcc=0.0))
fastaspec<-1-inp$mcc
inp <- scan("5/fasta_z11_top_89.fps",list(id1="",id2="",id3="",mcc=0.0))
fastaspec<-1-inp$mcc
inp <- scan("5/fasta_U_z11_152.fps",list(id1="",id2="",id3="",mcc=0.0))
fastaUspec<-1-inp$mcc
inp <- scan("5/fasta_z11_R9540_1010_top_6.fps",list(id1="",id2="",id3="",mcc=0.0))
fastaRspec<-1-inp$mcc
inp <- scan("5/fasta_z11_sdf_1010_top_392.fps",list(id1="",id2="",id3="",mcc=0.0))
fastaFspec<-1-inp$mcc

inp <- scan("5/WUblast_W7_541010_top_86.fps",list(id1="",id2="",id3="",mcc=0.0))
WUblastWUspec<-1-inp$mcc

inp <- scan("5/WUblast_W7_pupy_top_92.fps",list(id1="",id2="",id3="",mcc=0.0))
WUblastPUPYspec<-1-inp$mcc

######################################################################
#ASR sequences

inp <- scan("5/fasta34_top_89.rnanc.seq10asr.mcc",list(mcc=0.0))
fastaASR<-inp$mcc
inp <- scan("5/ssearch34_top_97.rnanc.seq10asr.mcc",list(mcc=0.0))
ssearchASR<-inp$mcc
inp <- scan("5/fasta34_top_89.rnanc.seq10asr.ranks",list(mcc=0.0))
fastaASRranks<-inp$mcc
inp <- scan("5/ssearch34_top_97.rnanc.seq10asr.ranks",list(mcc=0.0))
ssearchASRranks<-inp$mcc
inp <- scan("5/fasta34_top_89.rnanc.seq10asr.tps",list(id1="",id2="",id3="",mcc=0.0))
fastaASRsens<-inp$mcc
inp <- scan("5/ssearch34_top_97.rnanc.seq10asr.tps",list(id1="",id2="",id3="",mcc=0.0))
ssearchASRsens<-inp$mcc
inp <- scan("5/fasta34_top_89.rnanc.seq10asr.fps",list(id1="",id2="",id3="",mcc=0.0))
fastaASRspec<-1-inp$mcc
inp <- scan("5/ssearch34_top_97.rnanc.seq10asr.fps",list(id1="",id2="",id3="",mcc=0.0))
ssearchASRspec<-1-inp$mcc

######################################################################
###########################
##20 seq
#MCC
inp <- scan("20/fasta_z11_top_89.mcc",list(mcc=0.0))
fasta20<-inp$mcc
inp <- scan("20/fasta_U_z11_152.mcc",list(mcc=0.0))
fastaU20<-inp$mcc
inp <- scan("20/fasta_z11_R9540_1010_top_6.mcc",list(mcc=0.0))
fastaR20<-inp$mcc
inp <- scan("20/fasta_z11_sdf_1010_top_392.mcc",list(mcc=0.0))
fastaF20<-inp$mcc

inp <- scan("20/WUblast_W7_541010_top_86.mcc",list(mcc=0.0))
WUblastWU20<-inp$mcc
inp <- scan("20/WUblast_W7_pupy_top_92.mcc",list(mcc=0.0))
WUblastPUPY20<-inp$mcc
#Ranks:
inp <- scan("20/fasta_z11_top_89.ranks",list(mcc=0.0))
fasta20ranks<-inp$mcc
inp <- scan("20/fasta_U_z11_152.ranks",list(mcc=0.0))
fastaU20ranks<-inp$mcc
inp <- scan("20/fasta_z11_R9540_1010_top_6.ranks",list(mcc=0.0))
fastaR20ranks<-inp$mcc
inp <- scan("20/fasta_z11_sdf_1010_top_392.ranks",list(mcc=0.0))
fastaF20ranks<-inp$mcc

inp <- scan("20/WUblast_W7_541010_top_86.ranks",list(mcc=0.0))
WUblastWU20ranks<-inp$mcc
inp <- scan("20/WUblast_W7_pupy_top_92.ranks",list(mcc=0.0))
WUblastPUPY20ranks<-inp$mcc
###############
#Sens
inp <- scan("20/fasta_z11_top_89.tps",list(id1="",id2="",id3="",mcc=0.0))
fasta20sens<-inp$mcc
inp <- scan("20/fasta_U_z11_152.tps",list(id1="",id2="",id3="",mcc=0.0))
fastaU20sens<-inp$mcc
inp <- scan("20/fasta_z11_R9540_1010_top_6.tps",list(id1="",id2="",id3="",mcc=0.0))
fastaR20sens<-inp$mcc
inp <- scan("20/fasta_z11_sdf_1010_top_392.tps",list(id1="",id2="",id3="",mcc=0.0))
fastaF20sens<-inp$mcc

inp <- scan("20/WUblast_W7_541010_top_86.tps",list(id1="",id2="",id3="",mcc=0.0))
WUblastWU20sens<-inp$mcc
inp <- scan("20/WUblast_W7_pupy_top_92.tps",list(id1="",id2="",id3="",mcc=0.0))
WUblastPUPY20sens<-inp$mcc

##########
#Spec
inp <- scan("20/fasta_z11_top_89.fps",list(id1="",id2="",id3="",mcc=0.0))
fasta20spec<-1-inp$mcc

inp <- scan("20/fasta_z11_top_89.fps",list(id1="",id2="",id3="",mcc=0.0))
fasta20spec<-1-inp$mcc
inp <- scan("20/fasta_U_z11_152.fps",list(id1="",id2="",id3="",mcc=0.0))
fastaU20spec<-1-inp$mcc
inp <- scan("20/fasta_z11_R9540_1010_top_6.fps",list(id1="",id2="",id3="",mcc=0.0))
fastaR20spec<-1-inp$mcc
inp <- scan("20/fasta_z11_sdf_1010_top_392.fps",list(id1="",id2="",id3="",mcc=0.0))
fastaF20spec<-1-inp$mcc

inp <- scan("20/WUblast_W7_541010_top_86.fps",list(id1="",id2="",id3="",mcc=0.0))
WUblastWU20spec<-1-inp$mcc
inp <- scan("20/WUblast_W7_pupy_top_92.fps",list(id1="",id2="",id3="",mcc=0.0))
WUblastPUPY20spec<-1-inp$mcc

##################################################

## x11()
## op <- par(mfrow=c(1,2),las=2,cex.axis=1.0)
## boxplot(WUblastWU,WUblastPUPY,fasta,fastaU,fastaR,fastaF,names=c("WU-BLAST","WU-BLAST (PUPY)","FASTA","FASTA (U)","FASTA (RIBOSUM)","FASTA (FOLDALIGN)"),varwidth=TRUE,notch=TRUE,col=c(7,7,5,5,5,5),ylab="MCC",ylim=c(0,1),main="5 sequences")
## boxplot(WUblastWU20,WUblastPUPY20,fasta20,fastaU20,fastaR20,fastaF20,names=c("WU-BLAST","WU-BLAST (PUPY)","FASTA","FASTA (U)","FASTA (RIBOSUM)","FASTA (FOLDALIGN)"),varwidth=TRUE,notch=TRUE,col=c(7,7,5,5,5,5),ylab="MCC",ylim=c(0,1),main="20 sequences")
## dev.copy2eps(file="RNA-scoring.eps")

x11()
op <- par(mfrow=c(1,2),las=2,cex.axis=1.0)
boxplot(WUblastWU,WUblastPUPY,fasta,fastaU,fastaR,fastaF,names=c("WU-BLAST","WU-BLAST (PUPY)","FASTA","FASTA (U)","FASTA (RIBOSUM)","FASTA (FOLDALIGN)"),varwidth=TRUE,notch=TRUE,col=c(0,0,0,0,0,0),ylab="MCC",ylim=c(0,1),main="5 sequences",outline=FALSE)
boxplot(WUblastWU20,WUblastPUPY20,fasta20,fastaU20,fastaR20,fastaF20,names=c("WU-BLAST","WU-BLAST (PUPY)","FASTA","FASTA (U)","FASTA (RIBOSUM)","FASTA (FOLDALIGN)"),varwidth=TRUE,notch=TRUE,col=c(0,0,0,0,0,0),ylab="MCC",ylim=c(0,1),main="20 sequences",outline=FALSE)
dev.copy2eps(file="RNA-scoring_nocol.eps")


mn.MCC <- c(median(WUblastWU),median(WUblastPUPY),median(fasta),median(fastaU),median(fastaF),median(fastaR))
mn.SENS <- c(median(WUblastWUsens),median(WUblastPUPYsens),median(fastasens),median(fastaUsens),median(fastaFsens),median(fastaRsens))
mn.SPEC <- c(median(WUblastWUspec),median(WUblastPUPYspec),median(fastaspec),median(fastaUspec),median(fastaFspec),median(fastaRspec))

x11()
op <- par(mfrow=c(1,1),las=2,cex.axis=1.0)
boxplot(WUblastWUranks,WUblastPUPYranks,fastaranks,fastaUranks,fastaFranks,fastaRranks,names=c("WU-BLAST","WU-BLAST (PUPY)","FASTA","FASTA (U)","FASTA (FOLDALIGN)","FASTA (RIBOSUM)"),varwidth=TRUE,notch=TRUE,col=c(0,0,0,0,0,0),ylab="rank(MCC)",ylim=c(35,1),main="RNA centric scoring schemes",outline=FALSE)
#boxplot(WUblastWU20ranks,WUblastPUPY20ranks,fasta20ranks,fastaU20ranks,fastaR20ranks,fastaF20ranks,names=c("WU-BLAST","WU-BLAST (PUPY)","FASTA","FASTA (U)","FASTA (RIBOSUM)","FASTA (FOLDALIGN)"),varwidth=TRUE,notch=TRUE,col=c(0,0,0,0,0,0),ylab="rank(MCC)",ylim=c(35,1),main="20 sequences",outline=FALSE)

X<-seq(1+0.45,6-0.45, length.out=6 )
op <- par(new=T,las=0)
plot(X,mn.MCC,axes=F,ylim=c(0,1),xlim=c(1,6),xlab="",ylab="",pch="*",col="red", cex = 2)
lines(X[1:2],mn.MCC[1:2],col="red", lwd = 2)
lines(X[3:6],mn.MCC[3:6],col="red", lwd = 2)

axis(side=3, labels = FALSE, tick = FALSE)
axis(side=4)
mtext("MCC",side=4,line=1,col="red")

par(new=T)
plot(X,mn.SENS,axes=F,ylim=c(0,1),xlim=c(1,6),xlab="",ylab="",pch="*",col="blue", cex = 2)
lines(X[1:2],mn.SENS[1:2],col="blue", lwd = 2)
lines(X[3:6],mn.SENS[3:6],col="blue", lwd = 2)

par(new=T)
plot(X,mn.SPEC,axes=F,ylim=c(0,1),xlim=c(1,6),xlab="",ylab="",pch="*",col="green", cex = 2)
lines(X[1:2],mn.SPEC[1:2],col="green", lwd = 2)
lines(X[3:6],mn.SPEC[3:6],col="green", lwd = 2)

dev.copy2eps(file="RNA-scoring_ranks_nocol.eps")

x11()
op <- par(mfrow=c(1,2),las=2,cex.axis=1.0)
boxplot(WUblastWUsens,WUblastPUPYsens,fastasens,fastaUsens,fastaRsens,fastaFsens,names=c("WU-BLAST","WU-BLAST (PUPY)","FASTA","FASTA (U)","FASTA (RIBOSUM)","FASTA (FOLDALIGN)"),varwidth=TRUE,notch=TRUE,col=c(0,0,0,0,0,0),ylab="Sensitivity",ylim=c(0,1),main="5 sequences",outline=FALSE)
boxplot(WUblastWU20sens,WUblastPUPY20sens,fasta20sens,fastaU20sens,fastaR20sens,fastaF20sens,names=c("WU-BLAST","WU-BLAST (PUPY)","FASTA","FASTA (U)","FASTA (RIBOSUM)","FASTA (FOLDALIGN)"),varwidth=TRUE,notch=TRUE,col=c(0,0,0,0,0,0),ylab="Sensitivity",ylim=c(0,1),main="20 sequences",outline=FALSE)
dev.copy2eps(file="RNA-scoring_sens_nocol.eps")

x11()
op <- par(mfrow=c(1,2),las=2,cex.axis=1.0)
boxplot(WUblastWUspec,WUblastPUPYspec,fastaspec,fastaUspec,fastaRspec,fastaFspec,names=c("WU-BLAST","WU-BLAST (PUPY)","FASTA","FASTA (U)","FASTA (RIBOSUM)","FASTA (FOLDALIGN)"),varwidth=TRUE,notch=TRUE,col=c(0,0,0,0,0,0),ylab="Specificity",ylim=c(0.83,1),main="5 sequences",outline=FALSE)
boxplot(WUblastWU20spec,WUblastPUPY20spec,fasta20spec,fastaU20spec,fastaR20spec,fastaF20spec,names=c("WU-BLAST","WU-BLAST (PUPY)","FASTA","FASTA (U)","FASTA (RIBOSUM)","FASTA (FOLDALIGN)"),varwidth=TRUE,notch=TRUE,col=c(0,0,0,0,0,0),ylab="Specificity",ylim=c(0.83,1),main="20 sequences",outline=FALSE)
dev.copy2eps(file="RNA-scoring_spec_nocol.eps")

###########
#ASR

mn.MCC <- c(median(ssearch),median(ssearchASR),median(fasta),median(fastaASR),median(erpin),median(erpinPSR))
mn.SENS <- c(median(ssearchsens),median(ssearchASRsens),median(fastasens),median(fastaASRsens),median(erpinsens),median(erpinPSRsens))
mn.SPEC <- c(median(ssearchspec),median(ssearchASRspec),median(fastaspec),median(fastaASRspec),median(erpinspec),median(erpinPSRspec))


x11()
op <- par(mfrow=c(1,1),las=2,cex.axis=1.0)
boxplot(ssearchranks,ssearchASRranks,fastaranks,fastaASRranks,erpinranks,erpinPSRranks,names=c("SSEARCH","SSEARCH [ASR]","FASTA","FASTA [ASR]","ERPIN","ERPIN [PSR]"),varwidth=TRUE,notch=TRUE,col=c(0,0,0,0),ylab="rank(MCC)",ylim=c(35,1),main="Phylogenetic sequence reconstructions",outline=FALSE)
#boxplot(fasta,fastaASR,ssearch,ssearchASR,names=c("FASTA","FASTA [ASR]","SSEARCH","SSEARCH [ASR]"),varwidth=TRUE,notch=TRUE,col=c(0,0,0,0),ylab="MCC",ylim=c(0,1),main="",outline=FALSE)

#boxplot(fastasens,fastaASRsens,ssearchsens,ssearchASRsens,names=c("FASTA","FASTA [ASR]","SSEARCH","SSEARCH [ASR]"),varwidth=TRUE,notch=TRUE,col=c(0,0,0,0),ylab="Sensitivity",ylim=c(0,1),main="",outline=FALSE)
#boxplot(fastaspec,fastaASRspec,ssearchspec,ssearchASRspec,names=c("FASTA","FASTA [ASR]","SSEARCH","SSEARCH [ASR]"),varwidth=TRUE,notch=TRUE,col=c(0,0,0,0),ylab="Specificity",ylim=c(0.99,1),main="",outline=FALSE)

X<-seq(1+0.45,6-0.45, length.out=6 )
op <- par(new=T,las=0)
plot(X,mn.MCC,axes=F,ylim=c(0,1),xlim=c(1,6),xlab="",ylab="",pch="*",col="red", cex = 2)
lines(X[1:2],mn.MCC[1:2],col="red", lwd = 2)
lines(X[3:4],mn.MCC[3:4],col="red", lwd = 2)
lines(X[5:6],mn.MCC[5:6],col="red", lwd = 2)

axis(side=3, labels = FALSE, tick = FALSE)
axis(side=4)
mtext("MCC",side=4,line=1,col="red")

par(new=T)
plot(X,mn.SENS,axes=F,ylim=c(0,1),xlim=c(1,6),xlab="",ylab="",pch="*",col="blue", cex = 2)
lines(X[1:2],mn.SENS[1:2],col="blue", lwd = 2)
lines(X[3:4],mn.SENS[3:4],col="blue", lwd = 2)
lines(X[5:6],mn.SENS[5:6],col="blue", lwd = 2)

par(new=T)
plot(X,mn.SPEC,axes=F,ylim=c(0,1),xlim=c(1,6),xlab="",ylab="",pch="*",col="green", cex = 2)
lines(X[1:2],mn.SPEC[1:2],col="green", lwd = 2)
lines(X[3:4],mn.SPEC[3:4],col="green", lwd = 2)
lines(X[5:6],mn.SPEC[5:6],col="green", lwd = 2)

dev.copy2eps(file="ASR_results_nocol.eps")











