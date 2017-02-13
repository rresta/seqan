inp<-scan("ave_pw_seq_id_noaln12_sorted.dat",list(n="",z=0.0))
id<-inp$z
inp<-scan("structural_rnaz_noaln12_sorted.dat",list(n="",z=0.0))
structZ<-inp$z

inp<-scan("clustal_rnaz_noaln12_sorted.dat",list(n="",z=0.0))
cluZ<-inp$z
inp<-scan("foldalign_rnaz_noaln12_sorted.dat",list(n="",z=0.0))
foldaZ<-inp$z

inp<-scan("foldalign2_scanning_fmat25_rnaz_noaln12_sorted.dat",list(n="",z=0.0))
folda2scanZ<-inp$z
inp<-scan("foldalign2_structure_fmat25_rnaz_noaln12_sorted.dat",list(n="",z=0.0))
folda2structZ<-inp$z

inp<-scan("lowSS_rnaz_noaln12_sorted.dat",list(n="",z=0.0))
folda2lowSSZ<-inp$z

inp<-scan("foldalign.global.fold_rnaz_noaln12_sorted.dat",list(n="",z=0.0))
folda2globalZ<-inp$z

inp<-scan("muscle_rnaz_noaln12_sorted.dat",list(n="",z=0.0))
musZ<-inp$z
inp<-scan("pmcomp-fast_rnaz_noaln12_sorted.dat",list(n="",z=0.0))
pcomfZ<-inp$z
inp<-scan("pmcomp_rnaz_noaln12_sorted.dat",list(n="",z=0.0))
pcomZ<-inp$z

inp<-scan("pcma_rnaz_noaln12_sorted.dat",list(n="",z=0.0))
pcmaZ<-inp$z

inp<-scan("poa_rnaz_noaln12_sorted.dat",list(n="",z=0.0))
poaZ<-inp$z

inp<-scan("proalign_rnaz_noaln12_sorted.dat",list(n="",z=0.0))
proalignZ<-inp$z

inp<-scan("prrn_rnaz_noaln12_sorted.dat",list(n="",z=0.0))
prrnZ<-inp$z

inp<-scan("dynalign_rnaz_noaln12_sorted.dat",list(n="",z=0.0))
dynalignZ<-inp$z

inp<-scan("stemloc-slow_rnaz_noaln12_sorted.dat",list(n="",z=0.0))
stemloc_slowZ<-inp$z
inp<-scan("stemloc-fast110_rnaz_noaln12_sorted.dat",list(n="",z=0.0))
stemloc_fastZ<-inp$z

inp<-scan("stral_rnaz_noaln12_sorted.dat",list(n="",z=0.0))
stralZ<-inp$z


X11()
plot(lowess(id,cluZ,0.4), xlab="Percent Identity", ylab="SCI", xlim=c(10,90),ylim=c(0.05,1.35),type="l",col=2,lty=1,main="", font=2)

lines(lowess(id,structZ,0.4),col=1,lty=1,lwd=2)
lines(lowess(id,musZ,0.4),col=4,lty=1)
lines(lowess(id,pcmaZ,0.4),col=3,lty=4,lwd=1)
lines(lowess(id,poaZ,0.4),col=8,lty=4)
lines(lowess(id,proalignZ,0.4),col=12,lty=4)
lines(lowess(id,prrnZ,0.4),col=6,lty=1)
lines(lowess(id,folda2globalZ,0.4),col=3,lty=1,lwd=3)
lines(lowess(id,pcomZ,0.4),col=5,lty=1,lwd=3)
lines(lowess(id,pcomfZ,0.4),col=5,lty=2,lwd=3)

lines(lowess(id,dynalignZ,0.4),col=7,lty=1,lwd=3)
lines(lowess(id,stemloc_slowZ,0.4),col=8,lty=1,lwd=3)
lines(lowess(id,stemloc_fastZ,0.4),col=8,lty=2,lwd=3)

legend(60,0.75,c("Reference","ClustalW", "MUSCLE", "PCMA", "POA (gp)", "ProAlign", "Prrn","Dynalign", "Foldalign2.0", "PMcomp", "PMcomp (fast)", "Stemloc (slow)", "Stemloc (fast)","Stral"), col = c(1,2,4,3,8,12,6,7,3,5,5,8,8,7),lty=c(1,1,1,4,4,4,1,1,1,1,2,1,2,2),lwd=c(2,1,1,1,1,1,1,3,3,3,3,3,3,3),ncol=1)

dev.copy2eps(file="tRNA_pw_rnaz_6.eps", height=10)

######################################################################

inp<-scan("ave_pw_seq_id_noaln12_sorted.dat",list(n="",z=0.0))
id<-inp$z

inp<-scan("clustal_bali_noaln12_sorted.dat",list(n="",z=0.0))
cluB<-inp$z
inp<-scan("foldalign_bali_noaln12_sorted.dat",list(n="",z=0.0))
foldaB<-inp$z

inp<-scan("foldalign.global.fold_bali_noaln12_sorted.dat",list(n="",z=0.0))
folda2globalB<-inp$z

inp<-scan("muscle_bali_noaln12_sorted.dat",list(n="",z=0.0))
musB<-inp$z
inp<-scan("pmcomp-fast_bali_noaln12_sorted.dat",list(n="",z=0.0))
pcomfB<-inp$z
inp<-scan("pmcomp_bali_noaln12_sorted.dat",list(n="",z=0.0))
pcomB<-inp$z

inp<-scan("pcma_bali_noaln12_sorted.dat",list(n="",z=0.0))
pcmaB<-inp$z

inp<-scan("poa_bali_noaln12_sorted.dat",list(n="",z=0.0))
poaB<-inp$z

inp<-scan("proalign_bali_noaln12_sorted.dat",list(n="",z=0.0))
proalignB<-inp$z

inp<-scan("prrn_bali_noaln12_sorted.dat",list(n="",z=0.0))
prrnB<-inp$z

inp<-scan("stemloc-slow_bali_noaln12_sorted.dat",list(n="",z=0.0))
stemloc_slowB<-inp$z
inp<-scan("stemloc-fast110_bali_noaln12_sorted.dat",list(n="",z=0.0))
stemloc_fastB<-inp$z

inp<-scan("dynalign_bali_noaln12_sorted.dat",list(n="",z=0.0))
dynalignB<-inp$z

X11()
plot(lowess(id,cluB,0.4), xlab="Percent Identity", ylab="SPS", xlim=c(10,90),ylim=c(0.0,1.0),type="l",col=2,lty=1,main="", font = 2)

structB<-seq(1,1,length=118)
lines(lowess(id,structB,0.4),col=1,lty=1,lwd=2)

lines(lowess(id,musB,0.4),col=4,lty=1)

lines(lowess(id,pcmaB,0.4),col=3,lty=4,lwd=1)

lines(lowess(id,poaB,0.4),col=8,lty=4)
lines(lowess(id,proalignB,0.4),col=12,lty=4)
lines(lowess(id,prrnB,0.4),col=6,lty=1)

lines(lowess(id,folda2globalB,0.4),col=3,lty=1,lwd=3)

lines(lowess(id,pcomB,0.4),col=5,lty=1,lwd=3)
lines(lowess(id,pcomfB,0.4),col=5,lty=2,lwd=3)

lines(lowess(id,dynalignB,0.4),col=7,lty=1,lwd=3)

lines(lowess(id,stemloc_slowB,0.4),col=8,lty=1,lwd=3)
lines(lowess(id,stemloc_fastB,0.4),col=8,lty=2,lwd=3)

legend(60,0.55,c("Reference","ClustalW", "MUSCLE", "PCMA", "POA (gp)", "Proalign", "Prrn", "Dynalign", "Foldalign2.0", "PMcomp", "PMcomp (fast)",  "Stemloc (slow)", "Stemloc (fast)", "Stral"), col = c(1,2,4,3,8,12,6,7,3,5,5,8,8,7),lty=c(1,1,1,4,4,4,1,1,1,1,2,1,2,2),lwd=c(2,1,1,1,1,1,1,3,3,3,3,3,3,3),ncol=1)

dev.copy2eps(file="tRNA_pw_bali_6.eps", height=10)

######################################################################

