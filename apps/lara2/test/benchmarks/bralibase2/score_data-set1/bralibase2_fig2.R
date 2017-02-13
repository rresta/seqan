inp<-scan("ALL_ave_perc_id_grouped_HML_recalculated_noSRP.dat", list(n="",id=0.0, h=""))
id<-inp$id
inp<-scan("ALL_clu_rnaz_noSRP.dat", list(n="", z=0.0))
cluZ<-inp$z
inp<-scan("ALL_cluqt_rnaz_noSRP.dat", list(n="", z=0.0))
cluqtZ<-inp$z
inp<-scan("ALL_dial_rnaz_noSRP.dat", list(n="", z=0.0))
dialZ<-inp$z
inp<-scan("ALL_dial_it_rnaz_noSRP.dat", list(n="", z=0.0))
dial_itZ<-inp$z
inp<-scan("ALL_dial_ito_rnaz_noSRP.dat", list(n="", z=0.0))
dial_itoZ<-inp$z
inp<-scan("ALL_dial_o_rnaz_noSRP.dat", list(n="", z=0.0))
dial_oZ<-inp$z
inp<-scan("ALL_k5_rnaz_noSRP.dat", list(n="", z=0.0))
k5Z<-inp$z
inp<-scan("ALL_muscle_rnaz_noSRP.dat", list(n="", z=0.0))
muscleZ<-inp$z
inp<-scan("ALL_muscle_mi32_rnaz_noSRP.dat", list(n="", z=0.0))
muscle_mi32Z<-inp$z
inp<-scan("ALL_muscle_mi32mt6_rnaz_noSRP.dat", list(n="", z=0.0))
muscle_mi32mt6Z<-inp$z
inp<-scan("ALL_muscle_mt6_rnaz_noSRP.dat", list(n="", z=0.0))
muscle_mt6Z<-inp$z
inp<-scan("ALL_muscle_nj_rnaz_noSRP.dat", list(n="", z=0.0))
muscle_njZ<-inp$z
inp<-scan("ALL_muscle_njmi32_rnaz_noSRP.dat", list(n="", z=0.0))
muscle_njmi32Z<-inp$z
inp<-scan("ALL_muscle_njmi32mt6_rnaz_noSRP.dat", list(n="", z=0.0))
muscle_njmi32mt6Z<-inp$z
inp<-scan("ALL_muscle_njmt6_rnaz_noSRP.dat", list(n="", z=0.0))
muscle_njmt6Z<-inp$z
inp<-scan("ALL_poa_rnaz_noSRP.dat", list(n="", z=0.0))
poaZ<-inp$z
inp<-scan("ALL_pcma_rnaz_noSRP.dat", list(n="", z=0.0))
pcmaZ<-inp$z
inp<-scan("ALL_poa_g_rnaz_noSRP.dat", list(n="", z=0.0))
poa_gZ<-inp$z
inp<-scan("ALL_poa_gp_rnaz_noSRP.dat", list(n="", z=0.0))
poa_gpZ<-inp$z
inp<-scan("ALL_poa_p_rnaz_noSRP.dat", list(n="", z=0.0))
poa_pZ<-inp$z
inp<-scan("ALL_proaln_bw400_rnaz_noSRP.dat", list(n="", z=0.0))
proalnZ<-inp$z
inp<-scan("ALL_prrn_rnaz_noSRP.dat", list(n="", z=0.0))
prrnZ<-inp$z
inp<-scan("ALL_prrn_S10_rnaz_noSRP.dat", list(n="", z=0.0))
prrn_S10Z<-inp$z
inp<-scan("ALL_tco_c_rnaz_noSRP.dat", list(n="", z=0.0))
tco_cZ<-inp$z
inp<-scan("ALL_tco_rnaz_noSRP.dat", list(n="", z=0.0))
tcoZ<-inp$z
inp<-scan("ALL_tco_f_rnaz_noSRP.dat", list(n="", z=0.0))
tco_fZ<-inp$z
inp<-scan("ALL_tco_s_rnaz_noSRP.dat", list(n="", z=0.0))
tco_sZ<-inp$z

inp<-scan("ALL_handle_rnaz_noSRP.dat", list(n="", z=0.0))
handleZ<-inp$z
inp<-scan("ALL_mafft_fftnsi_rnaz_noSRP.dat", list(n="", z=0.0))
mafft_fftnsiZ<-inp$z
inp<-scan("ALL_mafft_fftns_rnaz_noSRP.dat", list(n="", z=0.0))
mafft_fftnsZ<-inp$z
inp<-scan("ALL_mafft_nwnsi_rnaz_noSRP.dat", list(n="", z=0.0))
mafft_nwnsiZ<-inp$z
inp<-scan("ALL_mafft_nwns_rnaz_noSRP.dat", list(n="", z=0.0))
mafft_nwnsZ<-inp$z

inp<-scan("ALL_alignm_mrna2_p2mfmin07_p2mnseqmin5_rnaz_noSRP.dat", list(n="", z=0.0))
alignm_mrna2_p2mfmin07_p2mnseqmin5Z<-inp$z
inp<-scan("ALL_alignm_mrna2_rnaz_noSRP.dat", list(n="", z=0.0))
alignm_mrna2Z<-inp$z
inp<-scan("ALL_alignm_mrna2_s2pgo10_s2pge1_p2mfmin07_p2mnseqmin5_rnaz_noSRP.dat", list(n="", z=0.0))
alignm_mrna2_s2pgo10_s2pge1_p2mfmin07_p2mnseqmin5Z<-inp$z
inp<-scan("ALL_alignm_mrna2_s2pgo10_s2pge1_rnaz_noSRP.dat", list(n="", z=0.0))
alignm_mrna2_s2pgo10_s2pge1Z<-inp$z
inp<-scan("ALL_alignm_mrna2_s2pw3_rnaz_noSRP.dat", list(n="", z=0.0))
alignm_mrna2_s2pw3Z<-inp$z



X11()

plot(lowess(id, k5Z, f=0.7), xlab = "Percent Identity", ylab = "SCI",xlim=c(35,95),ylim=c(-0.1,1.3), type="l",lty=1,col=1,lwd=2,main="", font = 2)
lines(lowess(id, cluZ, f=0.7),col=2,lty=1)

lines(lowess(id, dialZ, f=0.7),col=3,lty=1,lwd=2)

lines(lowess(id, muscleZ, f=0.7),col=4,lty=1,lwd=2)

lines(lowess(id, pcmaZ, f=0.7),col=3,lty=4,lwd=2)

lines(lowess(id, poaZ, f=0.7),col=5,lty=1,lwd=2)
lines(lowess(id, poa_gpZ, f=0.7),col=5,lty=2,lwd=2)
lines(lowess(id, proalnZ, f=0.7),col=12,lty=4,lwd=2)
lines(lowess(id, prrnZ, f=0.7),col=6,lty=1,lwd=2)
lines(lowess(id, prrn_S10Z, f=0.7),col=6,lty=2,lwd=2)

lines(lowess(id, tcoZ, f=0.7),col=8,lty=1,lwd=2)
lines(lowess(id, tco_cZ, f=0.7),col=8,lty=2,lwd=2)
lines(lowess(id, handleZ, f=0.7),col=9,lty=2,lwd=2)
lines(lowess(id, mafft_nwnsiZ, f=0.7),col=7,lty=1,lwd=2)

lines(lowess(id, alignm_mrna2Z, f=0.7),col=10,lty=2,lwd=2)

legend(62.5, 0.37,c("Reference","ClustalW","DIALIGN","MUSCLE","PCMA","POA","POA (g,p)","ProAlign","Prrn","Prrn (S10)","T-Coffee","T-Coffee (c)","Handel","MAFFT", "Align-m (1)"), col = c(1,2,3,4,3,5,5,12,6,6,8,8,9,7,10),lty = c(1,1,1,1,4,1,2,4,1,2,1,2,2,1,2),lwd = c(2,2,2,2,2,2,2,2,2,2,2,2,2,2,2), ncol=2)

dev.copy2eps(file="ALL_rnaz_sci_noSRP_5.eps", height=10)
 

######################################################################

inp<-scan("ALL_ave_perc_id_grouped_HML_recalculated_noSRP.dat", list(n="",id=0.0, h=""))
id<-inp$id
inp<-scan("ALL_clu_baliscore_noSRP.dat", list(n="", z=0.0))
cluB<-inp$z
inp<-scan("ALL_cluqt_baliscore_noSRP.dat", list(n="", z=0.0))
cluqtB<-inp$z
inp<-scan("ALL_dial_baliscore_noSRP.dat", list(n="", z=0.0))
dialB<-inp$z
inp<-scan("ALL_dial_it_baliscore_noSRP.dat", list(n="", z=0.0))
dial_itB<-inp$z
inp<-scan("ALL_dial_ito_baliscore_noSRP.dat", list(n="", z=0.0))
dial_itoB<-inp$z
inp<-scan("ALL_dial_o_baliscore_noSRP.dat", list(n="", z=0.0))
dial_oB<-inp$z
inp<-scan("ALL_k5_baliscore_noSRP.dat", list(n="", z=0.0))
k5B<-inp$z
inp<-scan("ALL_muscle_baliscore_noSRP.dat", list(n="", z=0.0))
muscleB<-inp$z
inp<-scan("ALL_muscle_mi32_baliscore_noSRP.dat", list(n="", z=0.0))
muscle_mi32B<-inp$z
inp<-scan("ALL_muscle_mi32mt6_baliscore_noSRP.dat", list(n="", z=0.0))
muscle_mi32mt6B<-inp$z
inp<-scan("ALL_muscle_mt6_baliscore_noSRP.dat", list(n="", z=0.0))
muscle_mt6B<-inp$z
inp<-scan("ALL_muscle_nj_baliscore_noSRP.dat", list(n="", z=0.0))
muscle_njB<-inp$z
inp<-scan("ALL_muscle_njmi32_baliscore_noSRP.dat", list(n="", z=0.0))
muscle_njmi32B<-inp$z
inp<-scan("ALL_muscle_njmi32mt6_baliscore_noSRP.dat", list(n="", z=0.0))
muscle_njmi32mt6B<-inp$z
inp<-scan("ALL_muscle_njmt6_baliscore_noSRP.dat", list(n="", z=0.0))
muscle_njmt6B<-inp$z
inp<-scan("ALL_pcma_baliscore_noSRP.dat", list(n="", z=0.0))
pcmaB<-inp$z
inp<-scan("ALL_poa_baliscore_noSRP.dat", list(n="", z=0.0))
poaB<-inp$z
inp<-scan("ALL_poa_g_baliscore_noSRP.dat", list(n="", z=0.0))
poa_gB<-inp$z
inp<-scan("ALL_poa_gp_baliscore_noSRP.dat", list(n="", z=0.0))
poa_gpB<-inp$z
inp<-scan("ALL_poa_p_baliscore_noSRP.dat", list(n="", z=0.0))
poa_pB<-inp$z
inp<-scan("ALL_proaln_bw400_baliscore_noSRP.dat", list(n="", z=0.0))
proalnB<-inp$z
inp<-scan("ALL_prrn_baliscore_noSRP.dat", list(n="", z=0.0))
prrnB<-inp$z
inp<-scan("ALL_prrn_S10_baliscore_noSRP.dat", list(n="", z=0.0))
prrn_S10B<-inp$z
inp<-scan("ALL_tco_c_baliscore_noSRP.dat", list(n="", z=0.0))
tco_cB<-inp$z
inp<-scan("ALL_tco_baliscore_noSRP.dat", list(n="", z=0.0))
tcoB<-inp$z
inp<-scan("ALL_tco_f_baliscore_noSRP.dat", list(n="", z=0.0))
tco_fB<-inp$z
inp<-scan("ALL_tco_s_baliscore_noSRP.dat", list(n="", z=0.0))
tco_sB<-inp$z

inp<-scan("ALL_handle_baliscore_noSRP.dat", list(n="", z=0.0))
handleB<-inp$z
inp<-scan("ALL_mafft_fftnsi_baliscore_noSRP.dat", list(n="", z=0.0))
mafft_fftnsiB<-inp$z
inp<-scan("ALL_mafft_fftns_baliscore_noSRP.dat", list(n="", z=0.0))
mafft_fftnsB<-inp$z
inp<-scan("ALL_mafft_nwnsi_baliscore_noSRP.dat", list(n="", z=0.0))
mafft_nwnsiB<-inp$z
inp<-scan("ALL_mafft_nwns_baliscore_noSRP.dat", list(n="", z=0.0))
mafft_nwnsB<-inp$z

inp<-scan("ALL_alignm_mrna2_p2mfmin07_p2mnseqmin5_baliscore_noSRP.dat", list(n="", z=0.0))
alignm_mrna2_p2mfmin07_p2mnseqmin5B<-inp$z
inp<-scan("ALL_alignm_mrna2_baliscore_noSRP.dat", list(n="", z=0.0))
alignm_mrna2B<-inp$z
inp<-scan("ALL_alignm_mrna2_s2pgo10_s2pge1_p2mfmin07_p2mnseqmin5_baliscore_noSRP.dat", list(n="", z=0.0))
alignm_mrna2_s2pgo10_s2pge1_p2mfmin07_p2mnseqmin5B<-inp$z
inp<-scan("ALL_alignm_mrna2_s2pgo10_s2pge1_baliscore_noSRP.dat", list(n="", z=0.0))
alignm_mrna2_s2pgo10_s2pge1B<-inp$z
inp<-scan("ALL_alignm_mrna2_s2pw3_baliscore_noSRP.dat", list(n="", z=0.0))
alignm_mrna2_s2pw3B<-inp$z


X11()
plot(lowess(id, cluB, f=0.5), xlab = "Percent Identity", ylab = "SPS",xlim=c(35,95), ylim=c(0.18,1.0), type="l",lty=1,col=2,lwd=2,main="", font = 2)
lines(c(35,95),c(1,1),col=1,lty=1,lwd=2)

lines(lowess(id, dialB, f=0.5),col=3,lty=1,lwd=2)
lines(lowess(id, muscleB, f=0.5),col=4,lty=1,lwd=2)

lines(lowess(id, pcmaB, f=0.5),col=3,lty=4,lwd=2)

lines(lowess(id, poaB, f=0.5),col=5,lty=1,lwd=2)
lines(lowess(id, poa_gpB, f=0.5),col=5,lty=2,lwd=2)

lines(lowess(id, proalnB, f=0.5),col=12,lty=4,lwd=2)

lines(lowess(id, prrnB, f=0.5),col=6,lty=1,lwd=2)
lines(lowess(id, prrn_S10B, f=0.5),col=6,lty=2,lwd=2)

lines(lowess(id, tcoB, f=0.5),col=8,lty=1,lwd=2)
lines(lowess(id, tco_cB, f=0.5),col=8,lty=2,lwd=2)

lines(lowess(id, handleB, f=0.5),col=9,lty=2,lwd=2)

lines(lowess(id, mafft_nwnsiB, f=0.5),col=7,lty=1,lwd=2)

lines(lowess(id, alignm_mrna2B, f=0.5),col=10,lty=2,lwd=2)

legend(62, 0.47,c("Reference","ClustalW","DIALIGN","MUSCLE","PCMA","POA","POA (g,p)","ProAlign","Prrn","Prrn (S10)","T-Coffee","T-Coffee (c)","Handel","MAFFT", "Align-m (1)"), col = c(1,2,3,4,3,5,5,12,6,6,8,8,9,7,10),lty = c(1,1,1,1,4,1,2,4,1,2,1,2,2,1,2),lwd = c(2,2,2,2,2,2,2,2,2,2,2,2,2,2,2),ncol=2)

dev.copy2eps(file="ALL_baliscore_sci_noSRP_5.eps", height = 10)
 
