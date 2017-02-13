#!/usr/bin/env python

import os, string, sys
from useful import unzip

def compacthits(hitpositions):
    #Sla ihop overlappande traffar
    if len(hitpositions)==0:
        return []
    hitpositions.sort()
    #Hits on the negative strand will be given negative coordinates
    hits = []
    hits.append(hitpositions[0])
    for i in range(1,len(hitpositions)):
        if hits[-1][1]>=hitpositions[i][0]:
            hits[-1][1] = max(hits[-1][1],hitpositions[i][1])
        else:
            hits.append(hitpositions[i])
    return hits

def optimise(targetfile,hitfile,grep,tgrep,genomelength=80000002,cutoffs=[],tpfp=0,old=0):
    #tpfp indicate if (tp,fp,tn,fn) should be the output or if (the
    #default) optimal (cutoff,sens,spec,mcc) should be the output.
    #start, stop positions and sign
    def hitsc(x):
        #Function for reading positions and scores from a line
        spl = x.split()
        pos = [int(spl[6]+spl[3]),int(spl[6]+spl[4])]
        pos.sort()
        return (pos,float(spl[5]))
    def readtarget(x):
        spl = x.split()
        pos = [int(spl[6]+spl[3]),int(spl[6]+spl[4])]
        pos.sort()
        return pos
    def hit(X,Y):
        d = min(X[1],Y[1]) - max(X[0],Y[0]) + 1
        if d <= 0:
            return 0
        return d
    targets = map(lambda x: readtarget(x),
                  os.popen("cat %s | egrep %s" % (targetfile,tgrep)))
    targets = compacthits(targets)
    (hitpositions,scores) = unzip(map(lambda x:
                                      hitsc(x),
                                      os.popen("cat %s | egrep %s | sort -n -k 6 -r" %
                                               (hitfile,grep)).readlines()))
    if hitfile.find("SAM35")>-1:#Adjust the positions
        print "Warning!! Have you adjusted the positions by taking the number - 90000001?"
        print "If not, then run the following command before computing optimal MCC"
        print "cat %s|awk '{print $1,$2,$3,$4-90000001,$5-90000001,$6,$7}'|tr ' ' '\\t' > tmp; mv tmp %s" % (hitfile,hitfile)
        print "AND remember to adjust coordinates in the reversed genome!"
        scores  = map(lambda x: -x, scores)
        hitpositions.reverse()
        scores.reverse()
        cutoffs = map(lambda x: -x, cutoffs)
    if len(cutoffs)==0:
        cutoffs = Unique_sorted(scores)
        if old==1:
            cutoffs.sort()
    targetslength = reduce(lambda y,x: y + x[1]-x[0]+1, targets, 0)
    tn=[];fn=[];tp=[];fp=[];mcc=[]
    hits = []
    index0 = 0
    for cutoff in cutoffs:
        index = scores.index(cutoff)
        hits = compacthits(hitpositions[index0:index]+hits)
        index0 = index
        hitslength = reduce(lambda y,x: y + x[1]-x[0]+1, hits, 0)
        if len(hits)==0 or len(targets)==0:
            tp.append(0)
            fp.append(hitslength)
            tn.append(genomelength-targetslength)
            fn.append(targetslength)
        else:
            t = 0
            for T in targets:
                for H in hits:
                    h = hit(T,H)
                    t += h
            tp.append(t)
            fn.append(targetslength - t)
            fp.append(hitslength - t)
            tn.append(genomelength - t - fp[-1] -fn[-1])
        try:
            mcc.append((tp[-1]*tn[-1]-fp[-1]*fn[-1]) /
                       ((tp[-1]+fp[-1])*(tp[-1]+fn[-1]) * 
                        (tn[-1]+fp[-1])*(tn[-1]+fn[-1]))**0.5) 
        except ZeroDivisionError:
            #print "ZeroDiv", tp[-1],tn[-1],fp[-1],fn[-1]
            mcc.append(0)
    if len(mcc)==0:
        print "No hits for %s" % grep
        return (0,0,0,0)
    i = mcc.index(max(mcc))
    sens = float(tp[i]) / (tp[i] + fn[i])
    spec = float(tn[i]) / (tn[i] + fp[i])
    print "%g\t%g\t%g" % (sens,spec,mcc[i])
    if tpfp==1:
        return (cutoffs,tp,fp,tn,fn)
    elif tpfp==2:
        return (cutoffs[i],tp[i],fp[i],tn[i],fn[i])
    return (cutoffs[i],sens,spec,mcc[i])

def optimiseall(hitfile,target="target",outfile=""):
    if len(outfile)==0:
        outfile = hitfile.replace(".gff","_%s.opt" % target)
    targetfile = "ncRNA.gff"
    X=[]
    for (id,RNA) in map(lambda x: x.split(".")[-2].split("/")[-2:],
                        os.popen("ls -1 ../*/id*/*_t*.fa").readlines()):
        grep = "%s|grep %s" % (RNA,id)
        if target == "target":
            tgrep = RNA
        elif target=="all":
            tgrep = RNA.split("_")[0]
        elif target=="all_np":
            tgrep = RNA.split("_")[0]+"|grep -v pseudo"
        X.append(((id,RNA),optimise(targetfile,hitfile,grep,tgrep)))
    file=open(outfile,'w')
    for ((id,RNA),(cutoff,sens,spec,mcc)) in X:
        file.write("%s\t%s\t%g\t%g\t%g\t%g\n" % (id,RNA,cutoff,sens,spec,mcc))
    file.close()
    return X

def Unique_sorted(X):
    uniqueX = [X[0]]
    for x in X:
        if x != uniqueX[-1]:
            uniqueX.append(x)
    return uniqueX

if __name__ == "__main__":

    def usage():
        sys.stderr.write("""
        Error: Usage:
        genome_optMCC.py gff-file target
        target is one of "target", "all", "w_pseudo"
        \n""")
        raise SystemExit
    
    if len(sys.argv)==6:
        cutoff = float(sys.argv[5])
    else:
        cutoff = 0.0
    if len(sys.argv)>=5:
        target = sys.argv[4]
    elif len(sys.argv)<2:
        usage()
    elif len(sys.argv)==1:
        target = "target"
    else:
        target = sys.argv[2]
    hitfile = sys.argv[1]

    x=optimiseall(hitfile,target)
    
