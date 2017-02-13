#!/usr/bin/env python
# -*- coding: UTF-8 -*-
import os, tempfile, string, sys, pickle, time, math

def wublast(fastafile,database,opt=""):
    """
    Run WU-blast on all the sequences in fastafile against database
    using the options defined in opt, the default is to use the
    program's default options.

    The output is a set of hits, the sequence identifier and the
    maximum score correlated to this sequence.
    """
    (id,seq) = readonelinefasta(fastafile)
    hits = []
    scoreindex = 1 # The score column in the result file
    #Run one query sequence against the database at the time.
    for k in range(len(id)):
        hits.append([])
        tmpfilename=tempfile.mktemp()
        tmpfile=open(tmpfilename,'w')
        tmpfile.write('>'+id[k]+'\n'+\
                      seq[k])
        tmpfile.close()
        lines = os.popen('blastn %s %s -noseqs %s B=9999999 V=9999999 S=1 S2=1' %
                         (database,tmpfilename,opt)).\
                         readlines()
        os.system("rm %s" % tmpfilename)
        #Find the list of hits
        try:
            i = lines.index("Sequences producing High-scoring Segment Pairs:              Score  P(N)      N\n")
            j = i+2+lines[i+2:].index("\n")
            hits[k] += map(lambda x: (x.split()[0],
                                      float(x.split()[scoreindex])),
                        lines[i+2:j])
        except ValueError:
            sys.stderr.write("No hits\n")
        #Find unique hits, if there are several hits to the same sequence,
        # save only the highest scoring one.
        hits[k] = Unique(hits[k],max)
    return hits
    
def ncbiblast(fastafile,database,opt=""):
    """
    Run NCBI-blast on all the sequences in fastafile against database
    using the options defined in opt, the default is to use the
    program's default options.

    The output is a set of hits, the sequence identifier and the
    maximum score correlated to this sequence.
    """
    (id,seq) = readonelinefasta(fastafile)
    hits = []
    scoreindex = 1 # The score column
    #Run blast on one query sequence at the time.
    for k in range(len(id)):
        hits.append([])
        tmpfilename=tempfile.mktemp()
        tmpfile=open(tmpfilename,'w')
        tmpfile.write('>'+id[k]+'\n'+\
                      seq[k])
        tmpfile.close()
        lines = os.popen('blastall -p blastn -d %s -i %s %s' %
                         (database,tmpfilename,opt) +
                         ' -e 1e10000 -b 0 -v 1000000').readlines()
        os.system('rm '+tmpfilename)
        #Find the list of hits
        try:
            i = lines.index("Sequences producing significant alignments:                      (bits) Value\n")
            j = i+2+lines[i+2:].index("\n")
            try:
                hits[k] += map(lambda x: (x.split()[0],
                                          float(x.split()[scoreindex])),
                               lines[i+2:j])
            except IndexError:
                hits[k] += map(lambda x: (x.split()[0],
                                          float(x.split()[scoreindex])),
                               lines[i+2:j-1])
        except ValueError:
            if len(filter(lambda x:
                          x.find("***** No hits found ******")>-1,lines))>0:
                sys.stderr.write("No hits!\n")
            else:
                sys.stderr.write(string.join(lines))
                raise ValueError
        hits[k] = Unique(hits[k],max)
    return hits

def fasta34(fastafile,database,opt=""):
    """
    Use fasta 3.4 to identify homologs to the sequences in fastafile
    in database using the options defined in opt, the default is to
    use the program's default options.

    The output is a set of hits, the sequence identifier and the
    maximum score correlated to this sequence.
    """
    hits = []
    opts = opt.split()
    if len(opts) % 2 == 1: #odd
        ktup = opts[-1]
        if ktup[0]=="-":#Not ktup
            ktup = ""
        else:
            opt = string.join(opts[:-1]," ")
    else:
        ktup = ""
    lines = os.popen('fasta34 -H -b 9999999 -n -d 0 -q %s %s %s %s' %
                     (opt,fastafile,database,ktup)).readlines()
    N = len(readonelinefasta(database)[0])
    hits = []
    k = 0
    scoreindex = -3 #Score column
    while 1:
        i=0
        try:
            while lines[i].find("The best scores are:")==-1:
                i += 1
        except IndexError:
            break
        hits.append([])
        lines = lines[i+1:]
        j = lines.index("\n")
        hits[k] = map(lambda x: (x.split()[0], float(x.split()[scoreindex])),
                      lines[:j])
        hits[k] = Unique(hits[k],max)
        k += 1
    return hits

def ssearch34(fastafile,database,opt=""):
    """
    Use ssearch 3.4 to identify homologs to the sequences in fastafile
    in database using the options defined in opt, the default is to
    use the program's default options.

    The output is a set of hits, the sequence identifier and the
    maximum score correlated to this sequence.
    """
    hits = []
    opts = opt.split()
    if len(opts) % 2 == 1: #odd
        ktup = opts[-1]
        if ktup[0]=="-":
            #Not ktup
            ktup = ""
        else:
            opt = string.join(opts[:-1]," ")
            print "ktup = " + ktup
    else:
        ktup = ""
    lines = os.popen('ssearch34 -H -b 9999999 -n -d 0 -q %s %s %s %s' %
                     (opt,fastafile,database,ktup)).readlines()
    N = len(readonelinefasta(database)[0])
    hits = []
    k = 0
    while 1:
        i=0
        try:
            while lines[i].find("The best scores are:")==-1:
                i += 1
        except IndexError:
            break
        hits.append([])
        lines = lines[i+1:]
        j = lines.index("\n")
        hits[k] += map(lambda x: (x.split()[0], float(x.split()[-3])), lines[:j])
        hits[k] = Unique(hits[k],max)
        k += 1
    return hits

def hmmer2(fastafile,database,opt=""):
    """
    Use HMMer 2.3.2 to identify homologs to the sequences in fastafile
    in database using the options defined in opt, the default is to
    use the program's default options.

    The output is a set of hits, the sequence identifier and the
    maximum score correlated to this sequence.
    """
    name = fastafile.split(".")
    if len(name)>1:
        name = string.join(name[:-1],".")
    else:
        name = name[0]
    if not os.path.exists(name+"_pro.aln"):
        os.system("java -jar ProAlign_0.5a1.jar -nogui  " +
                  "-seqfile=%s -newtree -outfile=%s.pir -outformat=pir" %
                  (fastafile,name))
        (id,seq) = readfasta(name+".pir")
        writeclustal(name+"_pro.aln",map(lambda x: x.split("DL;")[-1],id),
                     map(lambda x: x.replace("*",""),seq))
    os.system("hmmbuild %s -F -n %s %s.hmm2 %s_pro.aln" %
              (opt,name,name,name))
    lines = os.popen("hmmsearch -E 99999999 %s.hmm2 %s" % (name,database)).readlines()
    i=0
    while lines[i].find("Sequence")==-1 or \
              lines[i].find("Description")==-1 or \
              lines[i].find("Score")==-1 or \
              lines[i].find("E-value")==-1:
        i += 1
    lines = lines[i+2:]
    j = lines.index("\n")
    hits = map(lambda x: (x.split()[0],float(x.split()[1])), lines[:j])
    return [Unique(hits,max)]

def hmmer1(fastafile,database,opt=""):
    """
    Use HMMer 1.8.4 to identify homologs to the sequences in fastafile
    in database using the options defined in opt, the default is to
    use the program's default options.

    The output is a set of hits, the sequence identifier and the
    maximum score correlated to this sequence.
    
    Use hmmls for global/local, multidomain
    Use hmms for global/local
    """
    name = fastafile.split(".")
    if len(name)>1:
        name = string.join(name[:-1],".")
    else:
        name = name[0]

    if type(opt)==type(""):
        buildmethod = "hmmt"
        searchmethod = "hmms"
        while 1:
            if len(opt)>3 and opt[:3]=="hmm":
                sp = opt.split()
                opt = string.join(sp[1:])
                if sp[0][3]=="t" or sp[0][3]=="b":
                    buildmethod = sp[0]
                    nexthmm = opt.find("hmm")
                    if nexthmm>-1:
                        buildmethod += " " + opt[:nexthmm]
                        opt = opt[nexthmm:]
                    else:
                        buildmethod += " " + opt
                        opt=""
                        break
                else:
                    searchmethod = sp[0]
                    nexthmm = opt.find("hmm")
                    if nexthmm>-1:
                        searchmethod += " " + opt[:nexthmm]
                        opt = opt[nexthmm:]
                    else:
                        searchmethod += " " + opt
                        opt = ""
                        break
            else:
                break
        os.system("%s %s %s.hmm1 %s" %
                  (buildmethod,opt,name,fastafile))
        lines = os.popen("%s -t -99999 %s.hmm1 %s" %
                         (searchmethod,name,database)).readlines()
    else:
        os.system("hmmt %s.hmm1 %s" % (name,fastafile))
        lines = os.popen("hmms -t -99999 %s.hmm1 %s" % (name,database)).readlines()
    i=-1
    while i<len(lines):
        i = i+1+lines[i+1:].index(" - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n")
        if lines[i+1]=="\n":
            break
    lines = lines[i+2:]
    hits = map(lambda x: (x.split()[-1],float(x.split()[0])), lines)
    return [Unique(hits,max)]

def sam35(fastafile,database,opt=""):
    """
    Use SAM 3.5 to identify homologs to the sequences in fastafile
    in database using the options defined in opt, the default is to
    use the program's default options.

    The output is a set of hits, the sequence identifier and the
    maximum score correlated to this sequence.
    """
    name = fastafile.split(".")
    if len(name)>1:
        name = string.join(name[:-1],".")
    else:
        name = name[0]
    name2 = name + "_" + database.split("/")[-1]
    os.system("buildmodel %s -a RNA -train %s " % (name,fastafile))
    os.system("hmmscore %s %s -i %s.mod -db %s" % (name2,opt,name,database))
    lines = open("%s.dist" % name2,"r").readlines()
    i=0
    while i<len(lines):
        if lines[i].find("Sequence ID")>-1:
            break
        i += 1
    lines = lines[i+1:]
    hits = map(lambda x: (x.split()[0],float(x.split()[3])), lines)
    return [Unique(hits,min)]

def samhmmer2(fastafile,database,opt=""):
    """
    Use SAM 3.5 to create the HMM model from the sequences in
    fastafile and use HMMer 2.3.2 to search the database using the
    cutoff cutoff and the options defined in opt. The default cutoff
    is 0 and the default is to use the program's default options.

    The output is a set of hits, the sequence identifier and the
    maximum score correlated to this sequence.
    """
    name = fastafile.split(".")
    if len(name)>1:
        name = string.join(name[:-1],".")
    else:
        name = name[0]
    os.system("buildmodel %s -a RNA -train %s " % (name,fastafile))
    os.system("sam2hmmer %s -i %s.mod; mv %s.hmmer %s.hmmersam35" %
              (name, name, name, name))
    lines = os.popen("hmmsearch -E 99999999 %s.hmmersam35 %s" % (name,database)).readlines()
    i=0
    while lines[i].find("Sequence")==-1 or \
              lines[i].find("Description")==-1 or \
              lines[i].find("Score")==-1 or \
              lines[i].find("E-value")==-1:
        i += 1
    lines = lines[i+2:]
    j = lines.index("\n")
    hits = map(lambda x: (x.split()[0],float(x.split()[1])), lines[:j])
    return [Unique(hits,max)]

def rsmatch12(fastafile,database,cutoff=0,opt=""):
    """
    Use RSmatch 1.2 to identify homologs to the sequences in fastafile
    in database using the options defined in opt, the default is to
    use the program's default options.

    The output is a set of hits, the sequence identifier and the
    maximum score correlated to this sequence.
    """
    #Save results so that they can be reused
    if opt!=0:
        optname = string.join(opt.split(),"_")
    else:
        optname = ""
    filename = fastafile + "_" + database.split("/")[-1] + ".RSmatch12"
    if os.path.exists(filename):
        sys.stderr.write("Load %s..." % filename)
        (name,score) = pickle.load(open(filename,"r"))
    else:
        #Make sure that both query and database are folded.
        structfile = fastafile.replace('.fa','.struct')
	if not os.path.exists(structfile):
	    os.system("RNAfold < %s | awk '{print $1}' > %s" %
		      (fastafile,structfile))
	if not os.path.exists("%s.db.struct" % database):
	    os.system("RNAfold < %s | awk '{print $1}' > %s.db.struct" %
		      (database,database))
        (id,seq,str) = readseqstr(structfile)
        name = []
        score=[]
        for k in range(len(id)):
            tmpfilename=tempfile.mktemp()
            tmpfile=open(tmpfilename,'w')
            tmpfile.write('>'+id[k]+' \n'+\
                          seq[k] + "\n" + str[k])
            tmpfile.close()
            resultfile = "result_%u.out" % time.time()
            os.system("RSmatch1.2 -p dsearch -d %s.db.struct -q %s -n 30000 -o %s" %
                      (database,tmpfilename,resultfile))
            file = open(resultfile, "r")
            line = file.readline()
            while (len(line)>0 and line.find("Hits")==-1):
                line = file.readline()
            line = file.readline()
            line = file.readline()
            while len(line)>5:
                try:
                    sp = line.split()
                    score.append(int(sp[1]))
                    sp = sp[3].split(":")
                    name.append(sp[0])
                    line = file.readline()
                except:
                   sys.stderr.write( "No hits?" )
            file.close()
            os.system("rm %s" % resultfile)
            os.system("rm %s*" % tmpfilename)

        pickle.dump((name,score),open(filename,"w"))

    hits = zip(name,score)

    return [Unique(hits,max)]


def paralign(fastafile,database,opt=""):
    """
    Use paralign to identify homologs to the sequences in fastafile
    in database using the options defined in opt, the default is to
    use the program's default options.

    The output is a set of hits, the sequence identifier and the
    maximum score correlated to this sequence.
    """
    (id,seq) = readonelinefasta(fastafile)
    if opt == 0:
	opt = ""
    hits = []
    for k in range(len(id)):
        hits.append([])
        tmpfilename=tempfile.mktemp()
        tmpfile=open(tmpfilename,'w')
        tmpfile.write('>'+id[k]+'\n'+\
                      seq[k])
        tmpfile.close()
        lines = os.popen('paralign -g sencel.lic %s -l 9999999999 -n %s %s' %
                         (opt,database,tmpfilename)).\
                         readlines()
        os.system('rm '+tmpfilename)
        #Find the list of hits
        try:
            i=0
	    while i<len(lines) and (lines[i].find("SEQUENCE")==-1 or lines[i].find("STRAND SCORE")==-1):
		i += 1
            j = lines[i+1:].index("\n")
            hits[k] += map(lambda x: (x.split()[1],float(x.split()[3])), lines[i+1:i+j+1])
        except ValueError:
            sys.stderr.write("No hits\n")
        hits[k] = Unique(hits[k],max)
    return hits


#Useful functions of different kinds
    
def Unique(hits,minmax = min):
    """
    Return unique hits with minimal or maximal scores.
    """
    uniquename = []
    uniqueE = []
    for x in hits:
        if x[0] not in uniquename:
            uniquename.append(x[0])
            uniqueE.append(x[1])
        else:
            i = uniquename.index(x[0])
            uniqueE[i] = minmax(uniqueE[i],x[1])
    return zip(uniquename,uniqueE)

def Uniquenames(hits):
    uniquehits = []
    for x in hits:
        if x not in uniquehits:
            uniquehits.append(x)
    return uniquehits

def readonelinefasta(fastafile):
    lines = open(fastafile,"r").readlines()
    id = map(lambda x: x[1:].strip(), filter(lambda x: x.find(">")==0, lines))
    seq = map(lambda x: x.strip(), filter(lambda x: x.find(">")!=0, lines))
    return (id,seq)

def readfasta(fastafile):
    lines = open(fastafile,"r").readlines()
    id=[]
    seq=[]
    k=0
    while k<len(lines):
        if lines[k][0]==">":
            id.append(lines[k][1:].strip())
            seq.append(lines[k+1].strip())
            k+=2
        else:
            seq[-1]+=lines[k].strip()
            k+=1
    return (id,seq)

def writeclustal(filename,id,seq):
    """
    Write the list of identities, id, and sequences, seq, in
    clustal format to the file filename.
    """
    file=open(filename,'w')
    file.write('CLUSTAL W (1.7) multiple sequence alignment\n\n\n')
    maxidlen = reduce(lambda y,x:max(len(x),y),id,0)
    seqlen = len(seq[0])
    for k in range(0,seqlen,60):
        for i in range(len(id)):
            n = maxidlen - len(id[i])
            file.write(id[i]+' '*(n+2)+seq[i][k:k+60]+'\n')
        file.write("\n\n")
    file.close()

def printhitscores(hits):
    n = len(hits)
    for h in hits:
        h.sort()
    all_hits = Uniquenames(reduce(lambda y,x: y+map(lambda z: z[0],x),hits,[]))
    all_hits.sort()
    hitsdict = {}
    for name in all_hits:
        hitsdict[name] = n*[-1] #-1 indicates no score at all
    for k in range(n):
        for hit in hits[k]:
            hitsdict[hit[0]][k] = hit[1]
            
    for name in all_hits:
        for i in range(n):
            sys.stdout.write("%s\t%g\t" % (name,hitsdict[name][i]))
        print
    
if __name__ == "__main__":

    def usage():
        sys.stderr.write("""
        Error: Usage:
        scores.py searchmethod fastafile database
        
        searchmethod is one of ncbiblast, wublast, fasta34, ssearch34,
        hmmer2, hmmer1, sam35, samhmmer2, rsmatch12 or paralign

        Runs the selected search method with fastafile as quesry and
        database as database. The output is a list of scores produced
        by the search algorithm.
        \n""")
        raise SystemExit
    
    if len(sys.argv) <4:
        usage()
                
    method = sys.argv[1]
    opt = ""
    if len(sys.argv)>=5:
        # specified optimisation, an option to send to
        #the search program
        opt = string.join(sys.argv[4:])
    if method=="ncbiblast" or method=="wublast" or \
           method =="fasta34" or method == "hmmer2" or \
           method == "hmmer1" or method=="samhmmer2" or\
           method == "rsmatch12" or\
           method == "paralign" or method == "sam35" or\
           method == "ssearch34":
        fastafile = sys.argv[2]
        database = sys.argv[3]
        exec("hits = " + method + "(fastafile,database,opt=opt)")
        #False positives
        fp = map(lambda y: filter(lambda x:x[0].find("_")>-1,y),hits)
        #True positives
        tp = map(lambda y: filter(lambda x:x[0].find("_")==-1,y),hits)
        printhitscores(tp)
        printhitscores(fp)
    else:
        usage()
