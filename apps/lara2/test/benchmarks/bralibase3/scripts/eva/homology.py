#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import os, tempfile, string, sys, pickle, time, math

def wublast(fastafile,database,cutoff=0,opt=""):
    """
    Run WU-blast on all the sequences in fastafile against database
    using the cutoff cutoff and the options defined in opt. The
    default cutoff is 0 and the default is to use the program's
    default options.

    The output is a list of hits.
    """
    (id,seq) = readonelinefasta(fastafile)
    hits = []
    scoreindex = 1 # The score column in the result file
    #Run one query sequence against the database at the time.
    for k in range(len(id)):
        tmpfilename=tempfile.mktemp()
        tmpfile=open(tmpfilename,'w')
        tmpfile.write('>'+id[k]+'\n'+\
                      seq[k])
        tmpfile.close()
        lines = os.popen('blastn %s %s -noseqs %s B=9999999 V=9999999 S=%d S2=1' %
                         (database,tmpfilename,opt,cutoff)).\
                         readlines()
        os.system("rm %s" % tmpfilename)
        #Find the list of hits
        try:
            i = lines.index("Sequences producing High-scoring Segment Pairs:              Score  P(N)      N\n")
            j = i+2+lines[i+2:].index("\n")                
            for line in lines[i+2:j]:
                if float(line.split()[scoreindex])>=cutoff:
                    hits.append(line.split()[0])
        except ValueError:
            sys.stderr.write("No hits or database incorrectly formatted (did you run pressdb?)\n")
    return Uniquenames(hits)
    
def ncbiblast(fastafile,database,cutoff=0,opt=""):
    """
    Run NCBI-blast on all the sequences in fastafile against database
    using the cutoff cutoff and the options defined in opt. The
    default cutoff is 0 and the default is to use the program's
    default options.

    The output is a list of hits.
    """
    (id,seq) = readonelinefasta(fastafile)
    hits = []
    scoreindex = 1 # The score column
    #Run blast on one query sequence at the time.
    for k in range(len(id)):
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
            for line in lines[i+2:j]:
                try:
                    if float(line.split()[scoreindex])>=cutoff:
                        hits.append(line.split()[0])
                except IndexError:
                    pass
        except ValueError:
            if len(filter(lambda x:
                          x.find("***** No hits found ******")>-1,lines))>0:
                sys.stderr.write("No hits!\n")
            else:
                sys.stderr.write("Is the database correctly formatted? If not run formatdb\n")
                raise ValueError
    return Uniquenames(hits)

def fasta34(fastafile,database,cutoff=0,opt=""):
    """
    Use fasta 3.4 to identify homologs to the sequences in fastafile
    in database using the cutoff cutoff and the options defined in
    opt. The default cutoff is 0 and the default is to use the
    program's default options.

    The output is a list of hits.
    """
    hits = []
    opts = opt.split()
    if len(opts)>0 and opts[-1][0]!="-":# ktup?
        ktup = opts[-1][0]
        if len(opts)>1 and opts[-2][0]=="-":#Maybe not ktup after all?
            if opts[-2]!="-3":#Not ktup
                ktup=""
    else:
        ktup = ""
    if len(ktup)>0:
        opt = string.join(opts[:-1]," ")
    lines = os.popen('fasta34 -H -b 9999999 -n -d 0 -q %s %s %s %s' %
                     (opt,fastafile,database,ktup)).readlines()
    N = len(readonelinefasta(database)[0])
    hits = []
    scoreindex = -3 #Score column
    while 1:
        i=0
        try:
            while lines[i].find("The best scores are:")==-1:
                i += 1
        except IndexError:
            break
        lines = lines[i+1:]
        j = lines.index("\n")
        for line in lines[:j]:
            if float(line.split()[scoreindex])>=cutoff:
                hits.append(line.split()[0])
    return Uniquenames(hits)

def ssearch34(fastafile,database,cutoff=0,opt=""):
    """
    Use ssearch 3.4 to identify homologs to the sequences in fastafile
    in database using the cutoff cutoff and the options defined in
    opt. The default cutoff is 0 and the default is to use the
    program's default options.

    The output is a list of hits.
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
    lines = os.popen('ssearch34 -H -b 9999999 -n -d 0 -q %s %s %s %s' %
                     (opt,fastafile,database,ktup)).readlines()
    N = len(readonelinefasta(database)[0])
    hits = []
    while 1:
        i=0
        try:
            while lines[i].find("The best scores are:")==-1:
                i += 1
        except IndexError:
            break
        lines = lines[i+1:]
        j = lines.index("\n")
        for line in lines[:j]:
            if float(line.split()[-3])>=cutoff:
                hits.append(line.split()[0])
    return Uniquenames(hits)

def hmmer2(fastafile,database,cutoff=0,opt=""):
    """
    Use HMMer 2.3.2 to identify homologs to the sequences in fastafile
    in database using the cutoff cutoff and the options defined in
    opt. The default cutoff is 0 and the default is to use the
    program's default options.

    The output is a list of hits.
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
    hits = []
    for line in lines[:j]:
        if float(line.split()[1])>=cutoff:
            hits.append(line.split()[0])
    return Uniquenames(hits)

def hmmer1(fastafile,database,cutoff=0,opt=""):
    """
    Use old  HMMer 1.8.4 to identify homologs to the sequences in fastafile
    in database using the cutoff cutoff and the options defined in
    opt. The default cutoff is 0 and the default is to use the
    program's default options.

    The output is a list of hits.
    
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
        os.system("hmmt %s.hmm1 %s " % (name,fastafile))
        lines = os.popen("hmms -t -99999 %s.hmm1 %s" % (name,database)).readlines()
    i=-1
    while i<len(lines):
        i = i+1+lines[i+1:].index(" - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n")
        if lines[i+1]=="\n":
            break
    lines = lines[i+2:]
    hits=[]
    for line in lines:
        try:
            if float(line.split()[0])>=cutoff:
                hits.append(line.split()[-1])
        except ValueError:
            pass
    return Uniquenames(hits)

def sam35(fastafile,database,cutoff=0,opt=""):
    """
    Use SAM 3.5 to identify homologs to the sequences in fastafile
    in database using the cutoff cutoff and the options defined in
    opt. The default cutoff is 0 and the default is to use the
    program's default options.

    The output is a list of hits.
    """
    name = fastafile.split(".")
    if len(name)>1:
        name = string.join(name[:-1],".")
    else:
        name = name[0]
    name2 = name + "_" + database.split("/")[-1]
    os.system("buildmodel %s -a RNA -train %s " % (name,fastafile))
    os.system("hmmscore %s %s -i %s.mod -db %s " % (name2,opt,name,database))
    lines = open("%s.dist" % name2,"r").readlines()
    i=0
    while i<len(lines):
        if lines[i].find("Sequence ID")>-1:
            break
        i += 1
    lines = lines[i+1:]
    hits=[]
    for line in lines:
        # Good values are negative!!
        spl=line.split()
        if len(spl)>3 and float(spl[3])<=cutoff:
            hits.append(spl[0])
        elif len(spl)<=3:
            sys.stderr.write(line)

    return Uniquenames(hits)

def samhmmer2(fastafile,database,cutoff=0,opt=""):
    """
    Use SAM 3.5 to create the HMM model from the sequences in
    fastafile and use HMMer 2.3.2 to search the database using the
    cutoff cutoff and the options defined in opt. The default cutoff
    is 0 and the default is to use the program's default options.

    The output is a list of hits.
    """
    name = fastafile.split(".")
    if len(name)>1:
        name = string.join(name[:-1],".")
    else:
        name = name[0]                
    os.system("buildmodel %s -a RNA -train %s " % (name,fastafile))
    os.system("sam2hmmer %s -i %s.mod; mv %s.hmmer %s.hmmersam35" %
              (name, name, name, name))
    lines = os.popen("hmmsearch -E 99999999 %s.hmmersam35 %s" %
                     (name,database)).readlines()
    i=0
    while lines[i].find("Sequence")==-1 or \
              lines[i].find("Description")==-1 or \
              lines[i].find("Score")==-1 or \
              lines[i].find("E-value")==-1:
        i += 1
    lines = lines[i+2:]
    j = lines.index("\n")
    hits = []
    for line in lines[:j]:
        if float(line.split()[1])>=cutoff:
            hits.append(line.split()[0])
    return Uniquenames(hits)

def rsmatch12(fastafile,database,cutoff=0,opt=""):
    """
    Use RSmatch 1.2 to identify homologs to the sequences in fastafile
    in database using the cutoff cutoff and the options defined in
    opt. The default cutoff is 0 and the default is to use the
    program's default options.

    The output is a list of hits.
    """
    #Save results so that they can be reused
    if opt!=0:
        optname = string.join(opt.split(),"_")
    else:
        optname = ""
    filename = fastafile + "_" + database.split("/")[-1] + "_" + optname + ".RSmatch12"
    if os.path.exists(filename):
        sys.stderr.write("Load %s..." % filename)
        (name,score) = pickle.load(open(filename,"r"))
    else:
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

    hits = map(lambda j: name[j],
               filter(lambda i: score[i]>=cutoff, range(len(score))))
    return Uniquenames(hits)

def paralign(fastafile,database,cutoff=0,opt=""):
    """
    Use paralign to identify homologs to the sequences in fastafile in
    database using the cutoff cutoff and the options defined in
    opt. The default cutoff is 0 and the default is to use the
    program's default options.

    The output is a list of hits.
    """
    if opt == 0:
	opt = ""
    (id,seq) = readonelinefasta(fastafile);
    hits = []
    for k in range(len(id)):
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
            for line in lines[i+1:i+j+1]:
                if float(line.split()[3])>=cutoff:
                    hits.append(line.split()[1])
        except ValueError:
            sys.stderr.write("No hits\n")
    return Uniquenames(hits)

#Useful functions of different kinds

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

def readseqstr(fastafile):
    lines = open(fastafile,"r").readlines()
    id = map(lambda x: x[1:].strip(), map(lambda i: lines[i], range(0,len(lines),3)))
    seq = map(lambda x: x.strip(), map(lambda i: lines[i], range(1,len(lines),3)))
    str = map(lambda x: x.strip(), map(lambda i: lines[i], range(2,len(lines),3)))
    return (id,seq,str)

def writeclustal(filename,id,seq):
    """
    Write the list of identities (id) and sequences (seq) in
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
    
if __name__ == "__main__":

    def usage():
        sys.stderr.write("""
        Error: Usage:
        homology.py searchmethod fastafile database cutoff opt
        
        searchmethod is one of ncbiblast, wublast, fasta34, ssearch34,
        hmmer2, hmmer1, sam35, samhmmer2, paralign, rsmatch12

        Using the query fastafile the database is searched with the
        specified search method. The specified cutoff is used and if
        any options to the search method are defined in opt these are
        used. The output is a number showing the fraction of the
        database entries that were among the hits above the cutoff.
        \n""")
        raise SystemExit
    
    if len(sys.argv) <4:
        usage()
                
    method = sys.argv[1]
    opt = ""
    if len(sys.argv)>=6:
        # specified optimisation, an option to send to the search
        #program
        opt = string.join(sys.argv[5:])
    if len(sys.argv)>=5:
        # cutoff level is specified, if not, use 0
        cutoff = float(sys.argv[4])
    else:
        cutoff = 0.0
    if method=="ncbiblast" or method=="wublast" or \
           method =="fasta34" or method == "ssearch34" or \
           method == "hmmer1" or method == "hmmer2" or \
           method == "sam35" or method=="samhmmer2" or \
           method == "rsmatch12" or method == "paralign":
        fastafile = sys.argv[2]
        database = sys.argv[3]
        exec("hits = " + method + "(fastafile,database,cutoff=cutoff,opt=opt)")
        N = int(os.popen("grep '>' %s|wc -l" % database).
                read().split()[0])
        print "%g" % (float(len(hits))/N)
    else:
        usage()
        
