#!/usr/bin/env python

import sys, os, string, operator

def wublast(queryfile,database,cutoff=0,opt=""):
    """
    Run WU-blast on all the sequences in fastafile against database
    using the cutoff cutoff and the options defined in opt. The
    default cutoff is 0 and the default is to use the program's
    default options.

    Prints a gff-like output to stdout
    """
    lines = os.popen('blastn %s %s -noseqs -gspmax 0 -hspmax 0 %s' %
                     (database,fastafile,opt) +
                     "|egrep '(Sbjct|Query=|Strand =|Score =)'").\
                     readlines()
    seqname = queryfile.split("/")[-1]
    try:
        source = queryfile.split("/")[-2]
    except IndexError:
        source = "id"
    for line in lines:
        if line.find("Sbjct")>-1:
            if score>=cutoff:
                positions = line.split()
                positions=[int(positions[1]),int(positions[3])]
                positions.sort()
                start_end = "%d\t%d" % (positions[0],positions[1])
                #gff format: <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments]
                print "%s\t%s\t%s\t%s\t%0.4g\t%s" %\
                      (seqname,source,query,start_end,score,strand)
        elif line.find("Score = ")>-1:
            score = float(line.split()[2])
        else:
            i = line.find("Strand = ")
            if i>-1:
                sp = line[i+9:].split()
                strand = (sp[0] == sp[2])
                if strand:
                    strand = "+"
                else:
                    strand = "-"
            elif line.find("Query=")>-1:
                query = string.join(line.split()[1:])

def fasta34(queryfile,database,cutoff=0,opt=""):
    """
    Use fasta34 to identify regions in the database that are
    homologous to the sequences in the queryfile.
    Use the score cutoff.
    opt is options to send to fasta, such as -3 to search only the
    top strand etc.

    Prints a gff-like output to stdout
    """
    # Sort out the options
    opts = opt.split()
    if len(opts)>0 and opts[-1][0]!="-":
        # ktup?
        ktup = opts[-1][0]
        if len(opts)>1 and opts[-2][0]=="-":
            #Maybe not ktup after all?
            if opts[-2]!="-3" and opts[-1]!="-U" and opts[-2]!="-z11" :
                #Not ktup
                ktup=""
    else:
        ktup = ""
    if len(ktup)>0:
        opt = string.join(opts[:-1]," ")
    seqname = queryfile.split("/")[-1]
    try:
        source = queryfile.split("/")[-2]
    except IndexError:
        source="id"
    # Do the fasta search with the given options, make sure to report
    # all (or at least most of) the hits.
    lines = os.popen("fasta34 -H -b 9999999 -n -q %s %s %s %s" %
                     (opt,queryfile,database,ktup) +
                     "| egrep '(>>>|overlap|initn)'").\
                     readlines()
    # Unless the option -3 is used the above will search both strands,
    # reporting them as [f] and [r]
    for line in lines:
        if line.find("opt")>-1:
            score = int(line.split("opt:")[1].split()[0])
        elif line.find("overlap")>-1:
            if score<cutoff:
                pass
            else:
                sp = line.split( )[-1].replace("(","").replace(")","").split(":")
                pos = map(lambda x: map(lambda y: int(y),x.split("-")),sp)
                if operator.xor(pos[0][0]-pos[0][1]<0,pos[1][0]-pos[1][1]<0):
                    sign = "-"
                else:
                    sign = "+"
                positions = "%s\t%s" % (min(pos[1]),max(pos[1]))
                #gff format: <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments]
                print "%s\t%s\t%s\t%s\t%0.4g\t%s" %\
                      (seqname,source,query,positions,score,sign)
        elif line.find(">>>")>-1:
            query = line.split()[0].split(">>>")[-1]

def ssearch34(queryfile,database,cutoff=0,opt=""):
    """
    Use ssearch34 to identify regions in the database that are
    homologous to the sequences in the queryfile.
    Use the score cutoff.
    opt is options to sent to ssearch, such as -3 to search only the
    top strand etc.
    
    Prints a gff-like output to stdout
    """
    # Sort out the options
    opts = opt.split()
    if len(opts)>0 and opts[-1][0]!="-":
        # ktup?
        ktup = opts[-1][0]
        if len(opts)>1 and opts[-2][0]=="-":
            #Maybe not ktup after all?
            if opts[-2]!="-3" and opts[-1]!="-U" and opts[-2]!="-z11" :
                #Not ktup
                ktup=""
    else:
        ktup = ""
    if len(ktup)>0:
        opt = string.join(opts[:-1]," ")
    seqname = queryfile.split("/")[-1]
    try:
        source = queryfile.split("/")[-2]
    except IndexError:
        source = "id"
    # Do the ssearch search with the given options, make sure to report
    # all (or at least most of) the hits.
    lines = os.popen("ssearch34 -H -b 9999999 -n -q %s %s %s %s" %
                     (opt,queryfile,database,ktup) +
                     "| egrep '(>>>|overlap|s-w opt)'").\
                     readlines()
    # Unless the option -3 is used the above will search both strands,
    # reporting them as [f] and [r]
    for line in lines:
        if line.find("s-w")>-1:
            if line.split()[0]=="rev-comp":
                sign = "-"
            else:
                sign = "+"
            score = int(line.split("opt:")[1].split()[0])
        elif line.find("overlap")>-1:
            if score<cutoff:
                pass
            else:
                sp = line.split( )[-1].replace("(","").replace(")","").split(":")
                pos = map(lambda x: map(lambda y: int(y),x.split("-")),sp)
                positions = "%s\t%s" % (min(pos[1]),max(pos[1]))
                #gff format: <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments]
                print "%s\t%s\t%s\t%s\t%0.4g\t%s" %\
                      (seqname,source,query,positions,score,sign)
        elif line.find(">>>")>-1:
            query = line.split()[0].split(">>>")[-1]

def hmmer2(queryfile,database,cutoff=0,opt=""):
    """
    Use HMMer 2.3.2 to identify homologs to the sequences in fastafile
    in database using the cutoff cutoff and the options defined in
    opt. The default cutoff is 0 and the default is to use the
    program's default options.

    Prints a gff-like output to stdout
    """
    name = fastafile.split(".")
    if len(name)>1:
        name = string.join(name[:-1],".")
    else:
        name = name[0]
    seqname = queryfile.split("/")[-1]
    try:
        source = queryfile.split("/")[-2]
    except IndexError:
        source = "id"
    #Reverse database?
    if database.find("rev")>-1:
        sign="-"
    else:
        sign="+"
    if not os.path.exists(name+"_pro.aln"):
        os.system("java -jar /users/lcb/evaf/bin/ProAlign/" +
                  "ProAlign_0.5a1.jar -nogui  " +
                  "-seqfile=%s -newtree -outfile=%s.pir -outformat=pir" %
                  (fastafile,name))
        (id,seq) = readfasta(name+".pir")
        writeclustal(name+"_pro.aln",map(lambda x: x.split("DL;")[-1],id),
                     map(lambda x: x.replace("*",""),seq))
    os.system("hmmbuild %s -F -n %s %s.hmm2 %s_pro.aln > /dev/null" %
              (opt,name,name,name))
    #Searches only top strand
    lines = os.popen("hmmsearch -A0 -E 99999999 %s.hmm2 %s" % (name,database)).readlines()
    i=0
    # A line
    #Sequence Domain  seq-f seq-t    hmm-f hmm-t      score  E-value
    #is seen before the result list
    while lines[i].find("Sequence")==-1 or \
              lines[i].find("Domain")==-1 or \
              lines[i].find("seq-f")==-1 or \
              lines[i].find("hmm-f")==-1:
        i += 1
    lines = lines[i+2:]
    j = lines.index("\n")
    hits = []
    for line in lines[:j]:
        spl = line.split()
        score = float(spl[8])
        if score<cutoff:
            break
        positions = "%s\t%s" % (spl[2],spl[3])
        #gff format: <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments]
	print "%s\t%s\t%s\t%s\t%0.4g\t%s" % (seqname,source,"-",positions,score,sign)

def sam35(queryfile,database,cutoff=0,opt=""):
    """
    Use SAM 3.5 to identify homologs to the sequences in fastafile
    in database using the cutoff cutoff and the options defined in
    opt. The default cutoff is 0 and the default is to use the
    program's default options.

    Output is sent to stdout.
    """
    name = queryfile.split(".")
    if len(name)>1:
        name = string.join(name[:-1],".")
    else:
        name = name[0]

    seqname = queryfile.split("/")[-1]
    try:
        source = queryfile.split("/")[-2]
    except IndexError:
        source = "id"
    #Reverse database?
    if database.find("rev")>-1:
        sign="-"
    else:
        sign="+"
    
    name2 = name + "_" + database.split("/")[-1]
    os.system("buildmodel %s -a RNA -train %s " %
              (name,queryfile))
    os.system("hmmscore %s -i %s.mod -db %s %s -select_mdalign 8 -select_md 8 " %
              (name2,name,database,opt))
    lines = open("%s.mstat" % name2,"r").readlines()
    i=0
    while i<len(lines):
        if lines[i].find("Sequence ID")>-1:
            break
        i += 1
    lines = lines[i+1:]
    hits=[]
    for line in lines:
        # Good values are negative!!
        score = float(line.split()[3])
        if score>cutoff:
            break
        spl=line.split()[0].split("_")
        start=10000*int(spl[0].split(":")[-1].split("-")[0])
        if sign=="+":
            positions = "%d\t%d" % tuple(map(lambda x: start + int(x),
                                             spl[-1].split(":")))
        else:#sign=="-"
            length=10000
            if start==129900000:
                length+=1
                sys.stderr.write(length)
            positions = map(lambda x: length - int(x) + 1,spl[-1].split(":"))
            positions.reverse()
            positions = "%d\t%d" % tuple(positions)
        print "%s\t%s\t%s\t%s\t%0.4g\t%s" % (seqname,source,"-",positions,score\
,sign)


#Some useful functions;

def readfasta(fastafile):
    lines = open(fastafile,"r").readlines()
    id = []
    seq = []
    for line in lines:
        if len(line)==0:
            pass
        if line[0]==">":
            id.append(line[1:].strip())
            seq.append("")
        else:
            seq[-1] += line.strip()
    if len(id)!=len(seq):
        print "Error: There is something wrong with the file %s." % fastafile
        raise SystemExit
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


#______________________________________________________________________#
# Main program

if __name__ == "__main__":

    def usage():
        sys.stderr.write("""
        Error: Usage:
        getgff.py searchmethod queryfile database cutoff opt
        \n""")
        raise SystemExit
    
    if len(sys.argv) <4:
        usage()
                
    method = sys.argv[1]
    opt = ""
    if len(sys.argv)>=6:
        # specified optimisation, an option to send to
        #the search program
        opt = string.join(sys.argv[5:])
        sys.stderr.write("opt = %s\n" % opt)
    if len(sys.argv)>=5:
        # cutoff level is specified, if not, use 0
        cutoff = float(sys.argv[4])
    else:
        cutoff = 0.0
    if method=="wublast" or method =="fasta34" or \
           method == "ssearch34" or method == "hmmer2" or \
           method == "sam35":
        fastafile = sys.argv[2]
        database = sys.argv[3]
        exec(method + "(fastafile,database,cutoff=cutoff,opt=opt)")
    else:
        usage()
        

