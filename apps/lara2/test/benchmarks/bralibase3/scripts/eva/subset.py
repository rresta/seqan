#!/usr/bin/env python

#subset.py produces N subsets (subalignments) of n sequences from an
#input alignment (in fasta format), N and n are user defined. The
#sequences in a subset all have pairwise sequence identities to each
#other within a user defined range. Additionaly a maximum number of
#overlapping sequences between two different subsets can be defined.
#
#Copyright (C) 2006  Eva Freyhult
#
#This program is free software; you can redistribute it and/or
#modify it under the terms of the GNU General Public License
#as published by the Free Software Foundation; either version 2
#of the License, or (at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program; if not, write to the Free Software
#Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#
# Contact information: Eva.Freyhult@lcb.uu.se, http://www.lcb.uu.se/~evaf

import random, string, sys, os

def extract(fastafile,ident,N,n,overlap=-1):
    """
    extract(fastafile,ident,N,n,overlap=-1)

    where ident is a 2x1 array with a lower and an upper bound on the
    pairwise sequence identity, N is the number of subsets wanted, and
    n is the number of sequences to be included in each subset
    overlap is the number of elements that can overlap between two sets,
    if overlap is a negative number (the default) two subsets can overlap
    with any number of sequences (but two subsets will never be completely
    identical)

    """

    if overlap<0:
        overlap = n-1
        
    (I,S) = readfasta(fastafile)
    P = pairwiseids(fastafile)
    M = len(P)
    
    count = 0
    sets = []
    possible = filter(lambda x: len(filter(lambda y:ident[0]<=y and
                                           y<=ident[1],P[x]))>=n,
                      range(M))
    if len(possible)<1:
        return sets
    while len(sets) < N and count<50*N:
        seq = possible
        s = []
        s.append(random.choice(seq))
        seq = filter(lambda x:ident[0]<=P[s[-1]][x] and
                     P[s[-1]][x]<=ident[1],seq)
        while len(s)<n and len(seq)>0:
            s.append(random.choice(seq))
            seq = filter(lambda x:ident[0]<=P[s[-1]][x] and
                         P[s[-1]][x]<=ident[1],seq)
        s.sort()
        if len(s)==n and reduce(lambda y,x:y and intersect(s,x)<=overlap,
                                sets,1):
            #Save set only if it contains n elements and is not
            #already among the sets
            sets.append(s)
        count += 1
    return sets

def intersect(S,T):
    #Number of elements not in common to S and T
    #The elements in the two sets must be sorted in increasing order
    m=len(S);n=len(T)
    d=0
    i=0;j=0
    while i<m and j<n:
        if S[i]<T[j]:
            i+=1
        elif S[i]>T[j]:
            j+=1
        else:
            i+=1
            j+=1
            d+=1
    return d
    
def print_sets(fastafile,sets):
    (I,S) = readfasta(fastafile)
    filename = string.join(fastafile.split('.')[:-1],'.')
    if len(filename)==0:
        filename = fastafile
    for i in range(len(sets)):
        writeonelinefasta(filename + "%.2d.fa" % i,
                    map(lambda x:I[x],sets[i]),
                    map(lambda x:S[x].replace('-',''),sets[i]))
        
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

def writeonelinefasta(fastafile,id,seq):
    file=open(fastafile,"w")
    map(lambda x: file.write(">%s\n%s\n" % x), zip(id,seq))
    file.close()

def pairwiseid(S1,S2):
    if len(S1)!=len(S2):
        print "Seqences have different lengths"
        print len(S1), len(S2)
        raise IndexError
    #Remove gap only positions
    S = filter(lambda x:x[0]!="-" or x[1]!="-",zip(S1,S2))
    id = reduce(lambda s,x:s + (x[0]==x[1]),S,0.0)/len(S)
    return id

def pairwiseids(fastafile):
    if os.path.exists(fastafile + ".pwi"):
        print "Loading file %s.pwi" % fastafile
        file = open(fastafile + ".pwi","r")
        Pi = map(lambda x:map(lambda y:float(y),x.split()),file.readlines())
        file.close()
        return Pi
    (I,S) = readfasta(fastafile)
    n = len(S)
    Pi = [n*[1.0] for i in range(n)]
    for i in range(n):
        for j in range(i,n):
            Pi[i][j] = pairwiseid(S[i],S[j])
            Pi[j][i] = Pi[i][j]
    file = open(fastafile + ".pwi","w")
    map(lambda x:map(lambda y:
                     file.write("%0.2g\t" % y),x) and file.write("\n"),Pi)
    file.close()
    return Pi
    
if __name__ == "__main__":

    if len(sys.argv)<6:
        print """
        Usage:
        subset.py fastafile lid uid N n o
        lid and uid = lower and upper bound on %pairwise identity
        N = number of sets wanted (might not be possible to get)
        n = number of sequences in each set
        o = number of elements allowed to overlap between two sets
            (default is n-1)
        """
    else:
        fastafile = sys.argv[1]
        ident = float(sys.argv[2]),float(sys.argv[3])
        if ident[1]>1:
            ident = map(lambda x:x/100,ident)
        N = int(sys.argv[4])
        n = int(sys.argv[5])
        if len(sys.argv)>6:
            o = sys.argv[6]
            sets = extract(fastafile,ident,N,n,o)
        else:
            sets = extract(fastafile,ident,N,n)
        print_sets(fastafile,sets)
