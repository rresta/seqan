#!/usr/bin/python

'''
    Test Lara 1 versus Lara 2 on Bralibase 2 data sets.
    
    Author: Joerg Winkler <j.winkler@fu-berlin.de>
'''

import os
import sys
import time
import subprocess
from Bio import AlignIO

def errorhandle(returncode, program='Program'):
  if returncode < 0:
    print >>sys.stderr, program + " was terminated by signal", -returncode
    return True
  else:
    return False

######################
##  SET FILE PATHS  ##
######################

work_dir = os.getcwd()
results_dir  = os.path.join(work_dir, "results")
# old lara and old tcoffee
oldlara_dir = os.path.join(work_dir, "lara-1.3.2")
oldlara_bin = os.path.join(oldlara_dir, "lara")
oldtcof_bin = os.path.join(oldlara_dir, "t_coffee", "t_coffee_5.05")


# new lara and new tcoffee
seqan_dir = os.path.join(work_dir, "..", "..", "..", "..", "build", "bin") # adapt this to your system!
newlara_bin = os.path.join(seqan_dir, "laragu")
newtcof_bin = os.path.join(seqan_dir, "seqan_tcoffee")
tc_tempfile = os.path.join(results_dir, "tcoffeLara.lib")

# compiled RNAz program ( http://www.tbi.univie.ac.at/~wash/RNAz/ )
rnaz_bin = os.path.join(work_dir, "..", "..", "..", "..", "..", "apps", "rnaz", "bin", "RNAz") # adapt this

########################
##  CREATE FILE LIST  ##
########################

files = []
for d in "g2intron", "rRNA", "tRNA", "U5":
  p = os.path.join(work_dir, "benchmarks", "bralibase2", "data-set1", d, "unaligned")
  files.extend([(os.path.join(p, f), os.path.join(results_dir, d + "-" + f)) for f in os.listdir(p)])

if not os.path.isdir(results_dir):
  os.mkdir(results_dir)
  
print(str(len(files)) + " alignments to compute.")

############################
##  CALCULATE ALIGNMENTS  ##
############################

lara1time = 0.0
lara2time = 0.0
oldtctime = 0.0
newtctime = 0.0
had_err = []
for (infile, outfile) in files:
  errors = [False] * 5
  print >>sys.stderr, "  processing", os.path.basename(outfile)
  
  # run Lara1
  t = time.time()
  proc = subprocess.Popen([oldlara_bin, "-i", infile, "-w", outfile + "1.fasta", "-c"],\
         bufsize=-1, executable=oldlara_bin, stdout=subprocess.PIPE, shell=False, cwd=oldlara_dir)
  proc.communicate()
  lara1time += time.time() - t
  errors[0] = errorhandle(proc.returncode, oldlara_bin + " " + infile)
      
  # run Lara2
  t = time.time()
  proc = subprocess.Popen([newlara_bin, "-i", infile, "-td", results_dir],\
         bufsize=-1, executable=newlara_bin, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=False)
  proc.communicate()
  lara2time += time.time() - t
  errors[1] = errorhandle(proc.returncode, oldlara_bin + " " + infile)
  
  if errors[1]:
    had_err.append(errors)
    continue
  
  # old tcoffee
  t = time.time()
  proc = subprocess.Popen([oldtcof_bin, "-in", tc_tempfile, "-case=upper", "-output fasta", "-clean_seq_name 1",\
         "-outfile", outfile + "2.fasta", "-newtree", outfile + "2.dnd"],\
         bufsize=-1, executable=oldtcof_bin, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=False)
  proc.communicate()
  oldtctime += time.time() - t
  errors[2] = errorhandle(proc.returncode, oldtcof_bin + " " + infile)
  
  # new tcoffee
  t = time.time()
  proc = subprocess.Popen([newtcof_bin, "-s", infile, "-l", tc_tempfile, "-m", "global", "-a", "iupac", "-o",\
         outfile + "3.fasta", "-b", "wavg"], bufsize=-1, executable=newtcof_bin, stdout=subprocess.PIPE, shell=False)
  proc.communicate()
  newtctime += time.time() - t
  errors[3] = errorhandle(proc.returncode, newtcof_bin + " " + infile)
  
  if any(errors):
    had_err.append(errors)
    continue

  ##########################
  ##  ALIGNMENT ANALYSIS  ##
  ##########################
  
  ali = AlignIO.read(outfile + "1.fasta", "fasta")
  AlignIO.write([ali[:-1]], outfile + "1.aln", "clustal")
  ali = AlignIO.read(outfile + "2.fasta", "fasta")
  AlignIO.write([ali], outfile + "2.aln", "clustal")
  ali = AlignIO.read(outfile + "3.fasta", "fasta")
  AlignIO.write([ali], outfile + "3.aln", "clustal")
  
  for r in "1","2","3":
    proc = subprocess.Popen([rnaz_bin, "-o", outfile + r + ".stat", "-n", outfile + r + ".aln"],\
           bufsize=-1, executable=rnaz_bin, shell=False)
    proc.communicate()
    errors[4] = errorhandle(proc.returncode, rnaz_bin + " " + outfile + r + ".aln") or errors[4]
    
  if any(errors):
    had_err.append(errors)
    continue

if len(had_err) > 0:
  print("There were errors: " + str(had_err))
  
print('Total time for Lara1:          {} seconds.'.format(lara1time))
print('Total time for Lara2:          {} seconds.'.format(lara2time))
print('Total time for TCoffee:        {} seconds.'.format(oldtctime))
print('Total time for SeqAn::TCoffee: {} seconds.'.format(newtctime))


