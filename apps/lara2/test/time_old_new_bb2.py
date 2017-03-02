#!/usr/bin/python

'''
    Test Lara 1 versus Lara 2 on Bralibase 2 data sets.
    
    Author: Joerg Winkler <j.winkler@fu-berlin.de>
'''

import re
import os
import sys
import time
import numpy
import subprocess
from Bio import AlignIO
from Bio.Statistics.lowess import lowess
from shutil import copyfile
from matplotlib import use
use('PDF')
import matplotlib.pyplot as plt

def errorhandle(returncode, program='Program'):
  if returncode < 0:
    print >>sys.stderr, program + " was terminated by signal", -returncode
    return True
  else:
    return False

######################
##  SET FILE PATHS  ##
######################

# define macros
SCORES = (SPS, SCI, MPI) = (0, 1, 2)
SCORELBL = ('Sum of Pairs Score (compalignp)', 'Structure Conservation Index (RNAz)', 'Mean Pairwise Identity (RNAz)')
PROGRAMS = (LA1, L2O, L2N, STC, MAF, REF) = (0, 1, 2, 3, 4, 5)
PROGLBL = ('Lara1+TC', 'Lara2+TC', 'Lara2+S::TC', 'SeqAn::TC', 'MAFFT   ', 'Reference')

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

# compiled MAFFT program
mafft_bin = os.path.join(work_dir, "..", "..", "..", "..", "..", "apps", "tcoffee", "plugins", "linux", "mafft")
# compiled RNAz program ( http://www.tbi.univie.ac.at/~wash/RNAz/ )
rnaz_bin = os.path.join(work_dir, "..", "..", "..", "..", "..", "apps", "rnaz", "bin", "RNAz") # adapt
# compiled compalignp program ( http://www.biophys.uni-duesseldorf.de/bralibase/compalignp.tgz )
compali_bin = os.path.join(work_dir, "..", "..", "..", "..", "..", "apps", "compalignp", "compalignp") # adapt

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

progtime = [0.0] * len(PROGRAMS)
stats = {}
had_err = []

for (infile, outfile) in files:
  errors = [False] * 7
  basename = os.path.basename(outfile)
  print >>sys.stderr, "  processing", basename
  
  # run Lara1
  t = time.time()
  proc = subprocess.Popen([oldlara_bin, "-i", infile, "-w", outfile + str(LA1) + ".fasta"],\
         bufsize=-1, executable=oldlara_bin, stdout=subprocess.PIPE, shell=False, cwd=oldlara_dir)
  proc.communicate()
  progtime[LA1] += time.time() - t
  errors[0] = errorhandle(proc.returncode, oldlara_bin + " " + infile)
      
  # run Lara2
  t = time.time()
  proc = subprocess.Popen([newlara_bin, "-i", infile, "-td", results_dir, "-t", "1"],\
         bufsize=-1, executable=newlara_bin, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=False)
  proc.communicate()
  lara2time = time.time() - t
  errors[1] = errorhandle(proc.returncode, oldlara_bin + " " + infile)
  
  if errors[1]:
    had_err.append(errors)
    continue
  
  # old tcoffee
  t = time.time()
  proc = subprocess.Popen([oldtcof_bin, "-in", tc_tempfile, "-case=upper", "-output fasta", "-clean_seq_name 1",\
         "-outfile", outfile + str(L2O) + ".fasta", "-newtree", outfile + str(L2O) + ".dnd"],\
         bufsize=-1, executable=oldtcof_bin, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=False)
  proc.communicate()
  progtime[L2O] += lara2time + time.time() - t
  errors[2] = errorhandle(proc.returncode, oldtcof_bin + " " + infile)
  
  # new tcoffee
  t = time.time()
  proc = subprocess.Popen([newtcof_bin, "-s", infile, "-l", tc_tempfile, "-m", "global", "-a", "iupac", "-b", "wavg",\
         "-o", outfile + str(L2N) + ".fasta"], bufsize=-1, executable=newtcof_bin, stdout=subprocess.PIPE, shell=False)
  proc.communicate()
  progtime[L2N] += lara2time + time.time() - t
  errors[3] = errorhandle(proc.returncode, newtcof_bin + " " + tc_tempfile)
  
  # SeqAn::TCoffee without Lara
  t = time.time()
  proc = subprocess.Popen([newtcof_bin, "-s", infile, "-a", "iupac", "-b", "wavg", "-o", outfile + str(STC) + ".fasta"],\
         bufsize=-1, executable=newtcof_bin, stdout=subprocess.PIPE, shell=False)
  proc.communicate()
  progtime[STC] += time.time() - t
  errors[4] = errorhandle(proc.returncode, newtcof_bin + " " + infile)

  # MAFFT
  t = time.time()
  f = open(outfile + str(MAF) + ".fasta", "w")
  proc = subprocess.Popen([mafft_bin, infile],\
         bufsize=-1, executable=mafft_bin, stderr=subprocess.PIPE, stdout=f, shell=False)
  proc.communicate()
  progtime[MAF] += time.time() - t
  f.close()
  errors[5] = errorhandle(proc.returncode, mafft_bin + " " + infile)
  
  if any(errors):
    had_err.append(errors)
    continue
  
  # Reference
  refalignment = infile.replace("unaligned", "structural-no_str")
  copyfile(refalignment, outfile + str(REF) + ".fasta")
   
  # transform into Clustal format and run RNAz analysis
  for file in [outfile + str(x) for x in PROGRAMS]:
    ali = AlignIO.read(file + ".fasta", "fasta")
    AlignIO.write([ali], file + ".aln", "clustal")

    proc = subprocess.Popen([rnaz_bin, "-o", file + ".stat", "-n", file + ".aln"],\
           bufsize=-1, executable=rnaz_bin, shell=False)
    proc.communicate()
    errors[6] = errorhandle(proc.returncode, rnaz_bin + " " + file + ".aln") or errors[4]

  if any(errors):
    had_err.append(errors)

for (pname, tm) in zip(PROGLBL, progtime)[:-1]:
  print('Total time for ' + pname + '\t {} seconds.'.format(tm))

if len(had_err) > 0:
  print("There were errors: " + str(had_err))
  exit(1)


##########################
##  ALIGNMENT ANALYSIS  ##
##########################

print ("\nAnalyze alignments...")
for (infile, outfile) in files:
  errors = [False] * 2
  basename = os.path.basename(outfile)

  stats[basename] = ([],[],[],[],[],[])
  
  # run compalignp
  refalignment = infile.replace("unaligned", "structural-no_str")
  for (i, file) in [(x, outfile + str(x) + ".aln") for x in PROGRAMS]:
    proc = subprocess.Popen([compali_bin, "-t", file, "-r", refalignment],\
           bufsize=-1, executable=compali_bin, stdout=subprocess.PIPE, shell=False)
    stats[basename][i].append(float(proc.communicate()[0]))
    errors[0] = errorhandle(proc.returncode, compali_bin + " " + file) or errors[0]
  
  if any(errors):
    had_err.append(errors)
    continue
  
  # extract statistics from files
  for (i, file) in [(x, outfile + str(x) + ".stat") for x in PROGRAMS]:
    (sci, mpi) = (None, None)
    for line in open(file, 'r'):
      if line.startswith(" Mean pairwise identity:"):
        mpi = float(re.findall(r"(-?\d+\.\d*)", line)[0])
      if line.startswith(" Structure conservation index:"):
        sci = float(re.findall(r"(-?\d+\.\d*)", line)[0])
    
    if sci == None or mpi == None:
      print >>sys.stderr, "Could not parse file " + file
      errors[1] = True
      continue
      
    stats[basename][i].extend([sci, mpi])
  
  if any(errors):
    had_err.append(errors)
if len(had_err) > 0:
  print("There were errors: " + str(had_err))
  exit(1)

###############
##  RESULTS  ##
###############

# provide data sorted by MPI
data = [zip(*x) for x in zip(*stats.values())]
for i in PROGRAMS:
  (data[i][MPI], data[i][SPS], data[i][SCI]) = (list(x) for x in zip(*sorted(zip(data[i][MPI], data[i][SPS], data[i][SCI]))))
  
# print cumulative values
def _make_stat_str(values):
  return "\t Mean "   + str(tuple([round(numpy.mean(x),2) for x in values]))\
       + "\t Median " + str(tuple([round(numpy.median(x),2) for x in values]))\
       + "\t StdDev " + str(tuple([round(numpy.std(x),2) for x in values]))
       
print "Values for (Sum of Pairs Score, Structure Conservation Index, Mean Pairwise Identity):"
for (pname,pval) in zip(PROGLBL,data):
  print "  " + pname + _make_stat_str(pval)

# set view area for plots
view = ([28, 100, 0.3, 1], [28, 100, 0, 1.4])

for score in (SPS, SCI):
  # calculate lowess function
  x = [numpy.array(prog[MPI], numpy.float) for prog in data]
  y = [prog[score] for prog in data]
  f = [lowess(x[i], numpy.array(y[i], numpy.float)) for i in PROGRAMS]
  
  # plot data
  for (prog, color) in ((LA1,'b'), (L2O,'g'), (L2N,'r'), (STC,'c'), (MAF,'m'), (REF,'k')):
    plt.plot(map(int, x[prog]), y[prog], "." + color, ms=3)
    plt.plot(x[prog], f[prog], color, label=PROGLBL[prog])
    
  plt.axis(view[score])
  plt.xlabel(SCORELBL[MPI])
  plt.ylabel(SCORELBL[score])
  plt.legend(loc="lower right")
  plt.savefig(os.path.join(results_dir, "figure" + str(score) + ".pdf"))
  print "Plot saved to results/figure" + str(score) + ".pdf"
  plt.clf()

