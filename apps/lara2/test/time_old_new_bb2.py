#!/usr/bin/python

'''
    Test Lara 1 versus Lara 2 on Bralibase 2 data sets.
    
    Author: Joerg Winkler <j.winkler@fu-berlin.de>
'''

import os
import sys
import time
import subprocess

# set file paths
work_dir = os.getcwd()
oldlara_dir = os.path.join(work_dir, "lara-1.3.2")
oldlara_bin = os.path.join(oldlara_dir, "lara")
oldlara_res = os.path.join(work_dir, "benchmarks", "results-lara1")
oldtcof_bin = os.path.join(oldlara_dir, "t_coffee", "t_coffee_5.05")

newlara_bin = os.path.join(work_dir, "..", "..", "..", "..", "build", "bin", "laragu")
newlara_res = os.path.join(work_dir, "benchmarks", "results-lara2")
newtcof_bin = os.path.join(work_dir, "..", "..", "..", "..", "build", "bin", "seqan_tcoffee")
tc_tempfile = os.path.join(newlara_res, "tcoffeLara.lib")

# create list of files
files = []
for d in "g2intron", "rRNA", "tRNA", "U5":
  p = os.path.join(work_dir, "benchmarks", "bralibase2", "data-set1", d, "unaligned")
  files.extend([(os.path.join(p, f),\
                 os.path.join(oldlara_res, d + "-" + f),\
                 os.path.join(newlara_res, d + "-" + f)) for f in os.listdir(p)])

if not os.path.isdir(oldlara_res):
  os.mkdir(oldlara_res)
if not os.path.isdir(newlara_res):
  os.mkdir(newlara_res)
  
print(str(len(files)) + " alignments to compute.")

# process all files
lara1time = 0.0
lara2time = 0.0
oldtctime = 0.0
newtctime = 0.0
errors = 0
for (infile, out1, out2) in files:
  print >>sys.stderr, "  processing", os.path.basename(out1)
  
  # run Lara1
  t = time.time()
  proc = subprocess.Popen([oldlara_bin, "-i", infile, "-w", out1, "-c"],\
         bufsize=-1, executable=oldlara_bin, stdout=subprocess.PIPE, shell=False, cwd=oldlara_dir)
  proc.communicate()
  lara1time += time.time() - t
  if proc.returncode < 0:
      print >>sys.stderr, "Lara1 was terminated by signal", -proc.returncode
      errors += 1
      
  # run Lara2
  t = time.time()
  proc = subprocess.Popen([newlara_bin, "-i", infile, "-td", newlara_res],\
         bufsize=-1, executable=newlara_bin, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=False)
  proc.communicate()
  lara2time += time.time() - t
  if proc.returncode < 0:
      print >>sys.stderr, "Lara2 was terminated by signal", -proc.returncode
      errors += 1
      continue
  
  # old tcoffee
  t = time.time()
  proc = subprocess.Popen([oldtcof_bin, "-in", tc_tempfile, "-case=upper", "-output fasta", "-clean_seq_name 1",\
         "-outfile", out2 + ".aln", "-newtree", out2 + ".dnd"],\
         bufsize=-1, executable=oldtcof_bin, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=False)
  proc.communicate()
  oldtctime += time.time() - t
  if proc.returncode < 0:
      print >>sys.stderr, "T-Coffee was terminated by signal", -proc.returncode
      errors += 1
  
  # new tcoffee
  t = time.time()
  proc = subprocess.Popen([newtcof_bin, "-s", infile, "-l", tc_tempfile, "-m", "global", "-a", "iupac",\
         "-o", out2 + ".msf", "-b", "wavg"], bufsize=-1, executable=newtcof_bin, stdout=subprocess.PIPE, shell=False)
  proc.communicate()
  newtctime += time.time() - t
  if proc.returncode < 0:
      print >>sys.stderr, "SeqAn::T-Coffee was terminated by signal", -proc.returncode
      errors += 1

if errors > 0:
  print('There were {} errors.'.format(errors))
print('Total time for Lara1:          {} seconds.'.format(lara1time))
print('Total time for Lara2:          {} seconds.'.format(lara2time))
print('Total time for TCoffee:        {} seconds.'.format(oldtctime))
print('Total time for SeqAn::TCoffee: {} seconds.'.format(newtctime))
