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

newlara_bin = os.path.join(work_dir, "..", "..", "..", "..", "build", "bin", "laragu")
newlara_res = os.path.join(work_dir, "benchmarks", "results-lara2")

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
for (infile, out1, out2) in files:
  print >>sys.stderr, "  processing", os.path.basename(out1)
  
  # run Lara1
  t = time.time()
  proc = subprocess.Popen([oldlara_bin, "-i", infile, "-w", out1, "-c"],\
         bufsize=-1, executable=oldlara_bin, stdout=subprocess.PIPE, shell=False, cwd=oldlara_dir)
  proc.communicate()
  lara1time += time.time() - t
  if proc.returncode < 0:
      print >>sys.stderr, "Lara was terminated by signal", -proc.returncode
      exit(proc.returncode)
      
  ''' 
  JW working on 2 problems: 
  a) Lara2 does not take absolute input file paths 
  b) Bralibase files contain DNA letter 'T' (parse error with Lara2)
  '''
      
  # run Lara2
  t = time.time()
  #proc = subprocess.Popen([newlara_bin, "-i", infile],\
  #       bufsize=-1, executable=newlara_bin, stdout=subprocess.PIPE, shell=False)
  #print proc.communicate()[0]
  lara2time += time.time() - t
  #if proc.returncode < 0:
  #    print >>sys.stderr, "Lara was terminated by signal", -proc.returncode
  #    exit(proc.returncode)

print('\nTotal time for Lara1: {} seconds.'.format(lara1time))
print('\nTotal time for Lara2: {} seconds.'.format(lara2time))
