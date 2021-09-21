#!/usr/bin/env python3
# Written by August Woerner
#
# In short: Shapeit4 (and other phasing algorithms) will drop some variants
# In particular, it will drop tri-allelic variants (annoyingly including the reference allele in that assessment).
# 
# This script is used to reconcile and recover variants lost by Shapeit4.
# If the SNPs have been phased already (eg, by WhatsHap), this program
# will further retain the (correct) phase in the final output
#
# This is *not* a general use program however. It is intended to be run in the exome-genome -> proteome pipeline
# (as it relies on bcftools to do a lot of the heavy lifting)



import sys
import os
import gzip
from collections import OrderedDict


if len(sys.argv) < 2:
  print("I need a filename!", file=sys.stderr)
  exit(1)

filename = sys.argv[1]
out = gzip.open(filename=filename, mode="wt", compresslevel=9)

finalDict = OrderedDict()


phasesetOrientation={}

lastIDX=-1
individLUT = {}
for line in sys.stdin:
  if line.startswith("##"):
    out.write(line)
    continue
  elif line.startswith("#"): # fix the vcf header to just contain each individual once

    header = line.rstrip().split("\t")
    LUT = [-1]*len(header) # when filled out, if !=-1, then it refers to the index of the other genotype needed
    
    for i in range(9, len(header)):
      who = header[i]
      if who[0] >= '0' and who[0] <= '9':
        who=who.split(":")[1]

        if who not in individLUT:
          print("Should never happen. Problem with individual:", who, file=sys.stderr)
          exit(1)
        idx = individLUT[who]

        if LUT[idx] != -1:
          print("Also hould never happen. Too many entries for ", who, file=sys.stderr)
          exit(1)
        
        LUT[idx]=i # gives the index of the secondary index column (for the same person)
        
      else:
        individLUT[who] = i # the primary genotype call column (for this person)
        lastIDX = i + 1
        
    for i in range(0, 9):
      if i > 0:
        out.write("\t") 
      out.write(header[i])

    # just print the IDs once.
    for i in range(9, len(header)):
      if LUT[i] != -1:
        out.write("\t")
        out.write(header[i])
        
    out.write("\n")
  else:
    s = line.rstrip().split("\t")

    outLine = []

    # print out the site information
    for i in range(0, 9):
      outLine.append(s[i])

    for i in range(9, lastIDX):
      if LUT[i] != -1:
        # all info on the genotype call
        rec1 = s[i]
        rec2 = s[ LUT[i] ]

        # and the geno itself
        rec1 = rec1.split(":")
        gt1 = rec1[0]
        
        rec2 = rec2.split(":")
        gt2 = rec2[0]

        # I need to know whether the statistical phasing and physical phasing agree or if they're flipped.
        phaseset = rec2[-1]

        if gt1.find("|") != -1: # the genotype has been phased; we want it
          outLine.append( s[i] )
          
          if gt2.find("|") != -1: # the genotype has been physically phased; we want it
            genos = gt2.split("|")
            or2 = genos[1] + "|" + genos[0] # the other possible phasing orientation


            # chromosome, position, index of who, and which phaseset
            key = (s[0], len(outLine)-1, phaseset)
            if key not in phasesetOrientation:
              phasesetOrientation[key] = [0,0]

            # orientation in statistical phasing is consistent with the physical phasing
            if gt1 == gt2:
              phasesetOrientation[key][0] += 1
#              print(s[0], s[1], gt1, gt2)
              # opposite orientation
              # if opposite count is more, those that are not statistically phased (dropped)
              # but are physically phased need to have their orientations flipped.
            elif gt1 == or2:
              phasesetOrientation[key][1] += 1
            else:
              print("Confusing..", gt1, gt2, or2, key, s, file=sys.stderr)
              
        else: # otherwise we want the original genotype

          # missing data is treated as homozygous reference...
          pipIdx = gt2.find("|")
          if pipIdx < 0:
            if gt2[0] == '.':
              gt2 = '0|0'
            else:
              genos = gt2.split("/")
              if len(genos)>1 and genos[0] == genos[1]:
                gt2 = genos[0] + "|" + genos[1] # homozygotes are already phased...
            
            rec1[0] = gt2
            outLine.append(":".join(rec1) )
          else:

            outLine.append( (s[i], s[LUT[i]], phaseset) )

    if (s[0], s[1]) in finalDict:
      print("Should never be. The same site shows up multiple times. You need to specify bcftools merge -m all", s , sep="\n", file=sys.stderr)
      exit(1)
      
    finalDict[ ( s[0], s[1], s[3], s[4] ) ] = outLine


for (idx, line) in finalDict.items():

  # print out the site information
  for i in range(0, 9):
    if i>0:
      out.write("\t")
    out.write(line[i])
    

  for i in range(9, len(line)):
    out.write("\t")
    rec1 = line[i]
    if type(rec1) is str: # common case. just print it!
      out.write(rec1)
    else:
      # physical phase with no statistical phase!
      (rec1, rec2, phaseset) = rec1

      rec1 = rec1.split(":")
      rec2 = rec2.split(":")
      
      key = (line[0], i, phaseset)
      swap=False
      if key in phasesetOrientation and phasesetOrientation[key][1] > phasesetOrientation[key][0]:
        swap=True

      gt = rec2[0]
      if swap:
        t = gt.split("|")
        gt = t[1] + "|" + t[0]

      rec1[0] = gt

      out.write(":".join(rec1) )
      

  out.write("\n")
  
out.close()



