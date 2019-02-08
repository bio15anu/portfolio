#!/usr/bin/python3

'''

Title: Reciprocal_BLAST.py
Date: 2017-01-22
Author: Adam Nunn

Description:
  This program takes an input file containing the standard output from a
  single blastp search with query fasta file containing sequences from N
  organisms against a database of the same fasta file. In the command line,
  the user will define sequence identifiers (common to all sequences) for
  N-1 of the organisms. The program will identify reciprocal best hit pairs
  for all pairwise comparisons within N, then output only those queries
  that share reciprocal best hits in all N organisms to "ortho_pairs+.txt"
  

List of functions:
  open()
  .startswith()
  .split()
  any()
  .append()
  len()
  all()
  print()
  
  
Procedure:
  1. Locate the current query identifier and define variables for later use.
  2. Locate the next best hit for the query in the database, setup for 3) and 4)
  3. If the query is from an organism referenced by the user, determine if hit is
      also referenced or not, then append to dictionary entry for query only if
      there does not exist an entry for that organism already.
  4. If the query is not from an organism referenced by the user, only consider
      hits that are from referenced organisms, then append to dictionary entry for
      query only if there does not exist an entry for that organism already.
  5. Iterate through the dictionary of best hit pairs for only the FIRST organism
      referenced by user (ie. REF(1)) and output to file only if reciprocal best
      hits in ALL other organisms AND those organisms also share reciprocal best
      hits among themselves.

Usage:
     blastp -query combined.faa -db combined.faa -evalue 1e-10 -out combined.blastp

    ./blastReciprocal+.py INPUT_FILE REF(1) ... REF(N-1)  # for N organisms
eg. ./blastReciprocal+.py combined.blastp Paxin Phchr2 

'''

import sys

list_of_refs = sys.argv[2:] # eg. ["Paxin1", "Phchr2"]
pairDict = {} # eg. {"Paxin1": ["Phchr2","yeast"], "Phchr2": ["Paxin1","yeast"], "yeast": ["Paxin1","Phchr2"]

foundQuery = False

with open(sys.argv[1], 'r') as infile, open("ortho_pairs+.txt", 'w') as outfile:
  
 for line in infile: # parse the input blastp file
  line = line.rstrip()
  
  ### 1) Locate the query identifier
  if line.startswith("Query="): # find the query
   line = line.split(" ")
   query = line[1]
   pairDict[query] = [] # make an entry for the query in pairDict
   
   foundQuery = True # query found, now look for hits
   foundElse = False # hit for non-referenced organism not found yet
  
  ### 2) Locate the hit identifier
  elif (foundQuery == True) and line.startswith("  "): # iterate through hits under each query
   line = line.split(" ")
   hit = line[2]
   
   ### 3) Generate entry for query if from organism referenced by the user
   if any(ref in query for ref in list_of_refs): # CONDITIONAL IF query is in list_of_refs eg. ["Paxin1", "Phchr2", etc.]
    for i in list_of_refs: # iterate through the user-specified references eg. ["Paxin1", "Phchr2", etc.]
     # rule out hits from same organism eg. "Paxin1 vs Paxin1", and make sure there isn't already a hit in pairDict
     if hit and (i in hit) and not (i in query) and not any(i in j for j in pairDict[query]): 
      pairDict[query].append(hit)
      break
     # rule out list_of_refs from hit, return only eg. "Paxin1 vs yeast" or "Phchr2 vs yeast"
     elif hit and (foundElse == False) and not any(k in hit for k in list_of_refs):
      pairDict[query].append(hit)
      foundElse = True # make sure to only capture the single best hit for non-referenced organisms
      break
   
   ### 4) Generate entry for query if from organism not referenced by the user
   else: # CONDITIONAL ELSE query is not in list_of_refs
    for refs in list_of_refs: # iterate through list of references eg. ["Paxin1", "Phchr2", etc.]
     if hit and (refs in hit) and not any(refs in elems for elems in pairDict[query]):
      pairDict[query].append(hit)
      break
  
  elif line.startswith("*") or line.startswith(">") and (foundQuery == True):
   foundQuery = False
   if (len(pairDict[query]) < len(list_of_refs)): del pairDict[query]
 
 ### 5) Iterate through reciprocal best hit dictionary and output to file
 for seq in pairDict: # eg. {"Paxin1": ["Phchr2","yeast"], "Phchr2": ["Paxin1","yeast"], "yeast": ["Paxin1","Phchr2"]
  # only iterate through one set of organisms, make sure all reciprocal organisms have complete sets
  if (list_of_refs[0] in seq) and all(keys in pairDict for keys in pairDict[seq]):
   # check for reciprocal best hits in all cases
   if all(seq in pairDict[x] for x in pairDict[seq]) and \
       all(a in pairDict[b] for a in pairDict[seq] for b in pairDict[seq] if a!=b):
    print("{}\t{}".format(seq, "\t".join(pairDict[seq])), file=outfile)

