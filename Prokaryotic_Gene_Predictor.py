#!/usr/bin/python3

'''
Title: Prokaryotic_Gene_Predictor.py
Date: 2017-02-11
Author: Adam Nunn

Description:
  This program takes an input fasta file containing a single PROKARYOTIC
  genome sequence, together with user-defined search parameters and scoring
  system, and outputs predicted genes in the forward and reverse complement
  direction, ranked according to log-likelihood score. The search algorithm
  is based on finiding the -35 promotor sequence TTGACA, followed by the
  -10 promotor sequence TATAAT, followed by first incidence of a known
  start codon sequence. The gene transcript is then given until the first
  incidence of a known stop codon sequence. The log-likelihood score system
  is calculated by taking the relative frequency of all possible kmers in
  the genome and comparing to the kmers in the predicted gene regions. The
  score for each gene sequence is then given as the sum of kmer scores for
  that gene.

List of functions:
  parser.add_argument()
  parser.parse_args()
  os.makedirs()
  int()
  open()
  .startswith()
  .rstrip()
  .upper()
  translate()
  str.maketrans()
  print()
  quit()
  list()
  range()
  len()
  .append()
  .join()
  re.findall()
  enumerate()
  dict()
  zip()
  .count()
  os.path.join()
  .format()
  sorted()
  set()
  
  
Procedure:
  1. Designate user-defined variables and define empty variables for later use.
  2. Interpret the input fasta sequence. Concatenate interleaved sequence
      lines to a single line, then find the reverse complement. ABORT the
      program if more than one sequence is found in the input file.
      Search through the sequence, base by base, for -35 -10 and start
      codon sequences. Once found, save transcript and start pos to lists
      ending at the first incidence of stop codon TAA, TAG or TGA. Lists
      are defined based on starting position and therefore Open Reading Frame.
  3. Output to STDOUT basic log information such as overall sequence length,
      GC content and number of genes predicted for each ORF. Output to new
      files "forward.fna" and "reverse.fna" all genes predicted on the
      forward strand and the reverse strand (with meta info), respectively.
  4. Define the scoring system by calculating frequencies of all possible
      k-mers across the entire input sequence, then calculate the log-likelihood
      of k-mer frequencies as a proportion of the total frequency between the
      overall sequence and those found in the gene sequences. Individual gene
      scores are then calculated by summation of kmer log-likelihood scores.
  5. Output fasta files contained all predicted genes, one unranked and one
      ranked according to overall log-likelihood scores. Output k-mer frequency
      data and log/log-likelihood scores for complete sequence and gene sequences.


Usage:

    ./Prokaryotic_Gene_Predictor.py INPUT_FILE [-m number_of_mismatches] [-p -35_window] \
    [-s -10_window] [-g min_gene_length] [-k kmer_length] [-S start_codons] OUTPUT

eg. ./Prokaryotic_Gene_Predictor.py -m2 -p25 -s10 -g100 -k6 -S1 input.fna projectname 

'''

import argparse
import re
import os

# define argparse
usage = '''  This program takes an input fasta file containing a single PROKARYOTIC genome sequence, together with user-defined search parameters and scoring system, and outputs predicted genes in the forward and reverse complement direction, ranked according to log-likelihood score. The search algorithm is based on finiding the -35 promotor sequence TTGACA, followed by the -10 promotor sequence TATAAT, followed by first incidence of a known start codon sequence. The gene transcript is then given until the first incidence of a known stop codon sequence. The log-likelihood score system is calculated by taking the relative frequency of all possible kmers in the genome and comparing to the kmers in the predicted gene regions. The score for each gene sequence is then given as the sum of kmer scores for that gene.  '''

parser = argparse.ArgumentParser(description=usage)
parser.add_argument('-m', '--mismatch', metavar='', help='Number of allowed mismatches in promotor sequences (default = 2)', type=int, default=2)
parser.add_argument('-p', '--promotor_window', metavar='', help='Window size for identifying -10 promotor sequence after successful identification of -35 sequence (default = 25)',  type=int, default=25)
parser.add_argument('-s', '--start_window', metavar='', help='Window size for identifying start codon after -10 promotor (default = 10)',  type=int, default=10)
parser.add_argument('-g', '--min_gene_length', metavar='', help='Minimum length for gene products included in output (default = 100)',  type=int, default=100)
parser.add_argument('-k', '--kmer_length', metavar='', help='K-mer length for computing test statistics (default = 6)',  type=int, default=6)
parser.add_argument('-S', '--start_codons', metavar='', choices=[1, 2, 3], help='Number of start codon variations to include (1: ATG; 2: ATG, GTG; 3: ATG, GTG, TTG) (default = 1)',  type=int, default=1)
parser.add_argument('infile', metavar='INPUT_FILE', help='Path for INPUT fasta file to analyse')
parser.add_argument('outfile', metavar='PROJECT_NAME', help='Project name for OUTPUT directory')

args = parser.parse_args()
parser.parse_args()


### 1) SETUP THE ENVIRONMENT

# create output directory
dir_path = args.outfile
os.makedirs(dir_path)

# define variables
if args.start_codons == 1: start = ('ATG',)
elif args.start_codons == 2: start = ('ATG','GTG')
else: start = ('ATG','GTG','TTG')

stop = ('TAA', 'TAG', 'TGA')
concatenate = ''
minus35 = ['T', 'T', 'G', 'A', 'C', 'A']
minus10 = ['T', 'A', 'T', 'A', 'A', 'T']
ORF1, ORF2, ORF3, ORF4, ORF5, ORF6 = [], [], [], [], [], []
ORF1pos, ORF2pos, ORF3pos, ORF4pos, ORF5pos, ORF6pos = [], [], [], [], [], []

# define user arguments
range35 = int(args.promotor_window)+5
range10 = int(args.start_window)
minlen = int(args.min_gene_length)


### 2) INITIAL GENE PREDICTION

# interpret the input file
with open(args.infile, 'r') as infile:
 oneseq = False
 for line in infile:
  if not line.startswith('>'):
   line = line.rstrip().upper()
   concatenate += line
  else:
   if oneseq == False: oneseq = True
   else: print('ERROR: Invalid File: More than one sequence present in file'); quit()

# reverse complement
comcat = concatenate.translate(str.maketrans('ACGT', 'TGCA'))
revcat = comcat[::-1]

# prepare strands for analysis
forw = list(concatenate)
revr = list(revcat)
ForRev = (forw, revr)

# search for sequences in strands
for listline in ForRev:
 found35 = False
 found10 = False
 for n in range(len(listline)-(range35+range10)): # iterate through bases in sequence
  if found35 == False and found10 == False: # SEARCH FOR -35 PROMOTOR
   mismatch35 = 0
   match35 = listline[n:(n+6)] # take sequence length 6 from current base
   for x in range(6):
    if match35[x] != minus35[x]: mismatch35 += 1 # iterate through potential -35 match base by base
   if mismatch35 <= int(args.mismatch): found35 = True; # FOUND the -35 promotor!
  elif found35 == True and found10 == False: # SEARCH FOR -10 PROMOTOR
   for n2 in range((n+5), (n+range35)): # search the next 30 bases from end of -35 promotor
    mismatch10 = 0
    match10 = listline[n2:(n2+6)] # take sequence length 6 from current base
    for x2 in range(6):
     if match10[x2] != minus10[x2]: mismatch10 += 1 # iterate through potential -10 match base by base
    if mismatch10 <= int(args.mismatch):
     found10 = True
     Pos10 = n2+5
     break # FOUND the -10 promotor!
   found35 = False
  elif found10 == True: #SEARCH FOR ATG CODON
   startmer = [] # NEW CODE
   for n3 in range((Pos10), (Pos10+range10)): # search the next 15 bases from end of -10 promotor
    startmer.append(listline[n3:(n3+3)]) # [[A,T,G],[T,G,G],[G,G,A],etc.]
   mismatchATG = 0
   for x3 in range(len(startmer)):
    mismatchATG += 1
    if "".join(startmer[x3]) in start:
     StartPos = (mismatchATG + Pos10)-1
     seqMatch = re.findall('...', (''.join(listline)[StartPos::])) # break the remaining sequence into codons
     for i,codon in enumerate(seqMatch): #list all codons together with index number
      if codon in stop: #is current codon a STOP codon?
       closeread = i+1 #note the stop index
       newSeq = ''.join(seqMatch[0:closeread])
       break #end the seqMatch search
     if (listline == forw) and (StartPos % 3 == 0) and (newSeq not in ORF1): #assign the read to an ORF list
      ORF1.append(newSeq) #add the sequence to ORF forward 1 list
      ORF1pos.append(StartPos)
     elif (listline == forw) and (StartPos % 3 == 1) and (newSeq not in ORF2):
      ORF2.append(newSeq) #add the sequence to ORF forward 2 list
      ORF2pos.append(StartPos)
     elif (listline == forw) and (StartPos % 3 == 2) and (newSeq not in ORF3):
      ORF3.append(newSeq) #add the sequence to ORF forward 3 list
      ORF3pos.append(StartPos)
     elif (listline == revr) and (StartPos % 3 == 0) and (newSeq not in ORF4): #assign the read to an ORF list
      ORF4.append(newSeq) #add the sequence to ORF reverse 1 list
      ORF4pos.append(StartPos)
     elif (listline == revr) and (StartPos % 3 == 1) and (newSeq not in ORF5):
      ORF5.append(newSeq) #add the sequence to ORF reverse 2 list
      ORF5pos.append(StartPos)
     elif (listline == revr) and (StartPos % 3 == 2) and (newSeq not in ORF6):
      ORF6.append(newSeq) #add the sequence to ORF reverse 3 list
      ORF6pos.append(StartPos)
     found10 = False # reset promotor search
     break #break the search and start entire process again from next base
   found10 = False # reset promotor search

#dictionaries for sequences and starting positions
dictorf1 = dict(zip(ORF1, ORF1pos))
dictorf2 = dict(zip(ORF2, ORF2pos))
dictorf3 = dict(zip(ORF3, ORF3pos))
dictorf4 = dict(zip(ORF4, ORF4pos))
dictorf5 = dict(zip(ORF5, ORF5pos))
dictorf6 = dict(zip(ORF6, ORF6pos))


### 3) OUTPUT LOG AND UNFILTERED SEQUENCES

# print log to STDOUT
content = []
for stdbp in (concatenate, revcat):
 stdbp0 = stdbp.upper()
 stdbp1 = stdbp[::3].upper()
 stdbp2 = stdbp[1::3].upper()
 stdbp3 = stdbp[2::3].upper()
 for GC in (stdbp0, stdbp1, stdbp2, stdbp3):
  G = GC.count("G")
  C = GC.count("C")
  A = GC.count("A")
  T = GC.count("T")
  content.append((G + C)/(G + C + T + A))

print("\n############\nSequence Length: {} bp\nGC Content: {}\n\n############".format(len(concatenate),content[0]))

# print the unscored output
with open(os.path.join(dir_path, 'forward.fna'), 'w') as outfor, open(os.path.join(dir_path, 'reverse.fna'), 'w') as outrev:
 ORFS = (ORF1, ORF2, ORF3, ORF4, ORF5, ORF6)
 dictorf = (dictorf1, dictorf2, dictorf3, dictorf4, dictorf5, dictorf6)
 ORFcount = 0
 tempCount = 0
 print("\nGENES PREDICTED (Forward Strand)")
 for orf in range(3):
  for seq in ORFS[orf]:
   if len(seq) >= minlen:
    ORFcount += 1
    print('>{}:GENE:{}_1\tORF= {}\tStartPos= {}\tStrand= +\n{}'.format(args.outfile, ORFcount, orf, dictorf[orf][seq], seq), file=outfor)
  # print log to STDOUT
  print("ORF{} genes predicted: {}\nGC Content: {}\n".format(orf+1, ORFcount-tempCount, content[orf+1]))
  tempCount = ORFcount
 ORFcount = 0
 tempCount = 0
 print("\nGENES PREDICTED (Reverse Strand)")
 for orf in range(3,6):
  for seq in ORFS[orf]:
   if len(seq) >= minlen:
    ORFcount += 1
    print('>{}:GENE:{}_2\tORF= {}\tStartPos= {}\tStrand= -\n{}'.format(args.outfile,ORFcount, orf-3, dictorf[orf][seq], seq), file=outrev)
  # print log to STDOUT
  print("ORF{} genes predicted: {}\nGC Content: {}\n".format(orf-2, ORFcount-tempCount, content[orf+2]))
  tempCount = ORFcount

print("############\n")


### 4) DEFINE AND PROCESS SCORING SYSTEM

import math
os.makedirs(dir_path + "/statistics")

# compute overall kmers
kmer = []
for sequences in (concatenate, revcat):
 for k in range(len(sequences)-(args.kmer_length-1)):
  kmer.append(sequences[k:k+(args.kmer_length)])

uniqWords = sorted(set(kmer))
logscore=[]

for word in uniqWords:
 logscore.append(math.log2(kmer.count(word)/len(kmer)))
 
overall = dict(zip(uniqWords, logscore))

# concatenate sequence meta data into list of lists
ORFS = (ORF1, ORF2, ORF3, ORF4, ORF5, ORF6)
dictorf = (dictorf1, dictorf2, dictorf3, dictorf4, dictorf5, dictorf6)

forwardLIST = []
reverseLIST = []

ORFcount = 0
for orf in range(3):
 for seq in ORFS[orf]:
  if len(seq) >= minlen:
   metaLIST = []
   ORFcount += 1
   metaLIST.append(str(ORFcount) + "_1")
   metaLIST.append(orf)
   metaLIST.append(dictorf[orf][seq])
   metaLIST.append("+")
   metaLIST.append(seq)
   forwardLIST.append(metaLIST)

ORFcount = 0
for orf in range(3,6):
 for seq in ORFS[orf]:
  if len(seq) >= minlen:
   metaLIST = []
   ORFcount += 1
   metaLIST.append(str(ORFcount) + "_2")
   metaLIST.append(orf-3)
   metaLIST.append(dictorf[orf][seq])
   metaLIST.append("-")
   metaLIST.append(seq)
   reverseLIST.append(metaLIST)

finalLIST = forwardLIST+reverseLIST
iterLIST = (forwardLIST, reverseLIST)

gene_forward = []
gene_reverse = []
geneLIST = (gene_forward, gene_reverse)

# retrieve co-ordinates for gene coding segments
for iter0 in range(2):
 gene_cover = []
 for iter1 in range(len(iterLIST[iter0])):
  x_loop = sorted(iterLIST[iter0], key=lambda x: x[2])
  x_finish = (x_loop[iter1][2]+len(x_loop[iter1][4]))
  co_ordinate = [x_loop[iter1][2], x_finish]
  gene_cover.append(co_ordinate)
 xcoor = gene_cover[0][0]
 ycoor = gene_cover[0][1]
 for coor in gene_cover:
  if (coor[0] >= xcoor) and (coor[0] < ycoor) and (coor[1] <= ycoor): continue;
  elif (coor[0] > xcoor) and (coor[0] < ycoor) and (coor[1] > ycoor): ycoor = coor[1];
  else: geneLIST[iter0].append([xcoor, ycoor]); xcoor = coor[0]; ycoor = coor[1];
 geneLIST[iter0].append([xcoor, ycoor])

# retrieve slices for gene coding segments from co-ordinates
slices = []
for seqs in range(2):
 for coords in geneLIST[seqs]:
  sequence = (concatenate, revcat)[seqs][coords[0]:coords[1]]
  slices.append(sequence)

# compute kmers and scores for predicted genes
kmer_score = []
gene_kmer = []

# compute total kmers for predicted genes
for genes in range(len(slices)):
 for gkmers in range(len(slices[genes])-(args.kmer_length-1)): 
  gene_kmer.append(slices[genes][gkmers:gkmers+(args.kmer_length)]) 

uniqGWords = sorted(set(gene_kmer))

# construct overall scoring system
for gword in uniqGWords:
 score = math.log2(gene_kmer.count(gword)/len(gene_kmer))
 kmer_score.append(score-overall[gword])

dictGWords = dict(zip(uniqGWords, kmer_score))

# scoring individual genes
gene_score = []
for genes2 in range(len(finalLIST)):
 sumscore = 0
 for ikmers in range(len(finalLIST[genes2][4])-(args.kmer_length-1)):
  current = finalLIST[genes2][4][ikmers:ikmers+(args.kmer_length)]
  sumscore += dictGWords[current]
 gene_score.append(sumscore)


### 5) OUTPUT RANKED GENES AND STATISTICS META DATA

# output unranked fasta file
with open(os.path.join(dir_path + '/statistics', 'unranked.fasta'), 'w') as unranked:
 element = 0
 for item in finalLIST:
  print(">{}:GENE:{}\tORF= {}\tStartPos= {}\tEndPos ={}\tStrand= {}\tScore= {}\n{}".format(args.outfile,item[0],item[1],item[2],(item[2]+len(item[4])),item[3],gene_score[element],item[4]), file=unranked)
  element += 1

# output ranked fasta file
with open(os.path.join(dir_path + '/statistics', 'ranked.fasta'), 'w') as ranked:
 for elem in sorted(range(len(gene_score)), key=gene_score.__getitem__, reverse=True):
  print(">{}:GENE:{}\tORF= {}\tStartPos= {}\tEndPos= {}\tStrand= {}\tScore= {}\n{}".format(args.outfile,finalLIST[elem][0],finalLIST[elem][1],finalLIST[elem][2],(finalLIST[elem][2]+len(finalLIST[elem][4])),finalLIST[elem][3],gene_score[elem],finalLIST[elem][4]), file=ranked)

# output overall kmer freq
with open(os.path.join(dir_path + '/statistics', 'kmer.complete.freq'), 'w') as stats:
 for word in uniqWords: 
  print("{}\t{}\t{}".format(kmer.count(word), word, overall[word]), file=stats)

# output gene kmer freq
with open(os.path.join(dir_path + '/statistics', 'kmer.genes.freq'), 'w') as gstats:
 for gword in uniqGWords: 
  print("{}\t{}\t{}".format(gene_kmer.count(gword), gword, dictGWords[gword]), file=gstats)
