### This program is a modified version of Wessim (originally developed from GemSim).

import sys
import random
import bisect
import gzip
import cPickle
import argparse
import math
from time import time, localtime, strftime
import os

inds={'A':0,'T':1,'G':2,'C':3,'N':4,'a':0,'t':1,'g':2,'c':3,'n':4}

def main(argv):
    t0 = time()
    parser = argparse.ArgumentParser(description='sub-wessim: a sub-program for w-Wessim. (NOTE!) Do not run this program. Use "w-Wessim.py" instead. ', prog='w-Wessim-sub', formatter_class=argparse.RawTextHelpFormatter)
    group1 = parser.add_argument_group('Mandatory input files')
    group1.add_argument('-R', metavar = 'FILE', dest='reference', required=True, help='faidx-indexed reference genome FASTA file')
    group1.add_argument('-B', metavar = 'FILE', dest='probeblat', required=True, help='Blat matched probe regions .PSL file')
    group1.add_argument('-S', metavar = 'FILE', dest='syser', required=True, help='Systematic errors file for reference from ReSeq')
    group1.add_argument('-T', metavar = 'FILE', dest='stats', required=True, help='Stats profile for ReSeq')
    group1.add_argument('-N', metavar = 'INT', type=int, dest='readnumber', required=True, help='Total number of read pairs')
    group1.add_argument('-O', metavar = 'INT', dest='outfile', required=True, help='Output file for input to ReSeq')

    group2 = parser.add_argument_group('Optional parameters for exome capture')
    group2.add_argument('-f', metavar = 'INT', type=int, dest='fragsize', required=False, help='mean (f)ragment size, corresponding to insert size between paired ends. [200]', default=200)
    group2.add_argument('-d', metavar = 'INT', type=int, dest='fragsd', required=False, help='standard (d)eviation of fragment size [50]', default=50)
    group2.add_argument('-m', metavar = 'INT', type=int, dest='fragmin', required=False, help='(m)inimum fragment length [read_length + 20]')
    group2.add_argument('-y', metavar = 'PERCENT',type=int, dest='bind', required=False, help='minimum required fraction of probe match to be h(y)bridized [50]', default=50)


    args = parser.parse_args()
    reffile = args.reference
    fref = None
    frefs = []
    metap = []
    metamode = False
    alignfile = args.probeblat

    isize = args.fragsize
    isd = args.fragsd
    imin = args.fragmin
    bind = args.bind

    readnumber = args.readnumber

    if isize < imin:
        print "too small mean fragment size (" + str(isize) + ") compared to minimum length (" + str(imin) + "). Increase it and try again."
        sys.exit(0)

    outfile = args.outfile

    print
    print "-------------------------------------------"
    print "Reference:", reffile
    print "Probematch:", alignfile
    print "Fragment:",isize, "+-", isd, ">", imin
    print "ReSeq stats profile:", args.stats
    print "Systematic errors file:", args.syser
    print "Read number:", readnumber
    print "Output File:", outfile
    print "Job started at:", strftime("%Y-%m-%d %H:%M:%S", localtime())
    print "-------------------------------------------"
    print

    matches=[]
    referencedict={}
    reflist=[]
    tendendict_f={}
    tendendict_r={}
    ratesdict_f={}
    ratesdict_r={}

    #read in reference file into dictionary
    chromo=''
    with open(reffile,'r') as reffi:
        for l in reffi:
            if l.startswith('>'):
                if chromo!='': referencedict[chromo]=''.join(reflist)
                reflist=[]
                chromo=l.strip().split(" ")[0][1:]
                referencedict[l.strip().split(" ")[0]]=''
            else:
                reflist.append(l.strip())
    referencedict[chromo]=''.join(reflist)
    reflist=''

    rat=False
    with open(args.syser,'r') as sysfi:
        for l in sysfi:
            if l.startswith('@') and len(l)<100:
                id=l.strip()[1:l.find(' ')]
                if 'forward' in l:
                    direction='f'
                elif 'reverse' in l:
                    direction='r'
                else:
                    print('ERROR')
            elif l.strip()=='+':
                rat=True
            else:
                if rat==True:
                    if direction=='f':
                        ratesdict_f[id]=l.strip()
                    elif direction=='r':
                        ratesdict_r[id]=l.strip()[::-1]
                    rat=False
                else:
                    if direction=='f':
                        tendendict_f[id]=l.strip()
                    elif direction=='r':
                        tendendict_r[id]=l.strip()[::-1]

    #Read in align and probe file

    with open(alignfile, 'r') as f2:
        for i in range(0, 6): ### Ignore first 5 lines of psl file (header)
            line = f2.readline()
        while line:
            values = line.split("\t")
            if len(values)<16:
                line = f1.readline()
                line = f1.readline()
                continue
            qgapsize = int(values[5])
            tgapsize = int(values[7])
            strand= values[8]
            pslchrom = values[13]
            pslstart = values[15]
            pslend= values[16]
            if qgapsize > 2 or tgapsize > 2:
                line = f2.readline()
                continue
            matches.append((pslchrom, int(pslstart), int(pslend),strand))
            line = f2.readline()

    ### Generate!
    with open(outfile,'w+') as ofile:
        i = 0
        while i < readnumber:
            match = pickonematch(matches)
            fragment = getFragment(match, isize, isd, imin, bind)
            fragment_start = int(fragment[1])
            if fragment_start < 0:
                continue
            ref = getSequence(referencedict, fragment)
            if len(ref)<imin:
                continue
            ref_r= comp(ref)[::-1][:180]
            ref_f=ref[:180]
            tendenseq_f=getSequence(tendendict_f, fragment)[:180]
            ratesseq_f=getSequence(ratesdict_f, fragment)[:180]
            tendenseq_r=getSequence(tendendict_r, fragment)[::-1][:180]
            ratesseq_r=getSequence(ratesdict_r, fragment)[::-1][:180]
            fragment_length=str(int(fragment[2])-int(fragment[1]))
            if fragment[3]=='+':
                ofile.write('>'+str(i+1)+' 1;'+fragment_length+';'+tendenseq_f+';'+ratesseq_f+'\n'+ref_f+'\n')
                ofile.write('>'+str(i+1)+' 2;'+fragment_length+';'+tendenseq_r+';'+ratesseq_r+'\n'+ref_r+'\n')
            elif fragment[3]=='-':
                ofile.write('>'+str(i+1)+' 1;'+fragment_length+';'+tendenseq_r+';'+ratesseq_r+'\n'+ref_r+'\n')
                ofile.write('>'+str(i+1)+' 2;'+fragment_length+';'+tendenseq_f+';'+ratesseq_f+'\n'+ref_f+'\n')

            i+=1
            if i+1 % 1000000 == 0:
                t1 = time()
                print str(i+1) + " reads have been processed... in %f secs" % (t1-t0)

    t1 = time()
    print "Done processing " + str(readnumber) + " reads in %f secs" % (t1 - t0)

def pickonematch(matches):
	match = random.choice(matches)
	return match

def getSequence(dict, fragment):
	chrom = fragment[0]
	start = int(fragment[1])
	end = int(fragment[2])
	seq = dict[chrom][start:end]
	return seq

def getFragment(match, mu, sigma, lower, bind):
	ins = getInsertLength(mu, sigma, lower)
	pickedfragment = pickFragment(match, ins, bind)
	return pickedfragment

def getInsertLength(mu, sigma, lower):
	while True:
		length = int(random.gauss(mu, sigma))
		if length >= lower:
			return length

def pickFragment(match, ins, bind):
	probechrom = match[0]
	probestart = int(match[1])
	probeend = int(match[2])
	probelength = probeend - probestart
	minimummatch = int(probelength*bind/100)
	overlap = int(random.triangular(minimummatch, probelength, probelength))
	margin = max(ins - overlap, 0)
	rangestart = probestart - margin
	rangeend = probeend + margin
	seqstart = random.randint(rangestart, rangeend - ins)
	strand=match[3]
	return probechrom, seqstart, seqstart + ins, strand

def comp(sequence):
	""" complements a sequence, preserving case. Function imported from GemSim"""
	d={'A':'T','T':'A','C':'G','G':'C','a':'t','t':'a','c':'g','g':'c','N':'N','n':'n'}
	cSeq=''
	for s in sequence:
		if s in d.keys():
			cSeq+=d[s]
		else:
			cSeq+='N'
	return cSeq

if __name__=="__main__":
	main(sys.argv[1:])
	sys.exit(0)

