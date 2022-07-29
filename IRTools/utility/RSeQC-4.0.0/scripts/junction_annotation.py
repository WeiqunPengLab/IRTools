#!/usr/bin/env python
'''
Annotate splicing reads against gene model in two levels: reads level and  juncion level.
Note:
1) A read, especially long read, can be spliced 2 or more times
2) Multiple splicing reads spanning the same intron can be consolidated into one splicing junction.
'''

#import built-in modules
import os,sys
if sys.version_info[0] != 3:
	print("\nYou are using python" + str(sys.version_info[0]) + '.' + str(sys.version_info[1]) + " This verion of RSeQC needs python3!\n", file=sys.stderr)
	sys.exit()	

import re
import string
from optparse import OptionParser
import warnings
import string
import collections
import math
from time import strftime
import subprocess

#import my own modules
from qcmodule import SAM
#changes to the paths

#changing history to this module


__author__ = "Liguo Wang"
__copyright__ = "Copyleft"
__credits__ = []
__license__ = "GPL"
__version__="4.0.0"
__maintainer__ = "Liguo Wang"
__email__ = "wang.liguo@mayo.edu"
__status__ = "Production"


def printlog (mesg):
        '''print progress into stderr and log file'''
        mesg="@ " + strftime("%Y-%m-%d %H:%M:%S") + ": " + mesg
        LOG=open('class.log','a')
        print(mesg, file=sys.stderr)
        print(mesg, file=LOG)

def generate_bed12(infile,size=1):
	''' 
	infile: input file. eg: chrX    66766604        66788677        348     partial_novel
	size: the block size representing exons
	'''
	
	outfile = infile.replace('.xls','.bed')
	OUT = open(outfile,'w')
	for line in open(infile,'r'):
		if line.startswith('chrom'):continue
		line = line.strip()
		f = line.split()
		if len(f) != 5:continue
		chrom = f[0]
		start = int(f[1]) - size
		end = int(f[2]) + size
		score = int(f[3])
		strand = '.'
		name = f[4]
		thick_st = start
		thick_end = end
		if name == 'annotated':
			color = '205,0,0'
		elif name == 'partial_novel':
			color = '0,205,0'
		elif name == 'complete_novel':
			color = '0,0,205'
		else:
			color = '0,0,0'
		blockCount = 2
		blockSizes = ','.join((str(size),str(size)))
		blockStarts = '0,' + str(end-size-start)
		print('\t'.join( [str(i) for i in [chrom, start, end, name, score, strand, thick_st, thick_end,color, blockCount, blockSizes,blockStarts]]), file=OUT)
	OUT.close()
		
def generate_interact(infile, bam_file, size=1):
	''' 
	infile: input file. eg: chrX    66766604        66788677        348     partial_novel
	size: the block size representing exons
	'''
	
	outfile = infile.replace('.xls','.Interact.bed')
	OUT = open(outfile,'w')
	print ('track type=interact name="Splice junctions" description="Splice junctions detected from %s" maxHeightPixels=200:200:50 visibility=full' % bam_file, file=OUT)
	for line in open(infile,'r'):
		if line.startswith('chrom'):continue
		line = line.strip()
		f = line.split()
		if len(f) != 5:continue
		chrom = f[0]
		chromStart = int(f[1]) - size
		chromEnd = int(f[2]) + size
		name1 = f[4]
		if name1 == 'annotated':
			color = '205,0,0'
		elif name1 == 'partial_novel':
			color = '0,205,0'
		elif name1 == 'complete_novel':
			color = '0,0,205'
		else:
			color = '0,0,0'

		name = chrom + ":" + str(chromStart) + '-' + str(chromEnd) + "_" + name1
		score = int(f[3])
		value = float(score)
		exp = 'RNAseq_junction'
		
		sourceChrom = chrom
		sourceStart = chromStart
		sourceEnd = chromStart + size
		sourceName = sourceChrom + ":" + str(sourceStart) + '-' + str(sourceEnd)
		sourceStrand = '.'
		
		targetChrom = chrom
		targetStart = chromEnd - size
		targetEnd = chromEnd
		targetName = targetChrom + ":" + str(targetStart) + '-' + str(targetEnd)
		targetStrand = '.'
		print ('\t'.join([str(i) for i in [chrom, chromStart, chromEnd, name, score, value, exp, color, sourceChrom, sourceStart, sourceEnd, sourceName, sourceStrand, targetChrom, targetStart, targetEnd, targetName, targetStrand ]]), file=OUT)
	OUT.close()

def main():
	usage="%prog [options]" + '\n' + __doc__ + "\n"
	parser = OptionParser(usage,version="%prog " + __version__)
	parser.add_option("-i","--input-file",action="store",type="string",dest="input_file",help="Alignment file in BAM or SAM format.")
	parser.add_option("-r","--refgene",action="store",type="string",dest="ref_gene_model",help="Reference gene model in bed format. This file is better to be a pooled gene model as it will be used to annotate splicing junctions [required]")
	parser.add_option("-o","--out-prefix",action="store",type="string",dest="output_prefix",help="Prefix of output files(s). [required]")
	parser.add_option("-m","--min-intron",action="store",type="int",dest="min_intron",default=50, help="Minimum intron length (bp). default=%default [optional]")
	parser.add_option("-q","--mapq",action="store",type="int",dest="map_qual",default=30,help="Minimum mapping quality (phred scaled) for an alignment to be considered as \"uniquely mapped\". default=%default")

	(options,args)=parser.parse_args()
		
	if not (options.output_prefix and options.input_file and options.ref_gene_model):
		parser.print_help()
		sys.exit(0)
	if not os.path.exists(options.ref_gene_model):
		print('\n\n' + options.ref_gene_model + " does NOT exists" + '\n', file=sys.stderr)
		sys.exit(0)
	if os.path.exists(options.input_file):
		obj = SAM.ParseBAM(options.input_file)
		obj.annotate_junction(outfile=options.output_prefix,refgene=options.ref_gene_model,min_intron=options.min_intron, q_cut = options.map_qual)
		try:
			subprocess.call("Rscript " + options.output_prefix + '.junction_plot.r', shell=True)
		except:
			print("Cannot generate pdf file from " + '.junction_plot.r', file=sys.stderr)
			pass
	else:
		print('\n\n' + options.input_file + " does NOT exists" + '\n', file=sys.stderr)
		sys.exit(0)
	try:
		print ('Create BED file ...', file=sys.stderr)
		generate_bed12(options.output_prefix + '.junction.xls')	
	except:
		pass	
	try:
		print ('Create Interact file ...', file=sys.stderr)
		generate_interact(options.output_prefix + '.junction.xls', options.input_file)	
	except:
		pass	


if __name__ == '__main__':
        main()
