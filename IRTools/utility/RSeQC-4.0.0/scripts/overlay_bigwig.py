#!/usr/bin/env python
'''Manipulate two bigwig files'''
import os,sys

if sys.version_info[0] != 3:
	print("\nYou are using python" + str(sys.version_info[0]) + '.' + str(sys.version_info[1]) + " This verion of RSeQC needs python3!\n", file=sys.stderr)
	sys.exit()	

from optparse import OptionParser

import pyBigWig
from qcmodule import BED
from qcmodule import twoList
import numpy
__author__ = "Liguo Wang"
__copyright__ = "Copyleft"
__credits__ = []
__license__ = "GPL"
__version__="4.0.0"
__maintainer__ = "Liguo Wang"
__email__ = "wang.liguo@mayo.edu"
__status__ = "Testing"


def main():
	usage="%prog [options]"
	parser = OptionParser(usage,version="%prog " + __version__)
	
	parser.add_option("-i","--bwfile1",action="store",type="string",dest="BigWig_File1",help="One BigWig file.")
	parser.add_option("-j","--bwfile2",action="store",type="string",dest="BigWig_File2",help="Another BigWig file. Both BigWig files should use the same reference genome.")
	parser.add_option("-a","--action",action="store",type="string",dest="action",help='After pairwise align two bigwig files, perform the follow actions (Only select one keyword):"Add" = add signals. "Average" = average signals. "Division"= divide bigwig2 from bigwig1. Add 1 to both bigwig. "Max" = pick the signal that is larger. "Min" = pick the signal that is smaller. "Product" = multiply signals. "Subtract" = subtract signals in 2nd bigwig file from the corresponiding ones in the 1st bigwig file. "geometricMean" = take the geometric mean of signals.')
	parser.add_option("-o","--output",action="store",type="string",dest="output_wig",help="Output wig file")
	parser.add_option("-c","--chunk",action="store",type="int",dest="chunk_size",default=100000,help="Chromosome chunk size. Each chomosome will be cut into samll chunks of this size. Decrease chunk size will save more RAM. default=%default (bp)")
	(options,args)=parser.parse_args()
	
	if not (options.BigWig_File1 and options.BigWig_File2  and options.output_wig):
		parser.print_help()
		sys.exit(0)
	OUT=open(options.output_wig,'w')
	bw1   = pyBigWig.open(options.BigWig_File1)
	bw2   = pyBigWig.open(options.BigWig_File2)
	
	print("Get chromosome sizes from BigWig header ...", file=sys.stderr)
	chrom_sizes = {}
	for chr,size in bw1.chroms().items():
		chrom_sizes[chr] = size
	for chr,size in bw2.chroms().items():
		chrom_sizes[chr] = size
		
	for chr_name, chr_size in list(chrom_sizes.items()):		#iterate each chrom
		print("Processing " + chr_name + " ...", file=sys.stderr)
		OUT.write('variableStep chrom='+chr_name+'\n')
		for interval in BED.tillingBed(chrName = chr_name,chrSize = chr_size,stepSize = options.chunk_size):
			if (bw1.stats(chr_name,interval[1],interval[2] )[0] is None) and (bw2.stats(chr_name,interval[1],interval[2] )[0] is None):
				continue
			coord = interval[1]
			try:
				bw_signal1 = bw1.values(chr_name,interval[1],interval[2])
			except:
				bw_signal1 = numpy.array()
			try:
				bw_signal2 = bw2.values(chr_name,interval[1],interval[2])
			except:
				bw_signal2 = numpy.array()
			if bw_signal1 is None and bw_signal2 is None:
				continue
			if numpy.isnan(numpy.nansum(bw_signal1)) and numpy.isnan(numpy.nansum(bw_signal2)):
				continue
			if len(bw_signal1) == 0 and len(bw_signal2) == 0:
				continue
			bw_signal1 = numpy.nan_to_num( bw_signal1 )
			bw_signal2 = numpy.nan_to_num( bw_signal2 )
		
			call_back = getattr(twoList,options.action)
			for v in call_back(bw_signal1,bw_signal2):
				coord +=1
				if v != 0 : print("%d\t%.2f" % (coord,v), file=OUT)

if __name__=='__main__':
	main()
