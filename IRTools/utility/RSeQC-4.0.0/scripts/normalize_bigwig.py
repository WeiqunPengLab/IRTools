#!/usr/bin/env python
'''Normalize bigwig signal to fixed wigsum (equivelent to total reads). Output wiggle file'''
import os,sys
import collections
from operator import itemgetter
from itertools import groupby

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
	
	parser.add_option("-i","--bwfile",action="store",type="string",dest="BigWig_File",help="Input BigWig file. [required]")
	parser.add_option("-o","--output",action="store",type="string",dest="output_wig",help="Output wig file. [required]")
	parser.add_option("-t","--wigsum",action="store",type="int",dest="total_wigsum",default=100000000,help="Specified wigsum. 100000000 equals to coverage of 1 million 100nt reads. default=%default  [optional]")
	parser.add_option("-r","--refgene",action="store",type="string",dest="refgene_bed",help="Reference gene model in bed format. [optional]")	
	parser.add_option("-c","--chunk",action="store",type="int",dest="chunk_size",default=500000,help="Chromosome chunk size. Each chomosome will be cut into samll chunks of this size. Decrease chunk size will save more RAM. default=%default (bp) [optional]")
	parser.add_option("-f","--format",action="store",type="string",dest="out_format",default="bgr",help="Output format. either \"wig\" or \"bgr\". \"bgr\" save disk space but make program slower. default=%default")
	(options,args)=parser.parse_args()
	
	if not (options.BigWig_File and options.output_wig):
		parser.print_help()
		sys.exit(0)

	OUT=open(options.output_wig,'w')
	bw = pyBigWig.open(options.BigWig_File)
	
	if bw.isBigWig():
		pass
	else:
		print("%s is not a bigwig file!" % options.BigWig_File,  file=sys.stderr)
		sys.exit(0)
	
	print("Get chromosome sizes from BigWig header ...", file=sys.stderr)
	chrom_sizes = {}
	for chr,size in bw.chroms().items():
		chrom_sizes[chr] = size
	
	exons=[]
	WIG_SUM=0.0
	if (options.refgene_bed):	
		print("Extract exons from " + options.refgene_bed, file=sys.stderr)
		obj = BED.ParseBED(options.refgene_bed)
		exons = obj.getExon()
		print("Merge overlapping exons ...", file=sys.stderr)
		exons = BED.unionBed3(exons)
		print("Calculate wigsum covered by " + options.refgene_bed + ' only', file=sys.stderr)
		for chrom,st,end in exons:
			if bw.stats(chrom, st, end)[0] is None:
				continue
			bw_signal = bw.values(chrom,st,end)
			tmp = numpy.nansum(bw_signal)			#nan will be ignored. but if all items are 'nan', the result summay is 'nan' NOT 0
			if numpy.isnan(tmp):continue	
			WIG_SUM += tmp
		print("Total wigsum is %.2f\n" % WIG_SUM, file=sys.stderr)
	else:
		print("Calculate wigsum from " + options.BigWig_File, file=sys.stderr)
		for chr_name, chr_size in list(chrom_sizes.items()):		#iterate each chrom
			if bw.stats(chr_name, 0, chr_size)[0] is None:
				print("Skip " + chr_name + "!", file=sys.stderr)
				continue

			print("Processing " + chr_name + " ...", file=sys.stderr)	
			for interval in BED.tillingBed(chrName = chr_name,chrSize = chr_size,stepSize = options.chunk_size):
				if bw.stats(interval[0],interval[1],interval[2])[0] is None:
					continue
				bw_signal  = bw.values(interval[0],interval[1],interval[2])
				tmp = numpy.nansum(bw_signal)
				if numpy.isnan(tmp):continue
				WIG_SUM += tmp
		print("\nTotal wigsum is %.2f\n" % WIG_SUM, file=sys.stderr)
	
	try:
		weight = options.total_wigsum/WIG_SUM
	except:
		"Error, WIG_SUM cannot be 0"
		sys.exit(1)

	#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	print("Normalizing bigwig file ...", file=sys.stderr)
	for chr_name, chr_size in list(chrom_sizes.items()):          #iterate each chrom
		
		if bw.stats(chr_name, 0, chr_size)[0] is None:
			print("Skip " + chr_name + "!", file=sys.stderr)
			continue
		
		if options.out_format.upper() == "WIG":
			print("Writing " + chr_name + " ...", file=sys.stderr)
			OUT.write('variableStep chrom='+chr_name+'\n')
			for interval in BED.tillingBed(chrName = chr_name,chrSize = chr_size,stepSize = options.chunk_size):
				coord = interval[1]
				bw_signal = bw.values(chr_name,interval[1],interval[2])
				tmp = numpy.nansum(bw_signal)
				if numpy.isnan(tmp):continue
				bw_signal = numpy.nan_to_num(bw_signal) * weight
				for v in bw_signal:
					coord +=1
					if v != 0: print("%d\t%.2f" % (coord,v), file=OUT)
		elif options.out_format.upper() == "BGR":
			print("Writing " + chr_name + " ...", file=sys.stderr)
			#OUT.write('variableStep chrom='+chr_name+'\n')
			for interval in BED.tillingBed(chrName = chr_name,chrSize = chr_size,stepSize = options.chunk_size):
				v2p = collections.defaultdict(list)     #value to position
				range2p={}      #coorindate range to value, bedgraph. #[start]=[len,value]
				coord = interval[1]
				bw_signal = bw.values(chr_name,interval[1],interval[2])
				tmp = numpy.nansum(bw_signal)
				if numpy.isnan(tmp):continue
				bw_signal = numpy.nan_to_num(bw_signal) * weight
				for v in bw_signal:
					coord +=1
					if v != 0: v2p[v].append(coord)
				for v in v2p:
					for k,g in groupby(enumerate(v2p[v]), lambda i_x:i_x[0]-i_x[1]):
						for l in [list(map(itemgetter(1), g))]:
							range2p[l[0]-1] = [len(l),v]
				for i in sorted(range2p):
					print(chr_name + '\t' + str(i) +'\t' + str(i + range2p[i][0]) + '\t' + str(range2p[i][1]), file=OUT)
		else:
			print("unknown output format", file=sys.stderr)
			sys.exit(1)
				
			
			
if __name__=='__main__':
	main()
