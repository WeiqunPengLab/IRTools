#!/usr/bin/env python
'''
Equally divide BAM file (m alignments) into n parts. Each part contains roughly m/n alignments
that are randomly sampled from total alignments. 
'''

#import built-in modules
import os,sys
if sys.version_info[0] != 3:
    print("\nYou are using python" + str(sys.version_info[0]) + '.' + str(sys.version_info[1]) + " This verion of RSeQC needs python3!\n", file=sys.stderr)
    sys.exit()	

import string
from optparse import OptionParser
import warnings
import string
import collections
from random import randrange
#import third-party modules
import pysam

__author__ = "Liguo Wang"
__copyright__ = "Copyleft"
__credits__ = []
__license__ = "GPL"
__version__="4.0.0"
__maintainer__ = "Liguo Wang"
__email__ = "wang.liguo@mayo.edu"
__status__ = "Production"

def main():
	usage="%prog [options]" + '\n' + __doc__ + "\n"
	parser = OptionParser(usage,version="%prog " + __version__)
	parser.add_option("-i","--input-file",action="store",type="string",dest="input_file",help="Alignment file in BAM format. BAM file should be sorted and indexed.")
	parser.add_option("-n","--subset-num",action="store",type="int",dest="subset_num",help="Number of small BAM files")
	parser.add_option("-o","--out-prefix",action="store",type="string",dest="output_prefix",help="Prefix of output BAM files. Output \"Prefix_num.bam\".")
	parser.add_option("-s","--skip-unmap",action="store_true",dest="skip_unmap", help="Skip unmapped reads.")
	(options,args)=parser.parse_args()
		
	if not (options.input_file and options.subset_num and options.output_prefix):
		parser.print_help()
		sys.exit(0)
	if not os.path.exists(options.input_file):
		print('\n\n' + options.input_file + " does NOT exists" + '\n', file=sys.stderr)
		sys.exit(0)		
	
	samfile = pysam.Samfile(options.input_file,'rb')
	
	sub_bam = {}
	count_bam={}
	for i in range(0,options.subset_num):
		sub_bam[i] = pysam.Samfile(options.output_prefix + '_' + str(i) +'.bam','wb',template=samfile)
		count_bam[i] = 0
		
	total_alignment = 0
	print("Dividing " + options.input_file + " ...", end=' ', file=sys.stderr)
	try:
		while(1):
			aligned_read = next(samfile)
			if aligned_read.is_unmapped and options.skip_unmap is True:
				continue
			total_alignment += 1
			tmp = randrange(0,options.subset_num)
			sub_bam[tmp].write(aligned_read)
			count_bam[tmp] += 1
				
	except StopIteration:
		print("Done", file=sys.stderr)

	for i in range(0,options.subset_num):
		print("%-55s%d" %  (options.output_prefix + '_' + str(i) +'.bam', count_bam[i]))
				
if __name__ == '__main__':
	main()
