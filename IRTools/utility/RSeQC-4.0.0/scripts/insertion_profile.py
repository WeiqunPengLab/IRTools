#!/usr/bin/env python
'''
Calculate the distributions of inserted nucleotides across reads
Note CIGAR strings within SAM/BAM file should have 'I' operation 
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


def main():
	usage="%prog [options]" + '\n' + __doc__ + "\n"
	parser = OptionParser(usage,version="%prog " + __version__)
	parser.add_option("-i","--input-file",action="store",type="string",dest="input_file",help="Alignment file in BAM or SAM format.")
	parser.add_option("-o","--out-prefix",action="store",type="string",dest="output_prefix",help="Prefix of output files(s).")
	parser.add_option("-q","--mapq",action="store",type="int",dest="map_qual",default=30,help="Minimum mapping quality (phred scaled) for an alignment to be considered as \"uniquely mapped\". default=%default")
	parser.add_option("-s","--sequencing",action="store",dest="layout",help="Sequencing layout. \"SE\"(single-end) or \"PE\"(pair-end). ")
	(options,args)=parser.parse_args()

	if not (options.input_file and options.output_prefix and options.layout):
		parser.print_help()
		sys.exit(0)
	for input_file in ([options.input_file]):
		if not os.path.exists(input_file):
			print('\n\n' + input_file + " does NOT exists" + '\n', file=sys.stderr)
			sys.exit(0)

	obj = SAM.ParseBAM(options.input_file)
	if options.layout == "SE":
		obj.insertion_profile(outfile=options.output_prefix, q_cut = options.map_qual,PE=False)
	elif options.layout == "PE":
		obj.insertion_profile(outfile=options.output_prefix, q_cut = options.map_qual,PE=True)
	else:
		print("unknow sequencing layout. Must be \"SE\" or \"PE\"", file=sys.stderr)
	try:
		subprocess.call("Rscript " + options.output_prefix + '.insertion_profile.r',shell=True)
	except:
		print("Cannot generate pdf file from " + options.output_prefix + '.insertion_profile.r', file=sys.stderr)
		pass

if __name__ == '__main__':
	main()
