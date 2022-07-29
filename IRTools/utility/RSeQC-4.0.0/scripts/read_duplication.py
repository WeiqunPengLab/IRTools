#!/usr/bin/env python
'''
Calculte reads' duplication rate. 
Sequence-based: Reads with identical sequence are considered as "duplicate reads".
Mapping-based: Reads mapped to the exact same location are considered as "duplicate reads".
'''

#import built-in modules
import os,sys
if sys.version_info[0] != 3:
	print("\nYou are using python" + str(sys.version_info[0]) + '.' + str(sys.version_info[1]) + " This verion of RSeQC needs python3!\n", file=sys.stderr)
	sys.exit()	

import re
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


def main():
	usage="%prog [options]" + '\n' + __doc__ + "\n"
	parser = OptionParser(usage,version="%prog " + __version__)
	parser.add_option("-i","--input-file",action="store",type="string",dest="input_file",help="Alignment file in BAM or SAM format.")
	parser.add_option("-o","--out-prefix",action="store",type="string",dest="output_prefix",help="Prefix of output files(s).")
	parser.add_option("-u","--up-limit",action="store",type="int",dest="upper_limit",default=500,help="Upper limit of reads' occurrence. Only used for plotting, default=%default (times)")
	parser.add_option("-q","--mapq",action="store",type="int",dest="map_qual",default=30,help="Minimum mapping quality (phred scaled) for an alignment to be considered as \"uniquely mapped\". default=%default")
	(options,args)=parser.parse_args()

	if not (options.output_prefix and options.input_file):
		parser.print_help()
		sys.exit(0)
	if os.path.exists(options.input_file):
		obj = SAM.ParseBAM(options.input_file)
		obj.readDupRate(outfile=options.output_prefix,up_bound=options.upper_limit, q_cut = options.map_qual)
		try:
			subprocess.call("Rscript " + options.output_prefix +  ".DupRate_plot.r", shell=True)
		except:
			pass
	else:
		print('\n\n' + options.input_file + " does NOT exists" + '\n', file=sys.stderr)
		#parser.print_help()
		sys.exit(0)
	



if __name__ == '__main__':
	main()
