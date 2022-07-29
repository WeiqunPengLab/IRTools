#!/usr/bin/env python
'''
Calculate the RNA-seq reads coverage over gene body.
This module uses bigwig file as input.
'''
import sys
from sys import version_info, stderr, exit
if sys.version_info[0] != 3:
	print("\nYou are using python" + str(sys.version_info[0]) + '.' + str(sys.version_info[1]) + " This verion of RSeQC needs python3!\n", file=sys.stderr)
	sys.exit()	

from os import path
from optparse import OptionParser
from collections import defaultdict
from subprocess import call
from numpy import nan_to_num
from pyBigWig import open as openBigWig
from qcmodule import mystat

__author__ = "Liguo Wang, Santiago Revale"
__copyright__ = "Copyleft"
__credits__ = []
__license__ = "GPL"
__version__="4.0.0"
__maintainer__ = "Liguo Wang"
__email__ = "wang.liguo@mayo.edu"
__status__ = "Testing"


def coverageGeneBody_bigwig(bigFile, refbed, outfile, gtype='png'):
    '''Calculate reads coverage over gene body, from 5'to 3'. each gene will be equally divided
    into 100 regsions. bigFile is bigwig format file'''

    if refbed is None:
        print ("You must specify a bed file representing gene model", file=sys.stderr)
        exit(0)

    bw = openBigWig(bigFile, 'r')
    # Get chromosomes present in the bigwig file 
    chroms = bw.chroms().keys()

    print ("calculating coverage over gene body ...", file=sys.stderr)
    coverage = defaultdict(int)
    flag=0
    gene_count = 0
    with open(refbed, 'r') as handle:
        # Loop through the genes in the BED file
        for line in handle:
            try:
                if line.startswith(('#','track','browser')):
                    continue

                # Parse fields from gene tabls
                fields = line.split()
                chrom     = fields[0]
                tx_start  = int( fields[1] )
                tx_end    = int( fields[2] )
                geneName  = fields[3]
                strand    = fields[5]

                # Skip chromosomes present in the bed file but not present in the bigwig file.
                # This could happen with PATCHES or Unplaced chromosomes.
                if chrom not in chroms:
                    continue

                exon_starts = list(map( int, fields[11].rstrip( ',\n' ).split( ',' ) ))
                exon_starts = list(map((lambda x: x + tx_start ), exon_starts))
                exon_ends = list(map( int, fields[10].rstrip( ',\n' ).split( ',' ) ))
                exon_ends = list(map((lambda x, y: x + y ), exon_starts, exon_ends))
            except:
                print ("[NOTE: input bed must be 12-column] skipped this line: " + line, end = '\n', file=sys.stderr)
                continue

            # Count gene if it was properly read
            gene_count += 1

            gene_all_base=[]
            mRNA_len =0
            flag=0
            for st, end in zip(exon_starts, exon_ends):
                gene_all_base.extend(range(st+1, end+1))    # 0-based coordinates on genome
                mRNA_len = len(gene_all_base)
                if mRNA_len < 100:
                    flag = 1
                    break
            if flag == 1:
                continue

            # Sort coordinates according to the strand
            if strand == '-':
                gene_all_base.sort(reverse=True)
            else:
                gene_all_base.sort(reverse=False)

            # Get 100 points from each gene's coordinates
            percentile_base = []
            percentile_base = mystat.percentile_list(gene_all_base)

            for i in range(0, len(percentile_base)):
                sig = bw.values(chrom,percentile_base[i]-1,percentile_base[i])
                coverage[i] += nan_to_num(sig[0])

            print (" \t%d genes finished\r" % gene_count, end = ' ', file=sys.stderr)

    # Close bigwig file
    bw.close()
    print ("\n",file=sys.stderr)

    x_coord=[]
    y_coord=[]
    with open(outfile + ".geneBodyCoverage.txt", 'w') as handle:
        handle.write("percentile\tcount\n")
        for i in coverage:
            x_coord.append(str(i))
            y_coord.append(str(coverage[i]))
            handle.write("%i\t%i\n" % (i, coverage[i]))

    with open(outfile + ".geneBodyCoverage_plot.r", 'w') as handle:
        handle.write("%s(\'%s\')\n" % (gtype, outfile + ".geneBodyCoverage." + gtype))
        handle.write("x=1:100\n")
        handle.write("y=c(%s)\n" % ','.join(y_coord))
        handle.write("plot(x, y/%s, xlab=\"percentile of gene body (5'->3')\", ylab='average wigsum', type='s')\n" % gene_count)
        handle.write("dev.off()\n")

def main():
    usage = "%prog [options]" + '\n' + __doc__ + "\n"
    parser = OptionParser(usage, version="%prog " + __version__)
    parser.add_option("-i", "--input-file", action="store", type="string", dest="input_file", help="Coverage signal file in bigwig format")
    parser.add_option("-r", "--refgene", action="store", type="string", dest="ref_gene_model", help="Reference gene model in bed format. [required]")
    parser.add_option("-o", "--out-prefix", action="store", type="string", dest="output_prefix", help="Prefix of output files(s). [required]")
    parser.add_option("-t", "--graph-type", action="store", type="string", dest="graph_type", default="pdf", help="Graphic file type in \"pdf\", \"jpeg\", \"bmp\", \"bmp\", \"tiff\" or \"png\".default=%default [optional]")
    (options, args) = parser.parse_args()

    gt = options.graph_type.lower()
    if gt not in ('pdf', 'png', 'bmp', 'jpeg', 'tiff'):
        print ("graphic file type must be 'pdf' or 'png'", end = '\n', file=sys.stderr)
        parser.print_help()
        exit(0)

    if not (options.output_prefix and options.input_file and options.ref_gene_model):
        parser.print_help()
        exit(0)

    if not path.exists(options.ref_gene_model):
        print ('\n\n' + options.ref_gene_model + " does NOT exists", end='\n', file=sys.stderr)
        #parser.print_help()
        exit(0)

    if path.exists(options.input_file):
        coverageGeneBody_bigwig(options.input_file, options.ref_gene_model, options.output_prefix, gtype=options.graph_type)
        try:
            call("Rscript " + options.output_prefix + '.geneBodyCoverage_plot.r', shell=True)
        except:
            print ( "Cannot generate plot from " + options.output_prefix + ".geneBodyCoverage_plot.r", end = '\n', file=sys.stderr)
            pass
    else:
        print ('\n\n' + options.input_file + " does NOT exists", end = '\n', file=sys.stderr)
        #parser.print_help()
        exit(0)

if __name__ == '__main__':
        main()   
