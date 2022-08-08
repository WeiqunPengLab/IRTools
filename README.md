# IRTools
IRTools is a computational toolset for detection and analysis of intron retention from RNA-Seq libraries.

## Table of Contents
- [Installation](#installation)
- [Usage](#usage)
- [Tutorial](#tutorial)

## Installation

To use `IRTools`, you will need `python 3.6.9`.

#### PIP

To install directly from PyPI:

```
pip install IRTools
```

If this fails, please install all dependencies first:

```
pip install numpy
pip install scipy
pip install pandas
pip install networkx
pip install bx-python
pip install HTSeq==0.13.5
pip install pysam==0.19.1
```


#### From source

To install from source:

1\. Install the dependencies listed above.

2\. Run:

```
python setup.py install
```


## Usage

```
IRTools [-h] [-v] {annotation,quant,diff} ...
```

There are three major functions available in IRTools serving as sub-commands.

| Command | Function |
| --- | --- |
| annotation | Generate annotation GTF file for intron retention analysis. |
| quant | Quantify intron retention in both gene and intron levels. |
| diff | Detection of differential intron retention from two samples with replicates in both gene and intron levels. |

<br>
<br>

### annotation

#### `Arguments`
**-g/--GTF-file GTFFILE**

Input annotation [GTF](http://mblab.wustl.edu/GTF22.html) File. GTF file for a specific species can be downloaded from [iGenome](https://support.illumina.com/sequencing/sequencing_software/igenome.html).

**-o/--annotation-file ANNOFILE**

Output annotation GTF file.

**--outdir OUTDIR**

If specified, all output files will be written to that directory. DEFAULT: the current working directory.

#### `Outputs`

`ANNOFILE` is the output GTF file that contains information for intron retention analysis, including the genomic coordinates of introns, exon-intron junctions, etc.

Sample lines are as follows.

```
chr1	IR_annotation	constitutive_intronic_region	3411983  3660632     .  	-	  .    downstream_constitutive_junction_number "002"; constitutive_intronic_region_number "001"; upstream_constitutive_junction_number "001"; gene_id "Xkr4"
chr1	IR_annotation	constitutive_junction	        3660632	 3660632	 . 	    -	  .	   constitutive_junction_type "5'_splice_junction"; constitutive_junction_number "001"; downstream "constitutive_intronic_region_number 001"; gene_id "Xkr4"; upstream "constitutive_exonic_region_number 001"
```

<br>
<br>

### quant

#### `Arguments`

**-q/--quant-type {IRI,IRC}**

IR quantification types: intron retention index (IRI), intron retention coefficient (IRC). DEFAULT: "IRI".

**-i/--alt-file ALTFILE**

Input RNA-Seq alignment file. If IR quantification type is "IRI", the input file can be BAM or BED file. If IR quantification type is "IRC", the input file can only be BAM file.

**-p/--read-type {paired,single}**

"paired" is for paired-end data and "single" is for single-end data. DEFAULT: "single".

**-s/--library-type {fr-unstranded,fr-firststrand,fr-secondstrand}**

Library type. DEFAULT: "fr-unstranded" (unstranded). Use "fr-firststrand" or "fr-secondstrand" for strand-specific data.

**-u/--map-file MAPFILE** (optional)

Mappability score bigWig file (depends on species,
                        sequence length of RNA-Seq library, etc.). Or specify a
                        species (i.e. hg19 or mm9) for which a default
                        annotation file (default for 50 bps of single end RNA-
                        Seq library) can be downloaded and used. If specified,
                        mappability will take into account.
                        
Note: To take into account mappability, download [RSeQC 4.0.0](IRTools/utility/RSeQC-4.0.0), and install: `python setup.py install`
                        
**-e/--species {hg19,mm9}** (exclusive with -g)

Specify a species for which integrated IR annotation
                        GTF file can be used.                    
<br>Note: -e and -g are mutually exclusive
                        and one is required.
                        

**-g/--annotation-file ANNOFILE** (exclusive with -e)

IR annotation GTF file user-built by "IRTools
                        annotation" command. -e and -g are mutually exclusive
                        and one is required.
                        
**-n/--name NAME**

Sample name, which will be used to generate output
                        file names. REQUIRED.

**--outdir OUTDIR**

If specified, all output files will be written to that directory. DEFAULT: the current working directory.

**-f/--format {BAM,BED}** (specified when -q IRI)

Set when IR quantifiation type is "IRI". Specify input
                        RNA-Seq alignment file format: "BAM", "BED". DEFAULT:
                        "BAM".
                        
**-m/--min_overlap MINOVERLAP** (specified when -q IRC)

Set when IR quantifiation type is "IRC". Minimum
                        length of overlap between the reads and each of the
                        exons or introns involved in splicing. DEFAULT: 8.

#### `Outputs`

**-q IRI**

1\. `NAME.quant.IRI.genes.txt` is the quantification of intron retention index for all genes from RNA-Seq library ALTFILE.

The file format is as follows.

| gene_id |	gene_iv | gene_CIR_length | gene_CER_length | gene_CIR_read_count | gene_CER_read_count | gene_CIR_RPKM | gene_CER_RPKM | gene_IRI |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A1BG | chr19:58858171-58864865 | 4071 | 1766 | 80.13888888888889 | 93.86111111111111 | 0.7394798014780853 | 1.9965449549590883 | 0.37037974008115343 |
| A1CF | chr10:52559168-52645435 | 72712 | 9221 | 24.0 | 59.0 | 0.012399074027101793 | 0.2403577285944172 | 0.051585917788498296 |

2\. `NAME.quant.IRI.introns.txt` is the quantification of intron retention index for all introns from RNA-Seq library ALTFILE.

The file format is as follows.

| CIR_id |	CIR_iv |	CIR_length |	adjacent_CER_length |	CIR_read_count |	adjacent_CER_read_count |	CIR_RPKM |	adjacent_CER_RPKM |	intron_IRI |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A1BG:004 | chr19:58863053-58863648 | 474 | 570 | 7.138888888888889 | 17.02777777777778 | 0.5657653978142747 | 1.122192132703333 | 0.5041609019761704 |
| A1BG:005| chr19:58862017-58862756| 739 | 579 | 18.13888888888889 | 25.22222222222222 | 0.9220412349334726 | 1.636397795045001 | 0.5634578815281992 |

<br>

**-q IRC**

1\. `NAME.quant.IRC.genes.txt` is the quantification of intron retention coefficient for all genes from RNA-Seq library ALTFILE.

The file format is as follows.

| gene_id |	gene_iv |	gene_retained_reads |	gene_spliced_reads |	gene_IRC |
| --- | --- | --- | --- | --- |
| AAAS | chr12:53701239-53715412 | 2.5 | 64 | 0.0375939849624 |
| AAGAB | chr15:67493012-67547536 | 0.5 | 29 | 0.0169491525424 |

2\. `NAME.quant.IRC.introns.txt` is the quantification of intron retention coefficient for all introns from RNA-Seq library ALTFILE.

The file format is as follows.

| CIR_id | CIR_iv | CIR_5'retained_reads | CIR_3'retained_reads | CIR_spliced_reads | intron_IRC |
| --- | --- | --- | --- | --- | --- |
| AAAS:012 | chr12:53702133-53702218 | 1 | 0 | 4 | 0.111111111111 |
| AAGAB:007 | chr15:67496486-67500899 | 0 | 1 | 5 | 0.0909090909091 |

3\. `NAME.quant.IRC.junctions.txt` is the quantification of intron retention coefficient for all exon-intron junctions from RNA-Seq library ALTFILE.

The file format is as follows.

| CJ_id | CJ_iv | CJ_type | CJ_retained_reads | CJ_spliced_reads | junction_IRC |
| --- | --- | --- | --- | --- | --- |
| A1BG:012 | chr19:58859005-58859006 | 3'_splice_junction | 3 | 0 | 1.0 |
| AAAS:001 | chr12:53715125-53715126 | 5'_splice_junction | 0 | 12 | 0.0 |

<br>
<br>

### diff

#### `Arguments`

Arguments that are same as `IRTools quant`: 

**-q/--quant-type {IRI,IRC}**

**-n/--name NAME**

**--outdir OUTDIR**

<br>

Additional arguments:

**--indir INDIR**

The directory containing the IR quantification result
                        files to be analyzed. Default: the current working
                        directory

**-s1/--s1-files S1FILES**

A comma-separated list of names for each replicate in
                        sample 1. IR quantification result files referenced by
                        each name must be followed by the
                        ".quant.{IRI,IRC}.{introns,genes,junctions}.txt"
                        extension.
 
**-s2/--s2-files S1FILES**

A comma-separated list of names for each replicate in
                        sample 2. IR quantification result files referenced by
                        each name must be followed by the
                        ".quant.{IRI,IRC}.{introns,genes,junctions}.txt"
                        extension.                       
 
**-t/--analysis-type {P,U}**

Type of analysis performed. "P" is for paired
                        replicates analysis and "U" is for unpaired replicates
                        analysis. DEFAULT: "U".
                        
#### `Outputs`

**-q IRI**

1\. `NAME.diff.IRI.genes.txt` is the detection of differential intron retention index from two RNA-Seq samples S1FILES and S2FILES for all genes.

The file format is as follows.

| gene_id |	PValue |	FDR |	gene_IRI_S1 |	gene_IRI_S2 |	ene_IRI_difference |
| --- | --- | --- | --- | --- | --- |
| Aaas | 0.189649931299 | 0.558970463129 | 0.008,0.01,0.007 | 0.012,0.008,0.019 | 0.005 |
| 9530051G07Rik | 0.00474936472952 | 0.0199711015431 | 0.0,0.0,0.0 | 0.183,0.366,0.0 0.183 |

2\. `NAME.diff.IRI.introns.txt` is the detection of differential intron retention index from two RNA-Seq samples S1FILES and S2FILES for all introns.

The file format is as follows.

| CIR_id |	PValue |	FDR |	intron_IRI_S1 |	intron_IRI_S2 |	intron_IRI_difference |
| --- | --- | --- | --- | --- | --- |
| Aacs:014 | 1.0 | 1.0 | 0.0,0.004,0.0 | 0.003,0.005,0.03 | 0.011 |
| Abca7:017 | 1.97072661512e-54 | 1.13601945194e-52 | 0.108,0.394,0.281 | 0.0,0.033,0.0 | -0.25 |


<br>

**-q IRC**

1\. `NAME.diff.IRC.genes.txt` is the detection of differential intron retention coefficient from two RNA-Seq samples S1FILES and S2FILES for all genes.

The file format is as follows.

| gene_id |	PValue |	FDR |	gene_IRC_S1 |	gene_IRC_S2 |	gene_IRC_difference |
| --- | --- | --- | --- | --- | --- |

2\. `NAME.diff.IRC.introns.txt` is the detection of differential intron retention coefficient from two RNA-Seq samples S1FILES and S2FILES for all introns.

The file format is as follows.

| CIR_id |	PValue |	FDR |	intron_IRC_S1 |	intron_IRC_S2 |	intron_IRC_difference |
| --- | --- | --- | --- | --- | --- |

3\. `NAME.diff.IRC.junctions.txt` is the detection of differential intron retention coefficient from two RNA-Seq samples S1FILES and S2FILES for all exon-intron junctions.

The file format is as follows.

| CJ_id |	PValue |	FDR |	junction_IRC_S1 |	junction_IRC_S2 |	junction_IRC_difference |
| --- | --- | --- | --- | --- | --- |


## Tutorial

The following tutorial will take you through how to use IRTools and allow you to check that the installation was successful. You will be working with chromosome 19 of mm9 activated/resting B cell data. All data files can be found in the [Tutorial](Tutorial) directory. Note that all command line code in this tutorial is written with the assumption that the Tutorial directory is the user's current directory.

### Create an annotation file

The first step in the IRTools pipeline is creating an annotation file to be used for further intron retention analysis. Using `IRTools annotation`, create an annotation file for the mm9 genome. The estimated runtime is 5 minutes.

```
IRTools annotation -g gtf/mm9_2015.gtf -o gtf/mm9_annotation.gtf 
```

Next, check that the output is correct by comparing the annotation file you just created to the provided annotation file for the mm9 genome.

```
diff gtf/mm9_annotation.gtf gtf/tutorial.gtf
```

If there is no output after running the above command, then you have successfully created an annotation file for the mm9 genome.

### Determine data parameters

Before proceeding with intron retention analysis, you must determine whether the data is pair-end or single-end and the data's strand specificity. If this information is not known, you can run `infer_experiment.py` to determine these parameters. 

Note: To run `infer_experiment.py`, download [RSeQC 4.0.0](IRTools/utility/RSeQC-4.0.0), and install: `python setup.py install`

```
infer_experiment.py -r ../IRTools/utility/RSeQC-4.0.0/gene_model/mm9_UCSC_knownGene.bed -i bam/B0h_chr19_R1.bam
```

After running `infer_experiment.py` on the resting B cell replicate 1 BAM file, you should get the following output.

```
This is PairEnd Data
Fraction of reads failed to determine: 0.0001
Fraction of reads explained by "1++,1--,2+-,2-+": 0.5003
Fraction of reads explained by "1+-,1-+,2++,2--": 0.4996
```

Make sure to keep a record of this output as the information will be needed when running `IRTools quant`. 

The first line of the output gives whether the data is pair-end or single-end. In this case, the data is pair-end, meaning that the BAM file must be sorted by name to run `IRTools quant` properly. For the purpose of this tutorial, all BAM files have already been sorted by name. If the data is single-end it does not matter whether or not the BAM file is sorted by name.

The second line of the output gives the fraction of total mapped reads that could not be 

Now, try repeating the steps above with the other replicates and see what you get.
