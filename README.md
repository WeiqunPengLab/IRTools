# IRTools
IRTools is a computational toolset for detection and analysis of intron retention from RNA-Seq libraries.

## Installation



#### PIP

To install directly from PyPI:

```
pip install IRTools
```

If this fails, please install all dependencies first:

```
pip install HTSeq
pip install pysam
pip install pandas
pip install networkx
pip install bx-python==0.7.3
```


#### From source

To install from source:

1\. Install the dependencies with your favorite tools (`pip`, `conda`, etc.).

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


>### annotation

#### `Arguments`
**-g/--GTF-file GTFFILE**

Input annotation [GTF](http://mblab.wustl.edu/GTF22.html) File. GTF file for a specific species can be downloaded from [iGenome](https://support.illumina.com/sequencing/sequencing_software/igenome.html).

**-o/--annotation-file ANNOFILE**

Output annotation GTF file.

**--outdir**

If specified, all output files will be written to that directory. DEFAULT: the current working directory.

#### `Outputs`

`ANNOFILE` is the output GTF file that contains information for intron retention analysis, including the genomic coordinates of introns, exon-intron junctions, etc.

Sample lines are as follows.

```
chr1	IR_annotation	constitutive_intronic_region	3411983	3660632    .  	-	  .    downstream_constitutive_junction_number "002"; constitutive_intronic_region_number "001"; upstream_constitutive_junction_number "001"; gene_id "Xkr4"
chr1	IR_annotation	constitutive_junction	        3660632	3660632	   . 	-	  .	   constitutive_junction_type "5'_splice_junction"; constitutive_junction_number "001"; downstream "constitutive_intronic_region_number 001"; gene_id "Xkr4"; upstream "constitutive_exonic_region_number 001"
```

>### quant

#### `Arguments`
> General arguments:

**-q/--quant-type {IRI,IRC}**

IR quantifiation types: intron retention index (IRI), intron retention coefficient (IRC). DEFAULT: "IRI".

**-i/--alt-file ALTFILE**

Input RNA-Seq alignment file. If IR quantifiation type is "IRI", the input file can be BAM or BED file. If IR quantification type is "IRC", the input file can only be BAM file.

**-p/--read-type {paired,single}**

"paired" is for paired-end data and "single" is for single-end data. DEFAULT: "single".

**-f/--library-type {fr-unstranded,fr-firststrand,fr-secondstrand}**

Library type. DEFAULT: "fr-unstranded" (unstranded). Use "fr-firststrand" or "fr-secondstrand" for strand-specific data.

**-u/--map-file MAPFILE**

Mappability score bigWig file (depends on species,
                        sequence length of RNA-Seq library, etc). Or specify a
                        species (i.e. hg19 or mm9) for which a default
                        annotation file (default for 50 bps of single end RNA-
                        Seq library) can be downloaded and used. If specified,
                        mappability will take into account.
                        
**-e/--species {hg19,mm9}**

Specify a species for which integrated IR annotation
                        GTF file can be used.                    
<br>Note: -e and -g are mutually exclusive
                        and one is required.
                        

**-g/--annotation-file ANNOFILE**

IR annotation GTF file user-built by "IRTools
                        annotation" command. -e and -g are mutually exclusive
                        and one is required.
                        
**-n/--name NAME**

Sample name, which will be used to generate output
                        file names. REQUIRED.

**--outdir**

If specified, all output files will be written to that directory. DEFAULT: the current working directory.

#### `Outputs`

**-q IRI**

`ANNOFILE` is the output GTF file that contains information for intron retention analysis, including the genomic coordinates of introns, exon-intron junctions, etc.

Sample lines are as follows.

```
chr1	IR_annotation	constitutive_intronic_region	3411983	3660632    .  	-	  .    downstream_constitutive_junction_number "002"; constitutive_intronic_region_number "001"; upstream_constitutive_junction_number "001"; gene_id "Xkr4"
chr1	IR_annotation	constitutive_junction	        3660632	3660632	   . 	-	  .	   constitutive_junction_type "5'_splice_junction"; constitutive_junction_number "001"; downstream "constitutive_intronic_region_number 001"; gene_id "Xkr4"; upstream "constitutive_exonic_region_number 001"
```