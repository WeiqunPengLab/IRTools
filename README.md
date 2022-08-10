# IRTools
IRTools is a computational toolset for detection and quantitative analysis of intron retention from RNA-Seq libraries.

## Table of Contents
- [Installation](#installation)
- [Usage](#usage)
- [Tutorial](#tutorial)

## Installation

To use `IRTools`, you will need `python 3.6.9` or potentially above. 

Note: This has only been tested in `python 3.6.9`.

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

Note: This has only been tested with `numpy 1.19.5`, `scipy 1.5.4`, `pandas 1.1.5`, `networkx  2.5.1`, and `bx-python 0.8.12`, but will likely work with any later version.

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

`ANNOFILE` is the output GTF file that contains information for intron retention analysis, including the genomic coordinates of constitutive intronic regions (CIRs), constitutive exonic regions (CERs), constitutive exon-intron junctions (CJs), etc.

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

**-i/--alt-file SAMPLELIB**

Input RNA-Seq alignment file. The input file must be a BAM file.

**-p/--read-type {paired,single}**

"paired" is for paired-end data and "single" is for single-end data. DEFAULT: "single".

Note: See [Tutorial](#tutorial) for more information.

**-s/--library-type {fr-unstranded,fr-firststrand,fr-secondstrand}**

Library type. DEFAULT: "fr-unstranded" (unstranded). Use "fr-firststrand" or "fr-secondstrand" for strand-specific data.

Note: See [Tutorial](#tutorial) for more information.

**-u/--map-file MAPFILE** (optional)

Mappability score bigWig file (depends on species,
                        sequence length of RNA-Seq library, etc.). Or specify a
                        species (i.e. hg19 or mm9) for which a default
                        annotation file (default for 50 bps of single end RNA-
                        Seq library) can be downloaded and used. If specified,
                        mappability will be taken into account.
                        
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
                        
**-m/--min_overlap MINOVERLAP** (specified when -q IRC)

Set when IR quantifiation type is "IRC". Minimum
                        length of overlap between the reads and each of the
                        exons or introns involved in splicing. DEFAULT: 8.

#### `Outputs`

**-q IRI**

1\. `NAME.quant.IRI.genes.txt` is the quantification of intron retention index for all genes from RNA-Seq library SAMPLELIB.

The file format is as follows.

| gene_id |	gene_iv | gene_CIR_length | gene_CER_length | gene_CIR_read_count | gene_CER_read_count | gene_CIR_RPKM | gene_CER_RPKM | gene_IRI |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A1BG | chr19:58858171-58864865 | 4071 | 1766 | 80.13888888888889 | 93.86111111111111 | 0.7394798014780853 | 1.9965449549590883 | 0.37037974008115343 |
| A1CF | chr10:52559168-52645435 | 72712 | 9221 | 24.0 | 59.0 | 0.012399074027101793 | 0.2403577285944172 | 0.051585917788498296 |

2\. `NAME.quant.IRI.introns.txt` is the quantification of intron retention index for all CIRs from RNA-Seq library SAMPLELIB.

The file format is as follows.

| CIR_id |	CIR_iv |	CIR_length |	adjacent_CER_length |	CIR_read_count |	adjacent_CER_read_count |	CIR_RPKM |	adjacent_CER_RPKM |	intron_IRI |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A1BG:004 | chr19:58863053-58863648 | 474 | 570 | 7.138888888888889 | 17.02777777777778 | 0.5657653978142747 | 1.122192132703333 | 0.5041609019761704 |
| A1BG:005| chr19:58862017-58862756| 739 | 579 | 18.13888888888889 | 25.22222222222222 | 0.9220412349334726 | 1.636397795045001 | 0.5634578815281992 |

<br>

**-q IRC**

1\. `NAME.quant.IRC.genes.txt` is the quantification of intron retention coefficient for all genes from RNA-Seq library SAMPLELIB.

The file format is as follows.

| gene_id |	gene_iv |	gene_retained_reads |	gene_spliced_reads |	gene_IRC |
| --- | --- | --- | --- | --- |
| AAAS | chr12:53701239-53715412 | 2.5 | 64 | 0.0375939849624 |
| AAGAB | chr15:67493012-67547536 | 0.5 | 29 | 0.0169491525424 |

2\. `NAME.quant.IRC.introns.txt` is the quantification of intron retention coefficient for all CIRs from RNA-Seq library SAMPLELIB.

The file format is as follows.

| CIR_id | CIR_iv | CIR_5'retained_reads | CIR_3'retained_reads | CIR_spliced_reads | intron_IRC |
| --- | --- | --- | --- | --- | --- |
| AAAS:012 | chr12:53702133-53702218 | 1 | 0 | 4 | 0.111111111111 |
| AAGAB:007 | chr15:67496486-67500899 | 0 | 1 | 5 | 0.0909090909091 |

3\. `NAME.quant.IRC.junctions.txt` is the quantification of intron retention coefficient for all CJs from RNA-Seq library SAMPLELIB.

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

2\. `NAME.diff.IRI.introns.txt` is the detection of differential intron retention index from two RNA-Seq samples S1FILES and S2FILES for all CIRs.

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

2\. `NAME.diff.IRC.introns.txt` is the detection of differential intron retention coefficient from two RNA-Seq samples S1FILES and S2FILES for all CIRs.

The file format is as follows.

| CIR_id |	PValue |	FDR |	intron_IRC_S1 |	intron_IRC_S2 |	intron_IRC_difference |
| --- | --- | --- | --- | --- | --- |

3\. `NAME.diff.IRC.junctions.txt` is the detection of differential intron retention coefficient from two RNA-Seq samples S1FILES and S2FILES for all CJs.

The file format is as follows.

| CJ_id |	PValue |	FDR |	junction_IRC_S1 |	junction_IRC_S2 |	junction_IRC_difference |
| --- | --- | --- | --- | --- | --- |


## Tutorial

The following tutorial will demonstrate how to use IRTools. We will be working with mm9 data. Due to constraints with file sizes, we do not provide the data for users to run the following commands on their own.

### Create an annotation file

The first step in the IRTools pipeline is creating an annotation file to be used for further intron retention analysis. The following is an example of using `IRTools annotation` to create the annotation file for the mm9 genome. 

```
IRTools annotation -g gtf/mm9_2015.gtf -o gtf/mm9_annotation.gtf 
```

Pre-prepared annotation files for hg19 and mm9 can be found in the [data](IRTools/data) directory. These annotation files were generated with GTF files downloaded from [iGenome](https://support.illumina.com/sequencing/sequencing_software/igenome.html).

### Determine data parameters

Before proceeding with intron retention analysis, you must determine whether the data is pair-end or single-end and the data's strand specificity. If this information is not known, you can run `infer_experiment.py` to determine these parameters. The following is an example of running `infer_experiment.py`.

Note: To run `infer_experiment.py`, download [RSeQC 4.0.0](IRTools/utility/RSeQC-4.0.0), and install: `python setup.py install`

```
infer_experiment.py -r ../IRTools/utility/RSeQC-4.0.0/gene_model/mm9_UCSC_knownGene.bed -i bam/B_0h_R1.bam
```

When the run is finished, we get the following output.

```
This is PairEnd Data
Fraction of reads failed to determine: 0.0001
Fraction of reads explained by "1++,1--,2+-,2-+": 0.5003
Fraction of reads explained by "1+-,1-+,2++,2--": 0.4996
```

It is important to keep a record of this output as the information will be needed when running `IRTools quant`. 

The first line of the output gives whether the data is pair-end or single-end. ***In this case, the data is pair-end, meaning that the BAM file must be sorted by name to run `IRTools quant` properly.*** If the data is single-end it does not matter whether or not the BAM file is sorted by name.

The second line of the output gives the fraction of total mapped reads where the strand specificity could not be determined. This number should be relatively low otherwise it could indicate that the data's quality is not great. 

The third and fourth lines give the fraction of mapped reads that can be attributed to a certain strand. In this example, half can be explained by "1++,1--,2+-,2-+" and half can be explained by "1+-,1-+,2++,2--". This indicates that the data is unstranded. If the large majority of reads were to be explained by "1++,1--,2+-,2-+", this would indicate that the data is strand specific to the second strand. If the large majority of reads were to be explained by "1+-,1-+,2++,2--", this would indicate that the data is strand specific to the first strand. The interpretation of the output is that same for single-end data except that rather than the mapped reads being explained by "1++,1--,2+-,2-+" and "1+-,1-+,2++,2--" they are explained by "++,--" and "+-,-+", respectively.

The steps above should be repeated with all replicates.

### Quantify intron retention

`IRTools quant` allows users to detect and quantify intron retention events in RNA-Seq data. The intron retention events can either be quantified as an intron retention index (IRI) or intron retention coefficient (IRC). The IRI of a CIR is defined as the ratio of its read density to the read density of its adjacent CERs and the IRI of a gene is defined as the ratio of the overall read density of CIRs in that gene to the overall read density of CERs. The IRC of a CIR is defined as the fraction of junction reads that are CJs (average of 5' CJ reads and 3' CJ reads) and the IRC of a gene is defined as the fraction of junction reads altogether in this gene that are CJs. Additionally, IRC also quantifies CJs for which the IRC is defined as the fraction of junction reads that are CJs. The following is an example of running `IRTools quant` in IRI mode with the parameters obtained above. Make sure to double check that the BAM file is sorted by name if the data is pair-end.

```
IRTools quant -q IRI -i bam/B_0h_R1.bam -p paired -s fr-unstranded -e mm9 -f BAM -n B_0h_R1 --outdir quant
```

Shown below are the first few lines of the output CIR IRI file and output gene IRI file, respectively.

```
CIR_id  CIR_iv  CIR_length      adjacent_CER_length     CIR_read_count  adjacent_CER_read_count CIR_RPKM        adjacent_CER_RPKM       intron_IRI
0610005C13Rik:001       chr7:52830147-52830496  349     220     0.0     0.9519439748852935      0.0     0.1353347974887704      0.0
0610005C13Rik:002       chr7:52829892-52829977  85      280     0.0     2.155574337921597       0.0     0.24078295617790033     0.0
0610005C13Rik:003       chr7:52826562-52829782  3093    110     0.0     1.2036303630363037      0.0     0.3422324752934783      0.0
0610005C13Rik:004       chr7:52823749-52826355  2467    585     0.0     0.0     0.0     0.0     NA
```

```
gene_id gene_iv gene_CIR_length gene_CER_length gene_CIR_read_count     gene_CER_read_count     gene_CIR_RPKM   gene_CER_RPKM   gene_IRI
0610005C13Rik   chr7:52823164-52830546  5994    915     0.0     2.155574337921597       0.0     0.0736822160981553      0.0
0610007P14Rik   chr12:87156404-87165495 7690    1185    143.0433093553559       748.9566906446437       0.581784274398963       19.76783572161043       0.02943085336160243
0610009B22Rik   chr11:51498886-51502136 2453    795     33.161477559187205      101.8385224408128       0.42282153013403        4.006505375049167       0.10553374838011824
0610009L18Rik   chr11:120209991-120212504       1894    619     31.046335986936107      23.953664013063893      0.5126856329027132      1.2103251933843449      0.4235932918731762
```

### Analyze differential intron retention

The final step of the IRTools pipeline is analyzing differential intron retention across samples using `IRTools diff`. The following is an example of running `IRTools diff` on data with three replicates each.

```
IRTools diff -q IRI --indir quant -s1 B_0h_R1,B_0h_R2,B_0h_R3 -s2 B_24h_R1,B_24h_R2,B_24h_R3 -n B_rest_active
```

Shown below are the first few lines of the output differential CIR IRI file and differential gene IRI file, respectively.

```
CIR_id  PValue  FDR     intron_IRI_S1   intron_IRI_S2   intron_IRI_difference
0610007P14Rik:001       0.01748570124638199     0.04794954044905743     0.007060475733529612,0.014128453952795113,0.014477353807235999  0.0015103179228320555,0.002840643279225302,0.0026121421370695623        -0.009567726718144601
0610007P14Rik:002       0.06787966337015987     0.12263988572369974     0.009760132590207683,0.010902219030707407,0.019954850865529797  0.005107327432614331,0.006744590566349023,0.004010156747295317  -0.008251709246728738
0610007P14Rik:003       0.005272071342986928    0.02435463738145186     0.04637347332930471,0.06178728248456752,0.07360228639288818     0.01638970338560156,0.017775747116409467,0.016974940064366832   -0.043540883880127515
0610007P14Rik:004       0.18363340347714693     0.2615404726932045      0.06460308673989847,0.06047582787680872,0.1100761804966993      0.0379698930660664,0.06959717944498957,0.033493469152358374     -0.031364851149997375
```

```
gene_id PValue  FDR     gene_IRI_S1     gene_IRI_S2     gene_IRI_difference
0610005C13Rik   0.3739009663000589      0.43369087370573584     0.0,0.0,0.0     0.03589700330832069,0.0,0.0     0.011965667769440231
0610007P14Rik   0.015597116208720256    0.029785444823973002    0.02943085336160243,0.04145926912655988,0.05402988190915748     0.011682591550220334,0.01425293641802157,0.01240164576089424    -0.028860943556061218
0610009B22Rik   0.004408064327194515    0.013804288139961773    0.10553374838011824,0.1268034455814951,0.15060043224605094      0.04221736494775407,0.050984225384648134,0.05459385479080759    -0.07838072702815149
0610009L18Rik   0.09138325586455809     0.12748981426961586     0.4235932918731762,0.37842494303340185,0.4278781983026458       0.598773864510687,1.282172938644743,1.9609292502639917  0.8706598734033992
```

An helpful step to check differential intron retention in individual CIRs/genes/junctions is to visualize them using [IGV](https://software.broadinstitute.org/software/igv/) or another genome viewer of your choice. The following is what the Gar1 gene looks like when visualized on IGV with this data.

<br>
<img src="IGV_example.PNG" width="850">
<br>
