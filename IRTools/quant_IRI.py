import os
import logging
import re
import urllib.request, urllib.parse, urllib.error
import collections
import itertools
import gzip
import hashlib
import networkx as nx
import numpy as np
import pandas as pd
import bisect
import HTSeq
import warnings
import pysam
import tempfile
from functools import reduce

from IRTools._resources import annotation_input_path

class IRI_quant(object):       
        MAPFILE_SPECIES_ALIASES = {
                "mm11": "mm39",
        }
        SPECIES_BIGWIG_URLS = {
                "hg19": "https://hgdownload.soe.ucsc.edu/gbdb/hg19/hoffmanMappability/k{read_length}.Umap.MultiTrackMappability.bw",
                "hg38": "https://hgdownload.soe.ucsc.edu/gbdb/hg38/hoffmanMappability/k{read_length}.Umap.MultiTrackMappability.bw",
                "mm9": "https://hgdownload.soe.ucsc.edu/gbdb/mm9/hoffmanMappability/k{read_length}.Umap.MultiTrackMappability.bw",
                "mm10": "https://hgdownload.soe.ucsc.edu/gbdb/mm10/hoffmanMappability/k{read_length}.Umap.MultiTrackMappability.bw",
                # Inference from the Hoffman Lab raw layout. If unavailable upstream, the download step raises a clear error.
                "mm39": "https://bismap.hoffmanlab.org/raw/mm39/k{read_length}.umap.bw",
        }

        def __init__(self, args):
                self.params = args.__dict__.copy()
                self.stranded = self.is_stranded(self.params['libtype'])
                self.use_mappability = False
                
                self.gtffile = self.load_gtffile()
                self.harmonize_gtffile_chromosomes()
                
                self.gene_id2iv = self.get_gene_iv()
                self.valid_genes = set(self.gene_id2iv.keys())
                
                self.CIR_id2iv = self.get_CIR_iv()

                self.CIR_effective_length = self.get_CIR_effective_length()
                self.CER_length = self.get_CER_length()
                
                self.genes, self.gene_region, self.counts = self.init_GenomicArrayOfSets_and_Counter_for_quant_IRI()
                                       
                # bin filter
                self.bin_filter = True
                if self.bin_filter == True:
                        self.num_bins = 10
                        self.bins, self.bin_counts = self.init_GenomicArrayOfSets_and_Counter_for_bin_filter()
                
                self.G = self.get_constitutive_junction_graph()
                
        @staticmethod
        def is_stranded(libtype):
                if libtype == "fr-firststrand" or libtype == "fr-secondstrand":
                        return True
                elif libtype == "fr-unstranded":
                        return False
                
        def load_gtffile(self):
                with annotation_input_path(self.params['annofile']) as annofile:
                        self.valid_annofile(str(annofile))
                        gtffile = list(HTSeq.GFF_Reader(str(annofile), end_included=True))
                return gtffile
                                          
        @staticmethod        
        def valid_annofile(annofile):
                gtffile = HTSeq.GFF_Reader(annofile, end_included=True)
                feature_types = set([feat.type for feat in itertools.islice(gtffile, 1000)])
                for feat in ['gene_region', 'constitutive_exonic_region', 'constitutive_intronic_region', 'constitutive_junction']:
                        if feat not in feature_types:
                                raise Exception("Annotations for \"{}\" are missed in {}. Please generate valid annotation GTF file by \"IRTools annotation\" command.".format(feat, annofile))
        
        @staticmethod
        def iv_to_str(iv):
                return iv.chrom + ':' + str(iv.start) + '-' + str(iv.end)                        

        @staticmethod
        def clean_attr_value(value):
                if isinstance(value, str):
                        return value.strip('"')
                return value

        def attr_value(self, feature, key):
                return self.clean_attr_value(feature.attr[key])

        def htseq_iv(self, iv):
                if self.stranded:
                        return iv
                return HTSeq.GenomicInterval(iv.chrom, iv.start, iv.end, ".")

        @staticmethod
        def canonical_chrom_name(chrom):
                return chrom[3:] if chrom.startswith("chr") else chrom

        def harmonize_gtffile_chromosomes(self):
                bam_chroms = self.get_bam_chromosomes()
                if not bam_chroms:
                        return

                gtf_chroms = {feature.iv.chrom for feature in self.gtffile}
                if gtf_chroms.intersection(bam_chroms):
                        return

                bam_chrom_by_canonical = {}
                for chrom in bam_chroms:
                        bam_chrom_by_canonical.setdefault(self.canonical_chrom_name(chrom), chrom)

                canonical_overlap = {
                        self.canonical_chrom_name(chrom)
                        for chrom in gtf_chroms
                        if self.canonical_chrom_name(chrom) in bam_chrom_by_canonical
                }
                if not canonical_overlap:
                        return

                remapped_features = 0
                for feature in self.gtffile:
                        bam_chrom = bam_chrom_by_canonical.get(self.canonical_chrom_name(feature.iv.chrom))
                        if bam_chrom and bam_chrom != feature.iv.chrom:
                                feature.iv.chrom = bam_chrom
                                remapped_features += 1

                logging.info(
                        "Harmonized annotation chromosome names with BAM references for %i features",
                        remapped_features,
                )

        def get_bam_chromosomes(self):
                with pysam.AlignmentFile(self.params['altfile'], "rb") as bamfile:
                        return set(bamfile.references)

        def get_gene_iv(self):
                gene_id2iv = {}
                for feature in self.gtffile:
                        gene_id = self.attr_value(feature, "gene_id")
                        if feature.type == "gene_region" and int(self.attr_value(feature, "constitutive_exonic_region_length")) != 0 and int(self.attr_value(feature, "constitutive_intronic_region_length")) != 0:
                                gene_id2iv[gene_id] = self.iv_to_str(feature.iv)
                return gene_id2iv
        
        def get_CIR_iv(self):
                CIR_id2iv = {}
                for feature in self.gtffile:
                        gene_id = self.attr_value(feature, "gene_id")
                        if feature.type == "constitutive_intronic_region" and gene_id in self.valid_genes:
                                CIR_number = self.attr_value(feature, "constitutive_intronic_region_number")
                                CIR_id = gene_id + ':' + CIR_number
                                CIR_id2iv[CIR_id] = self.iv_to_str(feature.iv)
                return CIR_id2iv
        
        @staticmethod
        def download_mappability_file_by_species_name(species):
                url = 'http://hgdownload.cse.ucsc.edu/goldenPath/{}/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign50mer.bigWig'.format(species)
                response = urllib.request.urlopen(url)
                chunk_size = 10 * 1024 * 1024
                fname = species + '_' + url.split('/')[-1]
                with open(fname, 'wb') as f:
                        while True:
                                chunk_data = response.read(chunk_size)
                                if not chunk_data:
                                        break
                                f.write(chunk_data)
                return fname

        @classmethod
        def normalize_map_species(cls, species_name):
                return cls.MAPFILE_SPECIES_ALIASES.get(species_name, species_name)

        @staticmethod
        def cache_dir():
                return os.path.join(os.path.dirname(__file__), "data", "mappability_cache")

        def annotation_cache_token(self):
                if self.params['species']:
                        return "builtin-{}".format(self.params['species'])

                annofile = os.path.abspath(self.params['annofile'])
                stat = os.stat(annofile)
                token_input = "{}|{}|{}".format(annofile, stat.st_size, stat.st_mtime_ns)
                return hashlib.md5(token_input.encode("utf-8")).hexdigest()[:12]

        def cached_blacklist_path(self, cache_stem):
                safe_stem = re.sub(r'[^A-Za-z0-9._-]+', "_", cache_stem)
                return os.path.join(self.cache_dir(), safe_stem + ".bed.gz")

        def report_cache_status(self, cache_path, exists):
                if exists:
                        logging.info("Using cached low-mappability BED: %s", cache_path)
                        print("Using cached low-mappability BED: {}".format(cache_path))
                else:
                        logging.info("Low-mappability BED cache is missing: %s", cache_path)
                        print("Low-mappability BED cache not found: {}".format(cache_path))
                        print("Generating it now. This may take some time on the first run.")

        def iter_constitutive_introns(self):
                for feature in self.gtffile:
                        gene_id = self.attr_value(feature, "gene_id")
                        if feature.type == "constitutive_intronic_region" and gene_id in self.valid_genes:
                                yield feature

        def build_blacklist_from_bigwig(self, bigwig_path, blacklist_path, map_score_cutoff):
                from bx.bbi.bigwig_file import BigWigFile

                os.makedirs(os.path.dirname(blacklist_path), exist_ok=True)
                with open(bigwig_path, "rb") as bigwig_handle:
                        bigwig = BigWigFile(bigwig_handle)
                        with gzip.open(blacklist_path, "wt", encoding="utf-8") as out_handle:
                                for feature in self.iter_constitutive_introns():
                                        chrom = feature.iv.chrom.encode("utf-8")
                                        entries = bigwig.get(chrom, feature.iv.start, feature.iv.end) or []
                                        for start, end, score in entries:
                                                if score < map_score_cutoff:
                                                        out_handle.write(
                                                                "{}\t{}\t{}\t.\t0\t{}\n".format(
                                                                        feature.iv.chrom,
                                                                        start,
                                                                        end,
                                                                        feature.iv.strand,
                                                                )
                                                        )

        def download_to_tempfile(self, url, suffix):
                temp_handle = tempfile.NamedTemporaryFile(delete=False, suffix=suffix)
                temp_path = temp_handle.name
                temp_handle.close()
                logging.info("Downloading mappability resource from %s", url)
                print("Downloading mappability resource from {}".format(url))
                try:
                        with urllib.request.urlopen(url) as response, open(temp_path, "wb") as out_handle:
                                while True:
                                        chunk = response.read(10 * 1024 * 1024)
                                        if not chunk:
                                                break
                                        out_handle.write(chunk)
                except Exception:
                        if os.path.exists(temp_path):
                                os.unlink(temp_path)
                        raise
                return temp_path

        def build_species_blacklist_cache(self, species_name, map_score_cutoff):
                normalized_species = self.normalize_map_species(species_name)
                if normalized_species not in self.SPECIES_BIGWIG_URLS:
                        raise Exception(
                                "\"{}\" is not a supported species shortcut for --map-file. Supported shortcuts are: {}.".format(
                                        species_name,
                                        ", ".join(sorted(list(self.SPECIES_BIGWIG_URLS.keys()) + list(self.MAPFILE_SPECIES_ALIASES.keys()))),
                                )
                        )

                read_length = self.params['mapreadlength']
                if read_length not in (24, 36, 50, 100):
                        raise Exception(
                                "Unsupported --map-read-length {}. Supported values are 24, 36, 50, and 100.".format(read_length)
                        )

                annotation_token = self.annotation_cache_token()
                cache_stem = "{}_k{}_cutoff{}_{}".format(normalized_species, read_length, str(map_score_cutoff).replace(".", "p"), annotation_token)
                blacklist_path = self.cached_blacklist_path(cache_stem)
                self.report_cache_status(blacklist_path, os.path.exists(blacklist_path))
                if os.path.exists(blacklist_path):
                        return blacklist_path

                source_url = self.SPECIES_BIGWIG_URLS[normalized_species].format(read_length=read_length)
                temp_bigwig = self.download_to_tempfile(source_url, ".bw")
                try:
                        self.build_blacklist_from_bigwig(temp_bigwig, blacklist_path, map_score_cutoff)
                except Exception as exc:
                        raise Exception(
                                "Failed to build the low-mappability BED cache from {}. If this genome/read-length combination is unavailable upstream, please provide a local bigWig or BED file via --map-file. Original error: {}".format(
                                        source_url,
                                        exc,
                                )
                        )
                finally:
                        if os.path.exists(temp_bigwig):
                                os.unlink(temp_bigwig)

                logging.info("Generated low-mappability BED cache: %s", blacklist_path)
                print("Generated low-mappability BED cache: {}".format(blacklist_path))
                return blacklist_path

        def build_local_bigwig_blacklist_cache(self, bigwig_path, map_score_cutoff):
                stat = os.stat(bigwig_path)
                source_token = hashlib.md5(
                        "{}|{}|{}".format(os.path.abspath(bigwig_path), stat.st_size, stat.st_mtime_ns).encode("utf-8")
                ).hexdigest()[:12]
                cache_stem = "{}_k{}_cutoff{}_{}".format(
                        os.path.splitext(os.path.basename(bigwig_path))[0],
                        self.params['mapreadlength'],
                        str(map_score_cutoff).replace(".", "p"),
                        "{}_{}".format(self.annotation_cache_token(), source_token),
                )
                blacklist_path = self.cached_blacklist_path(cache_stem)
                self.report_cache_status(blacklist_path, os.path.exists(blacklist_path))
                if os.path.exists(blacklist_path):
                        return blacklist_path

                self.build_blacklist_from_bigwig(bigwig_path, blacklist_path, map_score_cutoff)
                logging.info("Generated low-mappability BED cache: %s", blacklist_path)
                print("Generated low-mappability BED cache: {}".format(blacklist_path))
                return blacklist_path

        def apply_blacklist_bed(self, gene_map_score, bed_path):
                opener = gzip.open if bed_path.endswith(".gz") else open
                with opener(bed_path, "rt", encoding="utf-8") as bed_handle:
                        for line in bed_handle:
                                if not line.strip() or line.startswith("#"):
                                        continue
                                fields = line.rstrip("\n").split("\t")
                                chrom = fields[0]
                                start = int(fields[1])
                                end = int(fields[2])
                                strand = fields[5] if len(fields) >= 6 else "."
                                gene_map_score[self.htseq_iv(HTSeq.GenomicInterval(chrom, start, end, strand))] = 0

        def resolve_blacklist_bed(self, map_score_cutoff):
                mapfile = self.params['mapfile']
                if not mapfile:
                        return None

                normalized_species = self.normalize_map_species(mapfile)
                if normalized_species in self.SPECIES_BIGWIG_URLS:
                        return self.build_species_blacklist_cache(normalized_species, map_score_cutoff)

                if re.search(r'\.bed(\.gz)?$', mapfile, re.IGNORECASE):
                        return mapfile

                if re.search(r'(bigWig|\.bw)$', mapfile, re.IGNORECASE):
                        return self.build_local_bigwig_blacklist_cache(mapfile, map_score_cutoff)

                raise Exception(
                        "\"{}\" is neither a supported species shortcut, nor a local BED file, nor a local bigWig file.".format(mapfile)
                )
        
        def init_mappability_GenomicArray(self, map_score_cutoff):
                gene_map_score = HTSeq.GenomicArray("auto", stranded=self.stranded, typecode="i")

                for feature in self.gtffile:
                        if feature.type == "constitutive_exonic_region" or feature.type == "constitutive_intronic_region":
                                gene_map_score[self.htseq_iv(feature.iv)] = 1

                blacklist_bed = self.resolve_blacklist_bed(map_score_cutoff)
                if blacklist_bed:
                        self.apply_blacklist_bed(gene_map_score, blacklist_bed)
                return gene_map_score
        
        def get_CIR_effective_length(self):
                CIR_effective_length = collections.defaultdict( lambda: collections.Counter() )	
                for feature in self.gtffile:
                        gene_id = self.attr_value(feature, "gene_id")
                        if feature.type == "constitutive_intronic_region" and gene_id in self.valid_genes:
                                CIR_number = self.attr_value(feature, "constitutive_intronic_region_number")
                                if "constitutive_intronic_region_effective_length" in feature.attr:
                                        CIR_effective_length[gene_id][CIR_number] = int(self.attr_value(feature, "constitutive_intronic_region_effective_length"))
                                else:
                                        CIR_effective_length[gene_id][CIR_number] = feature.iv.length
                return CIR_effective_length
        
        def get_CER_length(self):
                CER_length = collections.defaultdict( lambda: collections.Counter() )
                for feature in self.gtffile:
                        gene_id = self.attr_value(feature, "gene_id")
                        if feature.type == "constitutive_exonic_region" and gene_id in self.valid_genes:
                                CER_number = self.attr_value(feature, "constitutive_exonic_region_number")
                                CER_length[gene_id][CER_number] = feature.iv.length
                return CER_length
        
        @staticmethod
        def get_effective_length(alt_iv_seq, gene_map_score=None):
                if gene_map_score is None:
                        if type(alt_iv_seq) == list:
                                return sum(alt_iv.length for alt_iv in alt_iv_seq)
                        return alt_iv_seq.length
                if type(alt_iv_seq) == list:
                        read_length = 0
                        for alt_iv in alt_iv_seq:
                                read_length += int(sum(gene_map_score[alt_iv]))
                        return read_length
                else:
                        return int(sum(gene_map_score[alt_iv_seq]))   

        def calc_rpkm(self, read_count, feature_length):
                if self.total_read_count == 0 or feature_length == 0:
                        return np.nan
                return read_count / (feature_length / 1000.0) / (self.total_read_count / 1000000.0)
                
        def init_GenomicArrayOfSets_and_Counter_for_quant_IRI(self):
                genes = HTSeq.GenomicArrayOfSets("auto", stranded=self.stranded)
                gene_region = HTSeq.GenomicArrayOfSets("auto", stranded=self.stranded)                
                counts = collections.defaultdict( lambda:  collections.defaultdict( lambda:  collections.Counter()))  
                
                for feature in self.gtffile:
                        gene_id = self.attr_value(feature, "gene_id")
                        if feature.type == "gene_region":
                                gene_region[self.htseq_iv(feature.iv)] += gene_id
                        elif feature.type == "constitutive_exonic_region" and gene_id in self.valid_genes:
                                CER_number = self.attr_value(feature, "constitutive_exonic_region_number")
                                genes[self.htseq_iv(feature.iv)] += (gene_id, feature.type, CER_number)
                                counts[gene_id]["constitutive_exonic_region"][CER_number] = 0
                        elif feature.type == "constitutive_intronic_region" and gene_id in self.valid_genes:
                                CIR_number = self.attr_value(feature, "constitutive_intronic_region_number")
                                if self.CIR_effective_length[gene_id][CIR_number] > 0: 
                                        genes[self.htseq_iv(feature.iv)] += (gene_id, feature.type, CIR_number) 
                                        counts[gene_id]["constitutive_intronic_region"][CIR_number] = 0
                return genes, gene_region, counts             
        
        @staticmethod
        def binnize(iv, effective_length, gene_map_score, num_bins):
                bin_iv_list = []
                iv_chrom, iv_strand = iv.chrom, iv.strand
                bin_length = effective_length // num_bins
                for i in range(num_bins):
                        if i == 0:
                                bin_start = iv.start
                        else:
                                bin_start = bin_end
                                
                        if i == num_bins - 1:
                                bin_iv_list.append(HTSeq.GenomicInterval(iv_chrom, bin_start, iv.end, iv_strand))
                        else:
                                bin_end = bin_start + bin_length
                                cum_bin_length = sum(gene_map_score[HTSeq.GenomicInterval(iv_chrom, bin_start, bin_end, iv_strand)])
                                
                                while cum_bin_length < bin_length:
                                        step_start = bin_end
                                        bin_end = step_start + bin_length - cum_bin_length
                                        cum_bin_length += sum(gene_map_score[HTSeq.GenomicInterval(iv_chrom, step_start, bin_end, iv_strand)])
                
                                bin_iv_list.append(HTSeq.GenomicInterval(iv_chrom, bin_start, bin_end, iv_strand))        
                return bin_iv_list      

        @staticmethod
        def binnize_without_mappability(iv, num_bins):
                bin_iv_list = []
                iv_chrom, iv_strand = iv.chrom, iv.strand
                iv_length = iv.end - iv.start
                for i in range(num_bins):
                        bin_start = iv.start + (iv_length * i) // num_bins
                        bin_end = iv.start + (iv_length * (i + 1)) // num_bins
                        bin_iv_list.append(HTSeq.GenomicInterval(iv_chrom, bin_start, bin_end, iv_strand))
                return bin_iv_list
        
        def init_GenomicArrayOfSets_and_Counter_for_bin_filter(self):
                bins = HTSeq.GenomicArrayOfSets("auto", stranded=self.stranded)                   
                bin_counts = collections.defaultdict( lambda:  collections.defaultdict( lambda:  collections.Counter()))
                        
                for feature in self.gtffile:
                        gene_id = self.attr_value(feature, "gene_id")
                        if feature.type == "constitutive_intronic_region" and gene_id in self.valid_genes:
                                CIR_number = self.attr_value(feature, "constitutive_intronic_region_number")
                                CIR_effective_length = self.CIR_effective_length[gene_id][CIR_number]
                                if CIR_effective_length >= self.num_bins:
                                        if self.use_mappability:
                                                bins_iv_list = self.binnize(feature.iv, CIR_effective_length, self.gene_map_score, self.num_bins)
                                        else:
                                                bins_iv_list = self.binnize_without_mappability(feature.iv, self.num_bins)
                                        for bin_number in range(self.num_bins):
                                                bins[self.htseq_iv(bins_iv_list[bin_number])] += (gene_id, CIR_number, bin_number) 
                                                bin_counts[gene_id][CIR_number][bin_number] = 0 
                                else:
                                        for bin_number in range(self.num_bins):
                                                bin_counts[gene_id][CIR_number][bin_number] = np.nan                                         
                return bins, bin_counts
        
        def get_constitutive_junction_graph(self):
                G = collections.defaultdict( lambda: nx.Graph() )
                for feature in self.gtffile:
                        gene_id = self.attr_value(feature, "gene_id")
                        if feature.type == "constitutive_junction":
                                upstream = self.attr_value(feature, 'upstream')
                                downstream = self.attr_value(feature, 'downstream')
                                cj_number = self.attr_value(feature, 'constitutive_junction_number')
                                G[gene_id].add_edge(upstream, downstream, CJ_number=cj_number)
                return G
        
        def init_Counter_for_quant(self):
                for gene_id in self.counts:
                        for region_type in self.counts[gene_id]:
                                for region_number in self.counts[gene_id][region_type]:
                                        self.counts[gene_id][region_type][region_number] = 0
                                        
                if self.bin_filter:
                        for gene_id in self.bin_counts:
                                for CIR_number in self.bin_counts[gene_id]:
                                        for bin_number in self.bin_counts[gene_id][CIR_number]:
                                                self.bin_counts[gene_id][CIR_number][bin_number] = 0           
        
        @staticmethod
        def unique_aligned(alt):
                return True if (alt.get_tag("NH") if alt.has_tag("NH") else 1) == 1 else False
        
        @staticmethod
        # iv's strand will be reversed.
        def reverse_strand(iv):
                if iv.strand == "+":
                        iv_reversed = HTSeq.GenomicInterval(iv.chrom, iv.start, iv.end, "-")
                elif iv.strand == "-":
                        iv_reversed = HTSeq.GenomicInterval(iv.chrom, iv.start, iv.end, "+")
                return iv_reversed        
        
        def read_blocks_to_intervals(self, alt):
                strand = "-" if alt.is_reverse else "+"
                return [
                        HTSeq.GenomicInterval(alt.reference_name, start, end, strand)
                        for start, end in alt.get_blocks()
                        if end > start
                ]

        def get_alt_iv(self, alt):
                libtype = self.params['libtype']
                read_intervals = self.read_blocks_to_intervals(alt)
                if libtype == "fr-secondstrand" or libtype == "fr-unstranded":
                        alt_iv_seq = read_intervals
                elif libtype == "fr-firststrand":
                        alt_iv_seq = [ self.reverse_strand(iv) for iv in read_intervals ]
                return alt_iv_seq
        
        def get_pair_alt_iv(self, alt1, alt2):
                libtype = self.params['libtype']
                alt1_intervals = self.read_blocks_to_intervals(alt1)
                alt2_intervals = self.read_blocks_to_intervals(alt2)
                if libtype == "fr-secondstrand" or libtype == "fr-unstranded":
                        alt_iv_seq1 = alt1_intervals
                        alt_iv_seq2 = [ self.reverse_strand(iv) for iv in alt2_intervals ]
                elif libtype == "fr-firststrand":
                        alt_iv_seq1 = [ self.reverse_strand(iv) for iv in alt1_intervals ]
                        alt_iv_seq2 = alt2_intervals
                return alt_iv_seq1, alt_iv_seq2      

        def iter_single_reads(self):
                with pysam.AlignmentFile(self.params['altfile'], "rb") as bamfile:
                        for alt in bamfile.fetch(until_eof=True):
                                yield alt

        def iter_paired_reads(self):
                pending_reads = {}
                with pysam.AlignmentFile(self.params['altfile'], "rb") as bamfile:
                        for alt in bamfile.fetch(until_eof=True):
                                if not alt.is_paired:
                                        continue
                                mate = pending_reads.pop(alt.query_name, None)
                                if mate is None:
                                        pending_reads[alt.query_name] = alt
                                        continue

                                if alt.is_read1:
                                        yield alt, mate
                                elif mate.is_read1:
                                        yield mate, alt
                                else:
                                        yield mate, alt
        
        def is_read_in_gene_region(self, alt_iv_seq):
                iset_intersection = None
                iset_union = None                                        
                for alt_iv in alt_iv_seq:			
                        for iv, gene_step_set in self.gene_region[self.htseq_iv(alt_iv)].steps():
                                if iset_intersection is None:
                                        iset_intersection = gene_step_set.copy()
                                else:
                                        iset_intersection.intersection_update(gene_step_set)					
                                if iset_union is None:
                                        iset_union = gene_step_set.copy()
                                else:
                                        iset_union.update(gene_step_set) 
                return len(iset_intersection) == 1 and len(iset_union) == 1 and iset_intersection == iset_union
        
        # In this case, step_set is (feature.attr["gene_id"], feature.type, feature.attr["constitutive_exonic_region_number"])
        # We're interested in only which gene_id the step_set belongs to. So we extract the gene_id (index = 0) from step_set as new_step_set.
        @staticmethod
        def set_sorted_in_feature(step_set, index):
                new_step_set = set()
                for item in step_set:
                        new_step_set.add(item[index])
                return new_step_set        
        
        def is_read_in_CIR_or_CER(self, alt_iv_seq):
                iset_union = None
                for alt_iv in alt_iv_seq:			
                        for iv, step_set in self.genes[self.htseq_iv(alt_iv)].steps():
                                gene_step_set = self.set_sorted_in_feature(step_set, 0)					
                                if iset_union is None:
                                        iset_union = gene_step_set.copy()
                                else:
                                        iset_union.update(gene_step_set)                             
                return len(iset_union) == 1

        def collect_read_assignments(self, alt_iv_seq):
                gene_intersection = None
                gene_union = None
                feature_gene_union = None
                region_hits = []
                bin_hits = []

                for alt_iv in alt_iv_seq:
                        htseq_alt_iv = self.htseq_iv(alt_iv)
                        for iv, gene_step_set in self.gene_region[htseq_alt_iv].steps():
                                if gene_intersection is None:
                                        gene_intersection = gene_step_set.copy()
                                else:
                                        gene_intersection.intersection_update(gene_step_set)
                                if gene_union is None:
                                        gene_union = gene_step_set.copy()
                                else:
                                        gene_union.update(gene_step_set)

                        for iv, step_set in self.genes[htseq_alt_iv].steps():
                                if not step_set:
                                        continue
                                gene_step_set = self.set_sorted_in_feature(step_set, 0)
                                if feature_gene_union is None:
                                        feature_gene_union = gene_step_set.copy()
                                else:
                                        feature_gene_union.update(gene_step_set)
                                region_hits.append((iv, next(iter(step_set))))

                        if self.bin_filter:
                                for iv, step_set in self.bins[htseq_alt_iv].steps():
                                        if step_set:
                                                bin_hits.append((iv, next(iter(step_set))))

                if gene_intersection is None or gene_union is None or feature_gene_union is None:
                        return None

                if len(gene_intersection) != 1 or len(gene_union) != 1 or gene_intersection != gene_union:
                        return None

                if len(feature_gene_union) != 1 or feature_gene_union != gene_union:
                        return None

                return self.get_effective_length(alt_iv_seq), region_hits, bin_hits
        
        def assign_read_to_region(self, alt_read_length, region_hits):
                if alt_read_length > 0:
                        for iv, feature_assignment in region_hits:
                                gene_id, feature_type, region_number = feature_assignment
                                overlap_length = iv.length if not self.use_mappability else self.get_effective_length(iv, self.gene_map_score)
                                self.counts[gene_id][feature_type][region_number] += overlap_length * 1.0 / alt_read_length 
                                        
        def assign_read_to_bin_filter(self, alt_read_length, bin_hits):
                if alt_read_length > 0:
                        for iv, bin_assignment in bin_hits:
                                gene_id, CIR_number, bin_number = bin_assignment
                                overlap_length = iv.length if not self.use_mappability else self.get_effective_length(iv, self.gene_map_score)
                                self.bin_counts[gene_id][CIR_number][bin_number] += overlap_length * 1.0 / alt_read_length  
        
        @staticmethod                                        
        def combine_pair_iv_seq(alt_first_iv_seq, alt_second_iv_seq):
                combine_alt_iv_seq = []
                alt_iv_seq = alt_first_iv_seq + alt_second_iv_seq
                alt_iv_seq.sort(key=lambda x: x.start)
                combine_alt_iv_seq.append(alt_iv_seq[0].copy())
                for alt in itertools.islice(alt_iv_seq, 1, None):
                        if alt.overlaps(combine_alt_iv_seq[-1]):
                                combine_alt_iv_seq[-1].extend_to_include(alt)
                        else:
                                combine_alt_iv_seq.append(alt.copy())
                return combine_alt_iv_seq              
        
        def quant(self):
                self.init_Counter_for_quant()
                self.total_read_count = 0
                
                logging.info("Counting number of reads that map to each individual constitutive intronic region (CIR) and constitutive exonic region (CER)")
                   
                # Input is bam file
                # Single end
                if self.params['readtype'] == "single":
                        for alt in self.iter_single_reads():
                                # Consider the alignments that are aligned and uniquely mapped.
                                if (not alt.is_unmapped) and self.unique_aligned(alt):
                                        self.total_read_count += 1
                                        alt_iv_seq = self.get_alt_iv(alt)

                                        # Eligible alignments are those mapped into one gene's either constitutive exonic region (CER) or constitutive intronic region (CIR).
                                        # For each eligible alignment, we count by fraction of length. i.e. If an alt has 50 bps, 30 bps in CER "001", 20 bps in CIR "001". Then, count in CER "001" is 0.6,
                                        # and count in CIR "001" is 0.4. (IRI is considered in intron level, so count is distributed in intron level)
                                        read_assignments = self.collect_read_assignments(alt_iv_seq)
                                        if read_assignments:
                                                alt_read_length, region_hits, bin_hits = read_assignments
                                                self.assign_read_to_region(alt_read_length, region_hits)
                                                if self.bin_filter:
                                                        self.assign_read_to_bin_filter(alt_read_length, bin_hits)

                elif self.params['readtype'] == "paired":
                        for alt_first, alt_second in self.iter_paired_reads():
                                if (not alt_first.is_unmapped) and self.unique_aligned(alt_first) and (not alt_second.is_unmapped) and self.unique_aligned(alt_second) and alt_first.reference_name == alt_second.reference_name:
                                        self.total_read_count += 1
                                        alt_first_iv_seq, alt_second_iv_seq = self.get_pair_alt_iv(alt_first, alt_second)
                                        alt_iv_seq = self.combine_pair_iv_seq(alt_first_iv_seq, alt_second_iv_seq)

                                        read_assignments = self.collect_read_assignments(alt_iv_seq)
                                        if read_assignments:
                                                alt_read_length, region_hits, bin_hits = read_assignments
                                                self.assign_read_to_region(alt_read_length, region_hits)
                                                if self.bin_filter:
                                                        self.assign_read_to_bin_filter(alt_read_length, bin_hits)

        @staticmethod
        def CIRs_in_consitutive_junction_graph(graph):
                return sorted([node for node in graph.nodes() if re.match('constitutive_intronic_region', node)])        
        
        @staticmethod
        def get_quantile_index(read_count_qantile_list, read_count):
                return bisect.bisect_right(read_count_qantile_list, read_count) - 1
        
        def empirical_bin_filter(self, row, read_count_qantile_list, bin_filter_cutoff_quantile_list):
                if row.CIR_read_count == 0 or np.isnan(row.bin_max_percentage):
                        return row.intron_IRI
                else:
                        if row.bin_max_percentage <= bin_filter_cutoff_quantile_list[self.get_quantile_index(read_count_qantile_list, row.CIR_read_count)]:
                                return row.intron_IRI
                        else:
                                five_ss_bin_column = "bin0_percentage"
                                three_ss_bin_column = "bin{}_percentage".format(self.num_bins-1)
                                if row.bin_max_percentage == row[five_ss_bin_column]:
                                        return "NA (5'AS)"
                                elif row.bin_max_percentage == row[three_ss_bin_column]:
                                        return "NA (3'AS)"
                                else:
                                        return "NA (unannotated exon)"        
                                                        
        def apply_bin_filter(self, df, outlier=0.01):
                logging.info("Appling filter to remove fake intron retention events")
                
                df = df.copy()
                df['gene_id'] = df.CIR_id.str[:-4]
                df['CIR_number'] = df.CIR_id.str[-3:]
                cir_keys = list(zip(df['gene_id'], df['CIR_number']))
                for i in range(self.num_bins):
                        bin_read_count_column = 'bin{}_read_count'.format(i)
                        df[bin_read_count_column] = [self.bin_counts[gene_id][cir_number][i] for gene_id, cir_number in cir_keys]
                        
                        bin_percentage_column = 'bin{}_percentage'.format(i)
                        df[bin_percentage_column] = df[bin_read_count_column] / df['CIR_read_count']
                
                bin_read_count_columns = ['bin{}_read_count'.format(i) for i in range(self.num_bins)]        
                bin_percentage_columns = ['bin{}_percentage'.format(i) for i in range(self.num_bins)]
                df['bin_max_percentage'] = df[bin_percentage_columns].max(axis=1)
                
                try: 
                        read_count_qantile_list = []
                        bin_filter_cutoff_quantile_list = []
                        bin_filter_df = df[pd.notnull(df.bin_max_percentage)]
                        positive_read_count_df = bin_filter_df[bin_filter_df.CIR_read_count > 0]
                        for q in [0, 25, 50, 75]:
                                lower = np.percentile(positive_read_count_df.CIR_read_count, q)
                                upper = np.percentile(positive_read_count_df.CIR_read_count, q+25)
                                bin_max_percentage_cutoff = np.percentile(bin_filter_df.loc[(bin_filter_df.CIR_read_count > lower) & (bin_filter_df.CIR_read_count <= upper), 'bin_max_percentage'], (1 - outlier) * 100)
                                
                                read_count_qantile_list.append(lower)
                                bin_filter_cutoff_quantile_list.append(bin_max_percentage_cutoff)

                        read_counts = df['CIR_read_count'].to_numpy()
                        bin_max = df['bin_max_percentage'].to_numpy()
                        intron_iri = df['intron_IRI'].astype(object).to_numpy(copy=True)
                        threshold_indexes = np.searchsorted(np.array(read_count_qantile_list), read_counts, side='right') - 1
                        threshold_indexes = np.clip(threshold_indexes, 0, len(bin_filter_cutoff_quantile_list) - 1)
                        cutoff_array = np.array(bin_filter_cutoff_quantile_list, dtype=float)[threshold_indexes]

                        eligible_mask = (read_counts > 0) & ~np.isnan(bin_max)
                        filtered_mask = eligible_mask & (bin_max > cutoff_array)
                        five_mask = filtered_mask & np.isclose(bin_max, df['bin0_percentage'].to_numpy(), equal_nan=False)
                        three_mask = filtered_mask & np.isclose(bin_max, df['bin{}_percentage'.format(self.num_bins - 1)].to_numpy(), equal_nan=False)
                        unannotated_mask = filtered_mask & ~five_mask & ~three_mask

                        intron_iri[five_mask] = "NA (5'AS)"
                        intron_iri[three_mask] = "NA (3'AS)"
                        intron_iri[unannotated_mask] = "NA (unannotated exon)"
                        df['intron_IRI'] = intron_iri
                        
                        self.filtered_CIR_id_list = list(df[df.intron_IRI.isin(["NA (5'AS)", "NA (3'AS)", "NA (unannotated exon)"])].CIR_id)
                except IndexError:
                        self.filtered_CIR_id_list = []

                logging.info("{} constitutive intronic regions (CIR) are unlikely to be intron retention events and are filtered".format(len(self.filtered_CIR_id_list)))
                
                return df.drop(labels=['gene_id', 'CIR_number', 'bin_max_percentage'] + bin_read_count_columns + bin_percentage_columns, axis=1)
                                                                              
        def output_IRI_intron_level(self):                 
                logging.info("Calculating CIR RPKM and CER RPKM in intron level")
                logging.info("Calculating intron retention index (IRI) in intron level")
                
                IRI_intron_level_data = []
                for gene_id in sorted(self.G.keys()):
                        consitutive_junction_graph = self.G[gene_id]
                        for CIR in self.CIRs_in_consitutive_junction_graph(consitutive_junction_graph):
                                CIR_number = CIR[-3:]
                                CIR_id = gene_id + ":" + CIR_number
                                CIR_effective_length = self.CIR_effective_length[gene_id][CIR_number]
                                if CIR_effective_length == 0: 
                                        continue
                                
                                CIR_read_count = self.counts[gene_id]["constitutive_intronic_region"][CIR_number]
                                CIR_RPKM = self.calc_rpkm(CIR_read_count, CIR_effective_length)
                                
                                adjacent_CER_list = consitutive_junction_graph.neighbors(CIR)
                                adjacent_CER_number_list = [item[-3:] for item in adjacent_CER_list]
                                adjacent_CER_read_count = 0
                                adjacent_CER_length = 0 
                                for adjacent_CER_number in adjacent_CER_number_list:
                                        adjacent_CER_read_count += self.counts[gene_id]["constitutive_exonic_region"][adjacent_CER_number]
                                        adjacent_CER_length += self.CER_length[gene_id][adjacent_CER_number]
                                adjacent_CER_RPKM = self.calc_rpkm(adjacent_CER_read_count, adjacent_CER_length)
                                # Note: numerical warning likely to be gererated here, but does not cause issues with results
                                # RuntimeWarning: invalid value encountered in true_divide
                                # RuntimeWarning: divide by zero encountered in true_divide
                                # "N/A" if there are no intron or exon reads (CIR_RPKM == adjacent_CER_RPKM == 0)
                                # "inf" if there are intron reads but no exon reads (CIR_RPKM==0 and adjacent_CER_RPKM > 0)
                                warnings.filterwarnings('ignore') # Ignores above warning
                                intron_IRI = np.divide(CIR_RPKM, adjacent_CER_RPKM)
                                warnings.resetwarnings()

                                IRI_intron_level_dict = {"CIR_id": CIR_id,
                                                         "CIR_iv": self.CIR_id2iv[CIR_id],
                                                         "CIR_length": CIR_effective_length,
                                                         "adjacent_CER_length": adjacent_CER_length,
                                                         "CIR_read_count": CIR_read_count,
                                                         "adjacent_CER_read_count": adjacent_CER_read_count,
                                                         "CIR_RPKM": CIR_RPKM,
                                                         "adjacent_CER_RPKM": adjacent_CER_RPKM,
                                                         "intron_IRI": intron_IRI}
                                
                                IRI_intron_level_data.append(IRI_intron_level_dict)
                                                                   
                self.IRI_intron_level_df = pd.DataFrame(IRI_intron_level_data, columns=["CIR_id", "CIR_iv", "CIR_length", "adjacent_CER_length", "CIR_read_count", "adjacent_CER_read_count", "CIR_RPKM", "adjacent_CER_RPKM", "intron_IRI"])                
                
                if self.bin_filter:
                        self.IRI_intron_level_df = self.apply_bin_filter(self.IRI_intron_level_df)
                
                outfile = self.params['name'] + ".quant.IRI.introns.txt" 
                outfile_fullpath = os.path.join(self.params['outdir'], outfile)
                logging.info("Writing intron level result to file: {}".format(outfile_fullpath))
                self.IRI_intron_level_df.to_csv(outfile_fullpath, index=None, sep='\t', na_rep="NA")
                
        def filter_CIR_id(self, filtered_CIR_id_list):
                for CIR_id in filtered_CIR_id_list:
                        gene_id, CIR_number = CIR_id[:-4], CIR_id[-3:]
                        try:
                                del self.CIR_effective_length[gene_id][CIR_number]
                        except KeyError:
                                pass
                        
                        try:
                                del self.counts[gene_id]["constitutive_intronic_region"][CIR_number]   
                        except KeyError:
                                pass
                
        def output_IRI_gene_level(self, filtered_CIR_id_list=None):                 
                logging.info("Counting number of reads that map to constitutive intronic regions (CIR) and constitutive exonic regions (CER) of each individual gene")
                logging.info("Calculating CIR RPKM and CER RPKM in gene level")
                logging.info("Calculating intron retention index (IRI) in gene level")   
                
                if filtered_CIR_id_list is None and self.bin_filter:
                        filtered_CIR_id_list = self.filtered_CIR_id_list
                self.filter_CIR_id(filtered_CIR_id_list)
        
                IRI_gene_level_data = []
                for gene_id in sorted(self.counts.keys()):   
                        gene_CIR_effective_length = sum(self.CIR_effective_length[gene_id].values())
                        gene_CER_length = sum(self.CER_length[gene_id].values())
                        
                        if gene_CIR_effective_length == 0:
                                continue
        
                        gene_CIR_read_count = sum(self.counts[gene_id]["constitutive_intronic_region"].values())
                        gene_CER_read_count = sum(self.counts[gene_id]["constitutive_exonic_region"].values())
                        
                        gene_CIR_RPKM = self.calc_rpkm(gene_CIR_read_count, gene_CIR_effective_length)
                        gene_CER_RPKM = self.calc_rpkm(gene_CER_read_count, gene_CER_length)
                        
                        # Note: numerical warning likely to be gererated here, but does not cause issues with results
                        # RuntimeWarning: invalid value encountered in true_divide
                        # RuntimeWarning: divide by zero encountered in true_divide
                        # "N/A" if there are no intron or exon reads (gene_CIR_RPKM == gene_CER_RPKM == 0)
                        # "inf" if there are intron reads but no exon reads (gene_CIR_RPKM==0 and gene_CER_RPKM > 0)
                        warnings.filterwarnings('ignore') # Ignores above warning
                        gene_IRI = np.divide(gene_CIR_RPKM, gene_CER_RPKM)
                        warnings.resetwarnings()

                        IRI_gene_level_dict = {"gene_id": gene_id,
                                               "gene_iv": self.gene_id2iv[gene_id],
                                               "gene_CIR_length": gene_CIR_effective_length,
                                               "gene_CER_length": gene_CER_length,
                                               "gene_CIR_read_count": gene_CIR_read_count,
                                               "gene_CER_read_count": gene_CER_read_count,
                                               "gene_CIR_RPKM": gene_CIR_RPKM,
                                               "gene_CER_RPKM": gene_CER_RPKM,
                                               "gene_IRI": gene_IRI}
                        
                        IRI_gene_level_data.append(IRI_gene_level_dict)
                        
                self.IRI_gene_level_df = pd.DataFrame(IRI_gene_level_data, columns=["gene_id", "gene_iv", "gene_CIR_length", "gene_CER_length", "gene_CIR_read_count", "gene_CER_read_count", "gene_CIR_RPKM", "gene_CER_RPKM", "gene_IRI"])   
                
                outfile = self.params['name'] + ".quant.IRI.genes.txt"
                outfile_fullpath = os.path.join(self.params['outdir'], outfile)
                logging.info("Writing gene level result to file: {}".format(outfile_fullpath))              
                self.IRI_gene_level_df.to_csv(outfile_fullpath, index=None, sep='\t', na_rep="NA")  
                
        def output_IRI_genome_wide(self):
                total_CIR_effective_length = reduce(lambda x,y: x + sum(y.values()), list(self.CIR_effective_length.values()), 0) 
                total_CER_length = reduce(lambda x,y: x + sum(y.values()), list(self.CER_length.values()), 0) 
                
                total_CIR_read_count = reduce(lambda x,y: x + sum(y["constitutive_intronic_region"].values()), list(self.counts.values()), 0) 
                total_CER_read_count = reduce(lambda x,y: x + sum(y["constitutive_exonic_region"].values()), list(self.counts.values()), 0)     

                if self.total_read_count == 0:
                        logging.warning("No eligible reads were counted for IRI quantification; reporting genome-wide statistics as NA.")
                        IRI_genome_wide = np.nan
                        cer_percentage = np.nan
                        cir_percentage = np.nan
                else:
                        cer_percentage = total_CER_read_count * 100.0 / self.total_read_count
                        cir_percentage = total_CIR_read_count * 100.0 / self.total_read_count
                        if total_CIR_effective_length == 0 or total_CER_length == 0 or total_CER_read_count == 0:
                                IRI_genome_wide = np.nan
                        else:
                                IRI_genome_wide = (total_CIR_read_count * 1.0 / total_CIR_effective_length) / (total_CER_read_count * 1.0 / total_CER_length)
                
                logging.info("Genome-wide intron retention statistics of the library:")
                print(("\tTotal CER read count: %f (%.2f%%)" % (total_CER_read_count, cer_percentage)))
                print(("\tTotal CIR read count: %f (%.2f%%)" % (total_CIR_read_count, cir_percentage)))
                print(("\tEligible read count: %i" % self.total_read_count))
                print(("\tGenome-wide intron retention index: %f" % IRI_genome_wide))

                self.IRI_genome_wide_stats = {
                        "sample_name": self.params['name'],
                        "Total CIR": total_CIR_read_count,
                        "Total CER": total_CER_read_count,
                        "CIR ratio": cir_percentage,
                        "CER ratio": cer_percentage,
                        "Eligible read count": self.total_read_count,
                        "Genome-wide IRI": IRI_genome_wide,
                }
                return self.IRI_genome_wide_stats
                
                
             
        
                        
                
                
                
                
