import bisect
import collections
import itertools
import logging
import os

import HTSeq
import numpy as np
import pandas as pd
import pysam

from IRTools._resources import annotation_input_path


class ONT_IRI_quant(object):
        def __init__(self, args):
                self.params = args.__dict__.copy()
                self.stranded = self.is_stranded(self.params["libtype"])

                self.gtffile = self.load_gtffile()
                self.harmonize_gtffile_chromosomes()

                self.gene_id2iv = self.get_gene_iv()
                self.valid_genes = set(self.gene_id2iv.keys())
                self.gene_CIR_database = self.summarize_gene_CIR()
                self.gene_CER_database = self.summarize_gene_CER()
                self.CIR_id2iv = self.get_CIR_iv()
                self.CIR_length = self.get_CIR_length()
                self.CER_length = self.get_CER_length()
                self.CIR_adjacent_CER = self.get_CIR_adjacent_CER()
                self.gene_region, self.CIR_region, self.coverage_region, self.CIR_counts, self.coverage_bases = self.init_regions_and_counts()
                self.splice_index, self.splice_index_starts = self.init_splice_index()
                self.total_read_count = 0
                self.assigned_read_count = 0
                self.total_aligned_bases = 0

        @staticmethod
        def is_stranded(libtype):
                if libtype == "fr-firststrand" or libtype == "fr-secondstrand":
                        return True
                elif libtype == "fr-unstranded":
                        return False

        def load_gtffile(self):
                with annotation_input_path(self.params["annofile"]) as annofile:
                        self.valid_annofile(str(annofile))
                        gtffile = list(HTSeq.GFF_Reader(str(annofile), end_included=True))
                return gtffile

        @staticmethod
        def valid_annofile(annofile):
                gtffile = HTSeq.GFF_Reader(annofile, end_included=True)
                feature_types = set([feat.type for feat in itertools.islice(gtffile, 1000)])
                for feat in ["gene_region", "constitutive_exonic_region", "constitutive_intronic_region", "constitutive_junction"]:
                        if feat not in feature_types:
                                raise Exception("Annotations for \"{}\" are missed in {}. Please generate valid annotation GTF/GFF3 file by \"IRTools annotation\" command.".format(feat, annofile))

        @staticmethod
        def clean_attr_value(value):
                if isinstance(value, str):
                        return value.strip('"')
                return value

        def attr_value(self, feature, key):
                return self.clean_attr_value(feature.attr[key])

        @staticmethod
        def iv_to_str(iv):
                return iv.chrom + ":" + str(iv.start) + "-" + str(iv.end)

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
                with pysam.AlignmentFile(self.params["altfile"], "rb") as bamfile:
                        return set(bamfile.references)

        def get_gene_iv(self):
                gene_id2iv = {}
                for feature in self.gtffile:
                        gene_id = self.attr_value(feature, "gene_id")
                        if feature.type == "gene_region" and int(self.attr_value(feature, "constitutive_exonic_region_length")) != 0 and int(self.attr_value(feature, "constitutive_intronic_region_length")) != 0:
                                gene_id2iv[gene_id] = self.iv_to_str(feature.iv)
                return gene_id2iv

        @staticmethod
        def CIR_has_both_upstream_and_downstream_CERs(feature):
                return "upstream_constitutive_junction_number" in feature.attr and "downstream_constitutive_junction_number" in feature.attr

        def summarize_gene_CIR(self):
                gene_CIR_database = collections.defaultdict(lambda: {})
                for feature in self.gtffile:
                        gene_id = self.attr_value(feature, "gene_id")
                        if feature.type == "constitutive_intronic_region" and gene_id in self.valid_genes and self.CIR_has_both_upstream_and_downstream_CERs(feature):
                                CIR_number = self.attr_value(feature, "constitutive_intronic_region_number")
                                gene_CIR_database[gene_id][CIR_number] = feature
                return gene_CIR_database

        def summarize_gene_CER(self):
                gene_CER_database = collections.defaultdict(lambda: {})
                for feature in self.gtffile:
                        gene_id = self.attr_value(feature, "gene_id")
                        if feature.type == "constitutive_exonic_region" and gene_id in self.valid_genes:
                                CER_number = self.attr_value(feature, "constitutive_exonic_region_number")
                                gene_CER_database[gene_id][CER_number] = feature
                return gene_CER_database

        def get_CIR_iv(self):
                CIR_id2iv = {}
                for gene_id, database in self.gene_CIR_database.items():
                        for CIR_number, feature in database.items():
                                CIR_id2iv[gene_id + ":" + CIR_number] = self.iv_to_str(feature.iv)
                return CIR_id2iv

        def get_CIR_length(self):
                CIR_length = collections.defaultdict(lambda: collections.Counter())
                for gene_id, database in self.gene_CIR_database.items():
                        for CIR_number, feature in database.items():
                                CIR_length[gene_id][CIR_number] = feature.iv.length
                return CIR_length

        def get_CER_length(self):
                CER_length = collections.defaultdict(lambda: collections.Counter())
                for gene_id, database in self.gene_CER_database.items():
                        for CER_number, feature in database.items():
                                CER_length[gene_id][CER_number] = feature.iv.length
                return CER_length

        @staticmethod
        def region_number(region_label, expected_type):
                prefix = expected_type + "_number "
                if not region_label.startswith(prefix):
                        return None
                return region_label[len(prefix):]

        def get_CIR_adjacent_CER(self):
                CIR_adjacent_CER = collections.defaultdict(lambda: collections.defaultdict(set))
                for feature in self.gtffile:
                        if feature.type != "constitutive_junction":
                                continue

                        gene_id = self.attr_value(feature, "gene_id")
                        if gene_id not in self.valid_genes:
                                continue

                        upstream = self.attr_value(feature, "upstream")
                        downstream = self.attr_value(feature, "downstream")
                        CIR_number = self.region_number(upstream, "constitutive_intronic_region")
                        CER_number = self.region_number(downstream, "constitutive_exonic_region")
                        if CIR_number is None or CER_number is None:
                                CIR_number = self.region_number(downstream, "constitutive_intronic_region")
                                CER_number = self.region_number(upstream, "constitutive_exonic_region")

                        if CIR_number in self.gene_CIR_database[gene_id] and CER_number in self.gene_CER_database[gene_id]:
                                CIR_adjacent_CER[gene_id][CIR_number].add(CER_number)
                return CIR_adjacent_CER

        def init_regions_and_counts(self):
                gene_region = HTSeq.GenomicArrayOfSets("auto", stranded=self.stranded)
                CIR_region = HTSeq.GenomicArrayOfSets("auto", stranded=self.stranded)
                coverage_region = HTSeq.GenomicArrayOfSets("auto", stranded=self.stranded)
                CIR_counts = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.Counter()))
                coverage_bases = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.Counter()))

                for feature in self.gtffile:
                        gene_id = self.attr_value(feature, "gene_id")
                        if feature.type == "gene_region":
                                gene_region[self.htseq_iv(feature.iv)] += gene_id

                for gene_id, database in self.gene_CIR_database.items():
                        for CIR_number, feature in database.items():
                                CIR_region[self.htseq_iv(feature.iv)] += (gene_id, CIR_number)
                                coverage_region[self.htseq_iv(feature.iv)] += (gene_id, "constitutive_intronic_region", CIR_number)
                                CIR_counts[gene_id][CIR_number]["retained_reads"] = 0
                                CIR_counts[gene_id][CIR_number]["spliced_reads"] = 0
                                coverage_bases[gene_id]["constitutive_intronic_region"][CIR_number] = 0

                for gene_id, database in self.gene_CER_database.items():
                        for CER_number, feature in database.items():
                                coverage_region[self.htseq_iv(feature.iv)] += (gene_id, "constitutive_exonic_region", CER_number)
                                coverage_bases[gene_id]["constitutive_exonic_region"][CER_number] = 0

                return gene_region, CIR_region, coverage_region, CIR_counts, coverage_bases

        def splice_index_key(self, chrom, strand):
                if self.stranded:
                        return chrom, strand
                return chrom

        def init_splice_index(self):
                splice_index = collections.defaultdict(list)
                splice_index_starts = {}

                for gene_id, database in self.gene_CIR_database.items():
                        for CIR_number, feature in database.items():
                                key = self.splice_index_key(feature.iv.chrom, feature.iv.strand)
                                splice_index[key].append((feature.iv.start, feature.iv.end, gene_id, CIR_number))

                for key in list(splice_index.keys()):
                        splice_index[key].sort(key=lambda item: (item[0], item[1], item[2], item[3]))
                        splice_index_starts[key] = [item[0] for item in splice_index[key]]

                return splice_index, splice_index_starts

        @staticmethod
        def unique_aligned(alt):
                return True if (alt.get_tag("NH") if alt.has_tag("NH") else 1) == 1 else False

        def usable_alignment(self, alt):
                if alt.is_unmapped:
                        return False
                if not self.params["include_secondary"] and alt.is_secondary:
                        return False
                if not self.params["include_supplementary"] and alt.is_supplementary:
                        return False
                if alt.mapping_quality < self.params["minmapq"]:
                        return False
                if self.params["require_unique"] and not self.unique_aligned(alt):
                        return False
                return True

        @staticmethod
        def reverse_strand(iv):
                if iv.strand == "+":
                        return HTSeq.GenomicInterval(iv.chrom, iv.start, iv.end, "-")
                elif iv.strand == "-":
                        return HTSeq.GenomicInterval(iv.chrom, iv.start, iv.end, "+")
                return iv

        def read_blocks_to_intervals(self, alt):
                strand = "-" if alt.is_reverse else "+"
                return [
                        HTSeq.GenomicInterval(alt.reference_name, start, end, strand)
                        for start, end in alt.get_blocks()
                        if end > start
                ]

        def get_alt_iv(self, alt):
                libtype = self.params["libtype"]
                read_intervals = self.read_blocks_to_intervals(alt)
                if libtype == "fr-secondstrand" or libtype == "fr-unstranded":
                        alt_iv_seq = read_intervals
                elif libtype == "fr-firststrand":
                        alt_iv_seq = [self.reverse_strand(iv) for iv in read_intervals]
                return alt_iv_seq

        def read_strand(self, alt):
                strand = "-" if alt.is_reverse else "+"
                if self.params["libtype"] == "fr-firststrand":
                        return "-" if strand == "+" else "+"
                return strand

        @staticmethod
        def skipped_intervals(alt):
                ref_pos = alt.reference_start
                skipped = []
                if alt.cigartuples is None:
                        return skipped

                for operation, length in alt.cigartuples:
                        if operation == 3:
                                skipped.append((ref_pos, ref_pos + length))
                                ref_pos += length
                        elif operation in (0, 2, 7, 8):
                                ref_pos += length

                return skipped

        def collect_read_gene(self, alt_iv_seq):
                gene_union = set()

                for alt_iv in alt_iv_seq:
                        for iv, gene_step_set in self.gene_region[self.htseq_iv(alt_iv)].steps():
                                if gene_step_set:
                                        gene_union.update(gene_step_set)

                if len(gene_union) != 1:
                        return None
                return next(iter(gene_union))

        def collect_retained_candidates(self, alt_iv_seq, gene_id):
                retained_overlap = collections.Counter()

                for alt_iv in alt_iv_seq:
                        for iv, step_set in self.CIR_region[self.htseq_iv(alt_iv)].steps():
                                if not step_set:
                                        continue
                                for candidate_gene_id, CIR_number in step_set:
                                        if candidate_gene_id == gene_id:
                                                retained_overlap[(candidate_gene_id, CIR_number)] += iv.length

                min_overlap = self.params["minoverlap"]
                return {
                        key
                        for key, overlap in retained_overlap.items()
                        if overlap >= min_overlap
                }

        def accumulate_coverage_bases(self, alt_iv_seq, gene_id):
                for alt_iv in alt_iv_seq:
                        for iv, step_set in self.coverage_region[self.htseq_iv(alt_iv)].steps():
                                for candidate_gene_id, feature_type, region_number in step_set:
                                        if candidate_gene_id == gene_id:
                                                self.coverage_bases[gene_id][feature_type][region_number] += iv.length

        def collect_spliced_candidates(self, alt, gene_id):
                spliced_candidates = set()
                tolerance = self.params["splicetolerance"]
                strand = self.read_strand(alt)
                key = self.splice_index_key(alt.reference_name, strand)
                introns = self.splice_index.get(key, [])
                starts = self.splice_index_starts.get(key, [])

                for gap_start, gap_end in self.skipped_intervals(alt):
                        left = bisect.bisect_left(starts, gap_start - tolerance)
                        right = bisect.bisect_right(starts, gap_start + tolerance)
                        for index in range(left, right):
                                intron_start, intron_end, candidate_gene_id, CIR_number = introns[index]
                                if candidate_gene_id == gene_id and abs(intron_start - gap_start) <= tolerance and abs(intron_end - gap_end) <= tolerance:
                                        spliced_candidates.add((candidate_gene_id, CIR_number))

                return spliced_candidates

        def iter_reads(self):
                with pysam.AlignmentFile(self.params["altfile"], "rb") as bamfile:
                        for alt in bamfile.fetch(until_eof=True):
                                yield alt

        def quant(self):
                logging.info("Accumulating aligned ONT bases in constitutive intronic and exonic regions")

                for alt in self.iter_reads():
                        if not self.usable_alignment(alt):
                                continue

                        alt_iv_seq = self.get_alt_iv(alt)
                        if not alt_iv_seq:
                                continue

                        self.total_read_count += 1
                        self.total_aligned_bases += sum(iv.length for iv in alt_iv_seq)

                        gene_id = self.collect_read_gene(alt_iv_seq)
                        if gene_id is None:
                                continue

                        self.assigned_read_count += 1
                        self.accumulate_coverage_bases(alt_iv_seq, gene_id)
                        spliced_candidates = self.collect_spliced_candidates(alt, gene_id)
                        retained_candidates = self.collect_retained_candidates(alt_iv_seq, gene_id)

                        for candidate_gene_id, CIR_number in spliced_candidates:
                                self.CIR_counts[candidate_gene_id][CIR_number]["spliced_reads"] += 1

                        for candidate_gene_id, CIR_number in retained_candidates - spliced_candidates:
                                self.CIR_counts[candidate_gene_id][CIR_number]["retained_reads"] += 1

        @staticmethod
        def divide_or_nan(numerator, denominator):
                if denominator == 0:
                        return np.nan
                return numerator * 1.0 / denominator

        @staticmethod
        def divide_preserve_inf(numerator, denominator):
                if np.isnan(numerator) or np.isnan(denominator):
                        return np.nan
                if denominator == 0:
                        return np.nan if numerator == 0 else np.inf
                return numerator * 1.0 / denominator

        def aligned_bases_to_read_count(self, aligned_bases):
                if self.total_read_count == 0 or self.total_aligned_bases == 0:
                        return 0.0
                mean_aligned_length = self.total_aligned_bases * 1.0 / self.total_read_count
                return aligned_bases * 1.0 / mean_aligned_length

        def calc_coverage_rpkm(self, aligned_bases, feature_length):
                if self.total_aligned_bases == 0 or feature_length == 0:
                        return np.nan
                return aligned_bases / (feature_length / 1000.0) / (self.total_aligned_bases / 1000000.0)

        def output_ONT_intron_level(self):
                logging.info("Calculating coverage-based ONT IRI in intron level")

                rows = []
                for gene_id in sorted(self.CIR_counts.keys()):
                        for CIR_number in sorted(self.CIR_counts[gene_id].keys()):
                                CIR_id = gene_id + ":" + CIR_number
                                CIR_length = self.CIR_length[gene_id][CIR_number]
                                if CIR_length == 0:
                                        continue

                                CIR_aligned_bases = self.coverage_bases[gene_id]["constitutive_intronic_region"][CIR_number]
                                adjacent_CER_numbers = self.CIR_adjacent_CER[gene_id][CIR_number]
                                adjacent_CER_length = sum(self.CER_length[gene_id][number] for number in adjacent_CER_numbers)
                                adjacent_CER_aligned_bases = sum(
                                        self.coverage_bases[gene_id]["constitutive_exonic_region"][number]
                                        for number in adjacent_CER_numbers
                                )
                                CIR_read_count = self.aligned_bases_to_read_count(CIR_aligned_bases)
                                adjacent_CER_read_count = self.aligned_bases_to_read_count(adjacent_CER_aligned_bases)
                                CIR_RPKM = self.calc_coverage_rpkm(CIR_aligned_bases, CIR_length)
                                adjacent_CER_RPKM = self.calc_coverage_rpkm(adjacent_CER_aligned_bases, adjacent_CER_length)
                                intron_IRI = self.divide_preserve_inf(CIR_RPKM, adjacent_CER_RPKM)

                                retained_reads = self.CIR_counts[gene_id][CIR_number]["retained_reads"]
                                spliced_reads = self.CIR_counts[gene_id][CIR_number]["spliced_reads"]
                                total_reads = retained_reads + spliced_reads
                                intron_IR_event_ratio = self.divide_or_nan(retained_reads, total_reads)

                                rows.append({
                                        "CIR_id": CIR_id,
                                        "CIR_iv": self.CIR_id2iv[CIR_id],
                                        "CIR_length": CIR_length,
                                        "adjacent_CER_length": adjacent_CER_length,
                                        "CIR_read_count": CIR_read_count,
                                        "adjacent_CER_read_count": adjacent_CER_read_count,
                                        "CIR_RPKM": CIR_RPKM,
                                        "adjacent_CER_RPKM": adjacent_CER_RPKM,
                                        "intron_IRI": intron_IRI,
                                        "CIR_aligned_bases": CIR_aligned_bases,
                                        "adjacent_CER_aligned_bases": adjacent_CER_aligned_bases,
                                        "retained_reads": retained_reads,
                                        "spliced_reads": spliced_reads,
                                        "total_IR_reads": total_reads,
                                        "intron_IR_event_ratio": intron_IR_event_ratio,
                                })

                self.ONT_intron_level_df = pd.DataFrame(
                        rows,
                        columns=[
                                "CIR_id", "CIR_iv", "CIR_length", "adjacent_CER_length",
                                "CIR_read_count", "adjacent_CER_read_count", "CIR_RPKM",
                                "adjacent_CER_RPKM", "intron_IRI", "CIR_aligned_bases",
                                "adjacent_CER_aligned_bases", "retained_reads", "spliced_reads",
                                "total_IR_reads", "intron_IR_event_ratio",
                        ],
                )

                outfile = self.params["name"] + ".ont.IRI.introns.txt"
                outfile_fullpath = os.path.join(self.params["outdir"], outfile)
                logging.info("Writing ONT intron level result to file: {}".format(outfile_fullpath))
                self.ONT_intron_level_df.to_csv(outfile_fullpath, index=None, sep="\t", na_rep="NA")

        def output_ONT_gene_level(self):
                logging.info("Calculating coverage-based ONT IRI in gene level")

                rows = []
                for gene_id in sorted(self.CIR_counts.keys()):
                        gene_CIR_length = sum(self.CIR_length[gene_id].values())
                        gene_CER_length = sum(self.CER_length[gene_id].values())
                        if gene_CIR_length == 0:
                                continue

                        gene_CIR_aligned_bases = sum(self.coverage_bases[gene_id]["constitutive_intronic_region"].values())
                        gene_CER_aligned_bases = sum(self.coverage_bases[gene_id]["constitutive_exonic_region"].values())
                        gene_CIR_read_count = self.aligned_bases_to_read_count(gene_CIR_aligned_bases)
                        gene_CER_read_count = self.aligned_bases_to_read_count(gene_CER_aligned_bases)
                        gene_CIR_RPKM = self.calc_coverage_rpkm(gene_CIR_aligned_bases, gene_CIR_length)
                        gene_CER_RPKM = self.calc_coverage_rpkm(gene_CER_aligned_bases, gene_CER_length)
                        gene_IRI = self.divide_preserve_inf(gene_CIR_RPKM, gene_CER_RPKM)

                        retained_reads = sum(self.CIR_counts[gene_id][CIR_number]["retained_reads"] for CIR_number in self.CIR_counts[gene_id])
                        spliced_reads = sum(self.CIR_counts[gene_id][CIR_number]["spliced_reads"] for CIR_number in self.CIR_counts[gene_id])
                        total_reads = retained_reads + spliced_reads
                        gene_IR_event_ratio = self.divide_or_nan(retained_reads, total_reads)

                        rows.append({
                                "gene_id": gene_id,
                                "gene_iv": self.gene_id2iv[gene_id],
                                "gene_CIR_length": gene_CIR_length,
                                "gene_CER_length": gene_CER_length,
                                "gene_CIR_read_count": gene_CIR_read_count,
                                "gene_CER_read_count": gene_CER_read_count,
                                "gene_CIR_RPKM": gene_CIR_RPKM,
                                "gene_CER_RPKM": gene_CER_RPKM,
                                "gene_IRI": gene_IRI,
                                "gene_CIR_aligned_bases": gene_CIR_aligned_bases,
                                "gene_CER_aligned_bases": gene_CER_aligned_bases,
                                "gene_retained_reads": retained_reads,
                                "gene_spliced_reads": spliced_reads,
                                "gene_total_IR_reads": total_reads,
                                "gene_IR_event_ratio": gene_IR_event_ratio,
                        })

                self.ONT_gene_level_df = pd.DataFrame(
                        rows,
                        columns=[
                                "gene_id", "gene_iv", "gene_CIR_length", "gene_CER_length",
                                "gene_CIR_read_count", "gene_CER_read_count", "gene_CIR_RPKM",
                                "gene_CER_RPKM", "gene_IRI", "gene_CIR_aligned_bases",
                                "gene_CER_aligned_bases", "gene_retained_reads", "gene_spliced_reads",
                                "gene_total_IR_reads", "gene_IR_event_ratio",
                        ],
                )

                outfile = self.params["name"] + ".ont.IRI.genes.txt"
                outfile_fullpath = os.path.join(self.params["outdir"], outfile)
                logging.info("Writing ONT gene level result to file: {}".format(outfile_fullpath))
                self.ONT_gene_level_df.to_csv(outfile_fullpath, index=None, sep="\t", na_rep="NA")

        def output_ONT_genome_wide(self):
                total_CIR_length = sum(sum(lengths.values()) for lengths in self.CIR_length.values())
                total_CER_length = sum(sum(lengths.values()) for lengths in self.CER_length.values())
                total_CIR_aligned_bases = sum(
                        sum(counts["constitutive_intronic_region"].values())
                        for counts in self.coverage_bases.values()
                )
                total_CER_aligned_bases = sum(
                        sum(counts["constitutive_exonic_region"].values())
                        for counts in self.coverage_bases.values()
                )
                total_CIR_read_count = self.aligned_bases_to_read_count(total_CIR_aligned_bases)
                total_CER_read_count = self.aligned_bases_to_read_count(total_CER_aligned_bases)
                cir_percentage = self.divide_or_nan(total_CIR_aligned_bases * 100.0, self.total_aligned_bases)
                cer_percentage = self.divide_or_nan(total_CER_aligned_bases * 100.0, self.total_aligned_bases)
                total_CIR_RPKM = self.calc_coverage_rpkm(total_CIR_aligned_bases, total_CIR_length)
                total_CER_RPKM = self.calc_coverage_rpkm(total_CER_aligned_bases, total_CER_length)
                genome_wide_IRI = self.divide_preserve_inf(total_CIR_RPKM, total_CER_RPKM)

                total_retained_reads = sum(
                        counts["retained_reads"]
                        for gene_counts in self.CIR_counts.values()
                        for counts in gene_counts.values()
                )
                total_spliced_reads = sum(
                        counts["spliced_reads"]
                        for gene_counts in self.CIR_counts.values()
                        for counts in gene_counts.values()
                )
                total_IR_reads = total_retained_reads + total_spliced_reads

                if total_IR_reads == 0:
                        retained_percentage = np.nan
                        spliced_percentage = np.nan
                        genome_wide_event_ratio = np.nan
                else:
                        retained_percentage = total_retained_reads * 100.0 / total_IR_reads
                        spliced_percentage = total_spliced_reads * 100.0 / total_IR_reads
                        genome_wide_event_ratio = total_retained_reads * 1.0 / total_IR_reads

                logging.info("Genome-wide ONT intron retention statistics of the library:")
                print(("\tTotal CER equivalent read count: %f (%.2f%% of aligned bases)" % (total_CER_read_count, cer_percentage)))
                print(("\tTotal CIR equivalent read count: %f (%.2f%% of aligned bases)" % (total_CIR_read_count, cir_percentage)))
                print(("\tTotal retained reads: %i (%.2f%%)" % (total_retained_reads, retained_percentage)))
                print(("\tTotal spliced reads: %i (%.2f%%)" % (total_spliced_reads, spliced_percentage)))
                print(("\tEligible alignments: %i" % self.total_read_count))
                print(("\tAlignments assigned to one gene: %i" % self.assigned_read_count))
                print(("\tTotal aligned bases: %i" % self.total_aligned_bases))
                print(("\tGenome-wide coverage-based ONT IRI: %f" % genome_wide_IRI))
                print(("\tGenome-wide ONT IR event ratio: %f" % genome_wide_event_ratio))

                self.ONT_genome_wide_stats = {
                        "sample_name": self.params["name"],
                        "Total CIR": total_CIR_read_count,
                        "Total CER": total_CER_read_count,
                        "CIR ratio": cir_percentage,
                        "CER ratio": cer_percentage,
                        "Eligible alignment count": self.total_read_count,
                        "Assigned alignment count": self.assigned_read_count,
                        "Total aligned bases": self.total_aligned_bases,
                        "Genome-wide ONT IRI": genome_wide_IRI,
                        "Total retained reads": total_retained_reads,
                        "Total spliced reads": total_spliced_reads,
                        "Retained ratio": retained_percentage,
                        "Spliced ratio": spliced_percentage,
                        "Total IR reads": total_IR_reads,
                        "Genome-wide ONT IR event ratio": genome_wide_event_ratio,
                }
                return self.ONT_genome_wide_stats
