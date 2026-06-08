import os
import logging
import collections
import itertools

import numpy as np
import pandas as pd
import HTSeq
import pysam

from IRTools._resources import annotation_input_path


class RT_quant(object):
        def __init__(self, args):
                self.params = args.__dict__.copy()
                self.stranded = self.is_stranded(self.params['libtype'])

                self.gtffile = self.load_gtffile()
                self.harmonize_gtffile_chromosomes()

                self.gene_features, self.transcript_features, self.cer_features = self.index_annotation_features()
                self.valid_genes = self.get_valid_genes()
                self.gene_id2iv = self.get_gene_iv()
                self.cer_length = self.get_cer_length()
                self.rt_regions, self.rt_length = self.build_rt_regions()
                self.regions, self.counts = self.init_region_arrays_and_counts()

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
                for feat in ['gene_region', 'transcript_region', 'constitutive_exonic_region']:
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
                return iv.chrom + ':' + str(iv.start) + '-' + str(iv.end)

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

        def index_annotation_features(self):
                gene_features = {}
                transcript_features = collections.defaultdict(list)
                cer_features = collections.defaultdict(list)
                for feature in self.gtffile:
                        gene_id = self.attr_value(feature, "gene_id")
                        if feature.type == "gene_region":
                                gene_features[gene_id] = feature
                        elif feature.type == "transcript_region":
                                transcript_features[gene_id].append(feature)
                        elif feature.type == "constitutive_exonic_region":
                                cer_features[gene_id].append(feature)
                return gene_features, transcript_features, cer_features

        def get_valid_genes(self):
                valid_genes = set()
                for gene_id, feature in self.gene_features.items():
                        if int(self.attr_value(feature, "constitutive_exonic_region_length")) != 0:
                                valid_genes.add(gene_id)
                return valid_genes

        def get_gene_iv(self):
                gene_id2iv = {}
                for gene_id, feature in self.gene_features.items():
                        if gene_id in self.valid_genes:
                                gene_id2iv[gene_id] = self.iv_to_str(feature.iv)
                return gene_id2iv

        def get_cer_length(self):
                cer_length = collections.Counter()
                for gene_id, features in self.cer_features.items():
                        if gene_id not in self.valid_genes:
                                continue
                        cer_length[gene_id] = sum(feature.iv.length for feature in features)
                return cer_length

        def build_rt_regions(self):
                rt_regions = {}
                rt_length = collections.Counter()
                genes_by_chrom = collections.defaultdict(list)

                for gene_id, feature in self.gene_features.items():
                        if gene_id in self.valid_genes:
                                genes_by_chrom[feature.iv.chrom].append((feature.iv.start, feature.iv.end, gene_id, feature.iv))

                for chrom in genes_by_chrom:
                        genes_by_chrom[chrom].sort(key=lambda item: (item[0], item[1], item[2]))

                for chrom, gene_records in genes_by_chrom.items():
                        for idx, (_, _, gene_id, gene_iv) in enumerate(gene_records):
                                transcript_ivs = [feature.iv for feature in self.transcript_features.get(gene_id, [])]
                                if not transcript_ivs:
                                        continue

                                if gene_iv.strand == "+":
                                        furthest_tes = max(iv.end for iv in transcript_ivs)
                                        rt_end = furthest_tes + 10000
                                        for next_idx in range(idx + 1, len(gene_records)):
                                                next_gene_iv = gene_records[next_idx][3]
                                                if next_gene_iv.start >= furthest_tes:
                                                        rt_end = min(rt_end, next_gene_iv.start)
                                                        break
                                        if rt_end > furthest_tes:
                                                rt_iv = HTSeq.GenomicInterval(chrom, furthest_tes, rt_end, gene_iv.strand)
                                        else:
                                                rt_iv = None
                                else:
                                        furthest_tes = min(iv.start for iv in transcript_ivs)
                                        rt_start = max(0, furthest_tes - 10000)
                                        for prev_idx in range(idx - 1, -1, -1):
                                                prev_gene_iv = gene_records[prev_idx][3]
                                                if prev_gene_iv.end <= furthest_tes:
                                                        rt_start = max(rt_start, prev_gene_iv.end)
                                                        break
                                        if furthest_tes > rt_start:
                                                rt_iv = HTSeq.GenomicInterval(chrom, rt_start, furthest_tes, gene_iv.strand)
                                        else:
                                                rt_iv = None

                                if rt_iv is not None:
                                        rt_regions[gene_id] = rt_iv
                                        rt_length[gene_id] = rt_iv.length
                                else:
                                        rt_length[gene_id] = 0
                return rt_regions, rt_length

        def init_region_arrays_and_counts(self):
                regions = HTSeq.GenomicArrayOfSets("auto", stranded=self.stranded)
                counts = collections.defaultdict(lambda: collections.Counter())

                for gene_id in self.valid_genes:
                        counts[gene_id]["constitutive_exonic_region"] = 0
                        counts[gene_id]["readthrough_region"] = 0

                        for feature in self.cer_features.get(gene_id, []):
                                regions[self.htseq_iv(feature.iv)] += (gene_id, "constitutive_exonic_region")

                        if gene_id in self.rt_regions:
                                regions[self.htseq_iv(self.rt_regions[gene_id])] += (gene_id, "readthrough_region")

                return regions, counts

        @staticmethod
        def unique_aligned(alt):
                return True if (alt.get_tag("NH") if alt.has_tag("NH") else 1) == 1 else False

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
                libtype = self.params['libtype']
                read_intervals = self.read_blocks_to_intervals(alt)
                if libtype == "fr-secondstrand" or libtype == "fr-unstranded":
                        alt_iv_seq = read_intervals
                elif libtype == "fr-firststrand":
                        alt_iv_seq = [self.reverse_strand(iv) for iv in read_intervals]
                return alt_iv_seq

        def get_pair_alt_iv(self, alt1, alt2):
                libtype = self.params['libtype']
                alt1_intervals = self.read_blocks_to_intervals(alt1)
                alt2_intervals = self.read_blocks_to_intervals(alt2)
                if libtype == "fr-secondstrand" or libtype == "fr-unstranded":
                        alt_iv_seq1 = alt1_intervals
                        alt_iv_seq2 = [self.reverse_strand(iv) for iv in alt2_intervals]
                elif libtype == "fr-firststrand":
                        alt_iv_seq1 = [self.reverse_strand(iv) for iv in alt1_intervals]
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

        @staticmethod
        def get_effective_length(alt_iv_seq):
                return sum(alt_iv.length for alt_iv in alt_iv_seq)

        def calc_rpkm(self, read_count, feature_length):
                if self.total_read_count == 0 or feature_length == 0:
                        return np.nan
                return read_count / (feature_length / 1000.0) / (self.total_read_count / 1000000.0)

        @staticmethod
        def set_sorted_in_feature(step_set, index):
                new_step_set = set()
                for item in step_set:
                        new_step_set.add(item[index])
                return new_step_set

        def collect_read_assignments(self, alt_iv_seq):
                gene_union = None
                region_hits = []

                for alt_iv in alt_iv_seq:
                        htseq_alt_iv = self.htseq_iv(alt_iv)
                        for iv, step_set in self.regions[htseq_alt_iv].steps():
                                if not step_set:
                                        continue
                                gene_step_set = self.set_sorted_in_feature(step_set, 0)
                                if gene_union is None:
                                        gene_union = gene_step_set.copy()
                                else:
                                        gene_union.update(gene_step_set)
                                for assignment in step_set:
                                        region_hits.append((iv, assignment))

                if gene_union is None or len(gene_union) != 1:
                        return None

                return self.get_effective_length(alt_iv_seq), region_hits

        def assign_read_to_region(self, alt_read_length, region_hits):
                if alt_read_length <= 0:
                        return

                for iv, feature_assignment in region_hits:
                        gene_id, feature_type = feature_assignment
                        self.counts[gene_id][feature_type] += iv.length * 1.0 / alt_read_length

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
                logging.info("Counting exonic reads and read-through reads for RT quantification")
                self.total_read_count = 0

                if self.params['readtype'] == "single":
                        for alt in self.iter_single_reads():
                                if (not alt.is_unmapped) and self.unique_aligned(alt):
                                        self.total_read_count += 1
                                        alt_iv_seq = self.get_alt_iv(alt)
                                        read_assignments = self.collect_read_assignments(alt_iv_seq)
                                        if read_assignments:
                                                alt_read_length, region_hits = read_assignments
                                                self.assign_read_to_region(alt_read_length, region_hits)

                elif self.params['readtype'] == "paired":
                        for alt_first, alt_second in self.iter_paired_reads():
                                if (not alt_first.is_unmapped) and self.unique_aligned(alt_first) and (not alt_second.is_unmapped) and self.unique_aligned(alt_second) and alt_first.reference_name == alt_second.reference_name:
                                        self.total_read_count += 1
                                        alt_first_iv_seq, alt_second_iv_seq = self.get_pair_alt_iv(alt_first, alt_second)
                                        alt_iv_seq = self.combine_pair_iv_seq(alt_first_iv_seq, alt_second_iv_seq)
                                        read_assignments = self.collect_read_assignments(alt_iv_seq)
                                        if read_assignments:
                                                alt_read_length, region_hits = read_assignments
                                                self.assign_read_to_region(alt_read_length, region_hits)

        def output_RT_gene_level(self):
                logging.info("Calculating CER RPKM and RT RPKM in gene level")
                logging.info("Calculating read-through index (RT) in gene level")

                RT_gene_level_data = []
                for gene_id in sorted(self.valid_genes):
                        gene_feature = self.gene_features[gene_id]
                        cer_length = self.cer_length[gene_id]
                        rt_length = self.rt_length[gene_id]
                        cer_read_count = self.counts[gene_id]["constitutive_exonic_region"]
                        rt_read_count = self.counts[gene_id]["readthrough_region"]

                        if cer_length == 0:
                                continue

                        gene_CER_RPKM = self.calc_rpkm(cer_read_count, cer_length)
                        RT_RPKM = self.calc_rpkm(rt_read_count, rt_length)
                        with np.errstate(divide="ignore", invalid="ignore"):
                                rt_index = np.divide(RT_RPKM, gene_CER_RPKM)

                        RT_gene_level_data.append({
                                "gene_id": gene_id,
                                "gene_iv": self.gene_id2iv[gene_id],
                                "gene_strand": gene_feature.iv.strand,
                                "RT_iv": self.iv_to_str(self.rt_regions[gene_id]) if gene_id in self.rt_regions else "NA",
                                "gene_CER_length": cer_length,
                                "RT_length": rt_length,
                                "gene_CER_read_count": cer_read_count,
                                "RT_read_count": rt_read_count,
                                "gene_CER_RPKM": gene_CER_RPKM,
                                "RT_RPKM": RT_RPKM,
                                "gene_RT_index": rt_index,
                        })

                self.RT_gene_level_df = pd.DataFrame(
                        RT_gene_level_data,
                        columns=[
                                "gene_id",
                                "gene_iv",
                                "gene_strand",
                                "RT_iv",
                                "gene_CER_length",
                                "RT_length",
                                "gene_CER_read_count",
                                "RT_read_count",
                                "gene_CER_RPKM",
                                "RT_RPKM",
                                "gene_RT_index",
                        ],
                )

                outfile = self.params['name'] + ".quant.RT.genes.txt"
                outfile_fullpath = os.path.join(self.params['outdir'], outfile)
                logging.info("Writing gene level result to file: {}".format(outfile_fullpath))
                self.RT_gene_level_df.to_csv(outfile_fullpath, index=None, sep='\t', na_rep="NA")

        def output_RT_genome_wide(self):
                total_rt_length = sum(self.rt_length.values())
                total_cer_length = sum(self.cer_length.values())
                total_rt_read_count = sum(gene_counts["readthrough_region"] for gene_counts in self.counts.values())
                total_cer_read_count = sum(gene_counts["constitutive_exonic_region"] for gene_counts in self.counts.values())

                if self.total_read_count == 0:
                        logging.warning("No eligible reads were counted for RT quantification; reporting genome-wide statistics as NA.")
                        rt_ratio = np.nan
                        cer_ratio = np.nan
                        RT_genome_wide = np.nan
                else:
                        rt_ratio = total_rt_read_count * 100.0 / self.total_read_count
                        cer_ratio = total_cer_read_count * 100.0 / self.total_read_count
                        if total_rt_length == 0 or total_cer_length == 0 or total_cer_read_count == 0:
                                RT_genome_wide = np.nan
                        else:
                                total_RT_RPKM = self.calc_rpkm(total_rt_read_count, total_rt_length)
                                total_CER_RPKM = self.calc_rpkm(total_cer_read_count, total_cer_length)
                                RT_genome_wide = np.divide(total_RT_RPKM, total_CER_RPKM)

                logging.info("Genome-wide read-through statistics of the library:")
                print(("\tTotal CER read count: %f (%.2f%%)" % (total_cer_read_count, cer_ratio)))
                print(("\tTotal RT read count: %f (%.2f%%)" % (total_rt_read_count, rt_ratio)))
                print(("\tEligible read count: %i" % self.total_read_count))
                print(("\tGenome-wide read-through index: %f" % RT_genome_wide))

                self.RT_genome_wide_stats = {
                        "sample_name": self.params['name'],
                        "Total RT": total_rt_read_count,
                        "Total CER": total_cer_read_count,
                        "RT ratio": rt_ratio,
                        "CER ratio": cer_ratio,
                        "Eligible read count": self.total_read_count,
                        "Genome-wide RT": RT_genome_wide,
                }
                return self.RT_genome_wide_stats
