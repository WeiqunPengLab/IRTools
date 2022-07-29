import os
import logging
import re
import collections
import itertools
import bisect
import numpy as np
import pandas as pd
import HTSeq

class IRC_quant(object): 
        def __init__(self, args):
                self.params = args.__dict__.copy()
                self.stranded = self.is_stranded(self.params['libtype'])
                
                self.gtffile = self.load_gtffile()
                
                self.gene_id2iv = self.get_gene_iv()
                self.valid_genes = set(self.gene_id2iv.keys())
        
                self.CIR_id2iv = self.get_CIR_iv()                
                self.CJ_id2iv = self.get_CJ_iv()
                
                self.gene_CJ_database = self.summarize_gene_CJ()
                self.gene_CIR_database, self.gene_CIR_associated_CJ_database = self.summarize_gene_CIR()
                self.genes, self.gene_region, self.CER_region, self.gene_counts, self.CIR_counts, self.CJ_counts = self.init_GenomicArrayOfSets_and_Counter_for_quant_IRC()
                
                self.filter = True

        @staticmethod
        def is_stranded(libtype):
                if libtype == "fr-firststrand" or libtype == "fr-secondstrand":
                        return True
                elif libtype == "fr-unstranded":
                        return False
                
        def load_gtffile(self):
                if self.params['species']:
                        import pkg_resources
                        annofile = pkg_resources.resource_filename('IRTools', "data/" + self.params['species'] + "_IR_annotation.gtf.gz")
                        gtffile = HTSeq.GFF_Reader(annofile, end_included=True)
                elif self.params['annofile']:
                        self.valid_annofile(self.params['annofile'])
                        gtffile = HTSeq.GFF_Reader(self.params['annofile'], end_included=True)   
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

        def get_gene_iv(self):
                gene_id2iv = {}
                for feature in self.gtffile:
                        gene_id = feature.attr["gene_id"]
                        if feature.type == "gene_region" and int(feature.attr["constitutive_exonic_region_length"]) != 0 and int(feature.attr["constitutive_intronic_region_length"]) != 0:
                                gene_id2iv[gene_id] = self.iv_to_str(feature.iv)
                return gene_id2iv
        
        def get_CIR_iv(self):
                CIR_id2iv = {}
                for feature in self.gtffile:
                        gene_id = feature.attr["gene_id"]
                        if feature.type == "constitutive_intronic_region" and gene_id in self.valid_genes:
                                CIR_number = feature.attr["constitutive_intronic_region_number"]
                                CIR_id = gene_id + ':' + CIR_number
                                CIR_id2iv[CIR_id] = self.iv_to_str(feature.iv)
                return CIR_id2iv
        
        def get_CJ_iv(self):
                CJ_id2iv = {}
                for feature in self.gtffile:
                        gene_id = feature.attr["gene_id"]
                        if feature.type == "constitutive_junction" and gene_id in self.valid_genes:
                                CJ_number = feature.attr["constitutive_junction_number"]
                                CJ_id = feature.attr['gene_id'] + ':' + CJ_number
                                CJ_id2iv[CJ_id] = self.iv_to_str(feature.iv) 
                return CJ_id2iv
        
        @staticmethod                     
        def CIR_has_both_upstream_and_downstream_CERs(feature):
                return 'upstream_constitutive_junction_number' in list(feature.attr.keys()) and 'downstream_constitutive_junction_number' in list(feature.attr.keys())        
        
        def summarize_gene_CJ(self):
                gene_CJ_database = collections.defaultdict( lambda: {} )
                for feature in self.gtffile:
                        gene_id = feature.attr["gene_id"]
                        if feature.type == "constitutive_junction" and gene_id in self.valid_genes:
                                CJ_number = feature.attr["constitutive_junction_number"]
                                gene_CJ_database[gene_id][CJ_number] = feature
                return gene_CJ_database
        
        def summarize_gene_CIR(self):
                gene_CIR_database = collections.defaultdict( lambda: {} )
                gene_CIR_associated_CJ_database = collections.defaultdict( lambda: collections.defaultdict( lambda: {} ))
                for feature in self.gtffile:
                        gene_id = feature.attr["gene_id"]
                        if feature.type == "constitutive_intronic_region" and gene_id in self.valid_genes and self.CIR_has_both_upstream_and_downstream_CERs(feature):
                                CIR_number = feature.attr["constitutive_intronic_region_number"]
                                gene_CIR_database[gene_id][CIR_number] = feature
                                
                                upstream_CJ_number, downstream_CJ_number = feature.attr['upstream_constitutive_junction_number'], feature.attr['downstream_constitutive_junction_number']
                                gene_CIR_associated_CJ_database[gene_id][CIR_number]['upstream_constitutive_junction'] = self.gene_CJ_database[gene_id][upstream_CJ_number]
                                gene_CIR_associated_CJ_database[gene_id][CIR_number]['downstream_constitutive_junction'] = self.gene_CJ_database[gene_id][downstream_CJ_number] 
                return gene_CIR_database, gene_CIR_associated_CJ_database
                                      
        def init_GenomicArrayOfSets_and_Counter_for_quant_IRC(self):
                genes = HTSeq.GenomicArrayOfSets("auto", stranded=self.stranded)
                gene_region = HTSeq.GenomicArrayOfSets("auto", stranded=self.stranded)
                CER_region = HTSeq.GenomicArrayOfSets("auto", stranded=self.stranded)                
                       
                gene_counts = collections.defaultdict( lambda:  collections.Counter()) 
                CIR_counts = collections.defaultdict( lambda:  collections.defaultdict( lambda:  collections.Counter() ))      
                CJ_counts = collections.defaultdict( lambda:  collections.defaultdict( lambda:  collections.Counter() ))
                for feature in self.gtffile:
                        gene_id = feature.attr["gene_id"]
                        if feature.type == "gene_region":
                                gene_region[feature.iv] += gene_id
                        elif feature.type == "constitutive_exonic_region" and gene_id in self.valid_genes:
                                CER_number = feature.attr["constitutive_exonic_region_number"]
                                genes[feature.iv] += (gene_id, feature.type, CER_number)
                                CER_region[feature.iv] += "constitutive_exonic_region"
                        elif feature.type == "constitutive_intronic_region" and gene_id in self.valid_genes:
                                CIR_number = feature.attr["constitutive_intronic_region_number"]
                                genes[feature.iv] += (gene_id, feature.type, CIR_number)
                                if self.CIR_has_both_upstream_and_downstream_CERs(feature):
                                        CIR_counts[gene_id][CIR_number]["CIR_5'retained_reads"] = 0
                                        CIR_counts[gene_id][CIR_number]["CIR_3'retained_reads"] = 0
                                        CIR_counts[gene_id][CIR_number]["CIR_spliced_reads"] = 0                                        
                        elif feature.type == "constitutive_junction" and gene_id in self.valid_genes:
                                CJ_number = feature.attr["constitutive_junction_number"]
                                CJ_counts[gene_id][CJ_number]["CJ_retained_reads"] = 0
                                CJ_counts[gene_id][CJ_number]["CJ_spliced_reads"] = 0                                
                                
                return genes, gene_region, CER_region, gene_counts, CIR_counts, CJ_counts
        
        def init_Counter_for_quant(self):
                for read_type in self.gene_counts:
                        for gene_id in self.gene_counts[read_type]:
                                self.gene_counts[read_type][gene_id] = 0
                                
                for gene_id in self.CIR_counts:
                        for CIR_number in self.CIR_counts[gene_id]:
                                for read_type in self.CIR_counts[gene_id][CIR_number]:
                                        self.CIR_counts[gene_id][CIR_number][read_type] = 0
                                        
                for gene_id in self.CJ_counts:
                        for CJ_number in self.CJ_counts[gene_id]:
                                for read_type in self.CJ_counts[gene_id][CJ_number]:
                                        self.CJ_counts[gene_id][CJ_number][read_type] = 0
         
        @staticmethod
        def unique_aligned(alt):
                return True if dict(alt.optional_fields).get("NH", 1) == 1 else False
        
        @staticmethod
        # iv's strand will be reversed.
        def reverse_strand(iv):
                if iv.strand == "+":
                        iv_reversed = HTSeq.GenomicInterval(iv.chrom, iv.start, iv.end, "-")
                elif iv.strand == "-":
                        iv_reversed = HTSeq.GenomicInterval(iv.chrom, iv.start, iv.end, "+")
                return iv_reversed        
        
        def get_alt_iv(self, alt):
                libtype = self.params['libtype']
                if libtype == "fr-secondstrand" or libtype == "fr-unstranded":
                        alt_iv_seq = [ co.ref_iv for co in alt.cigar if co.type == "M" and co.size > 0 ]
                elif libtype == "fr-firststrand":
                        alt_iv_seq = [ self.reverse_strand(co.ref_iv) for co in alt.cigar if co.type == "M" and co.size > 0 ]   
                return alt_iv_seq
        
        def get_pair_alt_iv(self, alt1, alt2):
                libtype = self.params['libtype']
                if libtype == "fr-secondstrand" or libtype == "fr-unstranded":
                        alt_iv_seq1 = [ co.ref_iv for co in alt1.cigar if co.type == "M" and co.size > 0 ]
                        alt_iv_seq2 = [ self.reverse_strand(co.ref_iv) for co in alt2.cigar if co.type == "M" and co.size > 0 ]
                elif libtype == "fr-firststrand":
                        alt_iv_seq1 = [ self.reverse_strand(co.ref_iv) for co in alt1.cigar if co.type == "M" and co.size > 0 ]  
                        alt_iv_seq2 = [ co.ref_iv for co in alt2.cigar if co.type == "M" and co.size > 0 ] 
                return alt_iv_seq1, alt_iv_seq2      
        
        def is_read_in_gene_region(self, alt_iv_seq):
                iset_intersection = None
                iset_union = None                                        
                for alt_iv in alt_iv_seq:			
                        for iv, gene_step_set in self.gene_region[alt_iv].steps():
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
                        for iv, step_set in self.genes[alt_iv].steps():
                                gene_step_set = self.set_sorted_in_feature(step_set, 0)					
                                if iset_union is None:
                                        iset_union = gene_step_set.copy()
                                else:
                                        iset_union.update(gene_step_set)  

                return len(iset_union) == 1 
        
        def read_associated_gene(self, alt_iv_seq):
                for alt_iv in alt_iv_seq:			
                        for iv, gene_step_set in self.gene_region[alt_iv].steps():
                                return list(gene_step_set)[0]
                        
        @staticmethod
        def get_junction_pos_list(junction, overlap):
                gene_strand = junction.iv.strand
                junction_to_pos = junction.iv.start
                if gene_strand == "+":
                        junction_from_pos = junction_to_pos - 1
                        junction_to_pos_list = list(range(junction_to_pos, junction_to_pos + overlap))
                        junction_from_pos_list = list(range(junction_from_pos - overlap + 1, junction_from_pos + 1))
                elif gene_strand == "-":
                        junction_from_pos = junction_to_pos + 1
                        junction_to_pos_list = list(range(junction_to_pos - overlap + 1, junction_to_pos + 1))
                        junction_from_pos_list = list(range(junction_from_pos, junction_from_pos + overlap))
                
                return junction_from_pos_list, junction_to_pos_list
        
        @staticmethod                
        def generate_bisect_list(alt_iv_seq):
                start_list = []
                end_list = []
                last_alt_end = -2
                for alt in alt_iv_seq:
                        if alt.start - last_alt_end > 1:
                                start_list.append(alt.start)
                                end_list.append(alt.end - 1)
                                last_alt_end = end_list[-1]
                        else:
                                end_list[-1] = alt.end - 1
                                last_alt_end = end_list[-1]
                return (start_list, end_list)
        
        # input argument pos_list is a list of continuos positions that can be either (junction_from_pos - minimum_overlap, junction_from_pos] or 
        # [junction_to_pos, junction_to_pos - minimum_overlap) (junction joins one constitutive exonic region and one constitutive intronic region).
        # junction_from_pos is the last position in upstream region and junction_to_pos is the first position in downstream region.
        # (start_list, end_list) is extracted from an alignment. By using binary search, we can know if the pos_list is in the alignment. -1 and -2 mean not found.
        @staticmethod
        def find_pos_in_bisect_list(pos_list, start_list, end_list):
                start_index_set = set()
                end_index_set = set()
                for pos in pos_list:
                        start_index = bisect.bisect_right(start_list, pos)
                        end_index = bisect.bisect_left(end_list, pos)
                        start_index_set.add(start_index)
                        end_index_set.add(end_index)
                if len(start_index_set) == len(end_index_set) == 1:
                        if list(start_index_set)[0] - list(end_index_set)[0] == 1:
                                return list(end_index_set)[0]
                        else:
                                return -1
                else:
                        return -2   
                
        def is_spliced_read_entirely_in_CER(self, alt_pos_list):
                flag = True
                for pos in alt_pos_list:
                        if self.CER_region[pos]:
                                flag &= True
                        else:
                                flag &= False 
                                break     
                return flag
                                
        def assign_read_to_CJ(self, alt_iv_seq, CJ):     
                gene_id, gene_chrom, gene_strand = CJ.attr["gene_id"], CJ.iv.chrom, CJ.iv.strand
                CJ_number, CJ_type = CJ.attr["constitutive_junction_number"], CJ.attr["constitutive_junction_type"]
                overlap = self.params['minoverlap']
                
                start_list, end_list = self.generate_bisect_list(alt_iv_seq)
                
                junction_from_pos_list, junction_to_pos_list = self.get_junction_pos_list(CJ, overlap)
                junction_from_pos_index = self.find_pos_in_bisect_list(junction_from_pos_list, start_list, end_list)
                junction_to_pos_index = self.find_pos_in_bisect_list(junction_to_pos_list, start_list, end_list)              
                
                # -2 means for this regionicular junction, this read is neither across junction read nor splicing junction read considering the minimum overlap requirement.
                if junction_from_pos_index == -2 or junction_to_pos_index == -2:                   
                        return
                else:
                        # Junction is not spliced. in other words, intron retention is happend in this alignment for this constitutive junction
                        if junction_from_pos_index != -1 and junction_to_pos_index != -1:
                                assert junction_from_pos_index == junction_to_pos_index
                                self.CJ_counts[gene_id][CJ_number]["CJ_retained_reads"] += 1
        
                        # Junction is spliced. The upstream region of constituive exonic region is covered by the alignment but 
                        # the downstream region of constitutive intronic region is not covered by the alignment.
                        # If the alignment splices the CIR and jumps to another CER, we consider this is a splicing junction alignment for this constitutive junction.
                        elif junction_from_pos_index != -1 and CJ_type == "5'_splice_junction":
                                if gene_strand == "+":
                                        if junction_from_pos_index != len(end_list) - 1:
                                                # alt_next_pos = HTSeq.GenomicPosition(gene_chrom, start_list[junction_from_pos_index + 1], gene_strand)
                                                alt_next_pos_list = []
                                                next_fragment_length = end_list[junction_from_pos_index + 1] - start_list[junction_from_pos_index + 1] + 1
                                                for i in range(min(overlap, next_fragment_length)):
                                                        pos = HTSeq.GenomicPosition(gene_chrom, start_list[junction_from_pos_index + 1] + i, gene_strand)
                                                        alt_next_pos_list.append(pos)
                                        else:
                                                return
                                elif gene_strand == "-":
                                        if junction_from_pos_index != 0:
                                                # alt_next_pos = HTSeq.GenomicPosition(gene_chrom, end_list[junction_from_pos_index - 1], gene_strand)
                                                alt_next_pos_list = []
                                                next_fragment_length = end_list[junction_from_pos_index - 1] - start_list[junction_from_pos_index - 1] + 1
                                                for i in range(min(overlap, next_fragment_length)):
                                                        pos = HTSeq.GenomicPosition(gene_chrom, end_list[junction_from_pos_index - 1] - i, gene_strand)
                                                        alt_next_pos_list.append(pos)
                                        else:
                                                return 
                                
                                # If all the pos in alt_next_pos_list are located in constituitive exonic region, we count this as a splicing junction read.
                                if self.is_spliced_read_entirely_in_CER(alt_next_pos_list):
                                        self.CJ_counts[gene_id][CJ_number]["CJ_spliced_reads"] += 1                                      
        
                        # Junction is spliced. The upstream region of constituive intronic region is not covered by the alignment and  
                        # the downstream region of constitutive exonic region is covered by the alignment.
                        # If the alignment splices the CIR and jumps from another CER in upstream direction, we consider this is a splicing junction alignment for this constitutive junction.
                        elif junction_to_pos_index != -1 and CJ_type == "3'_splice_junction":
                                if gene_strand == "+":
                                        if junction_to_pos_index != 0:
                                                # alt_prev_pos = HTSeq.GenomicPosition(gene_chrom, end_list[junction_to_pos_index - 1], gene_strand)
                                                alt_prev_pos_list = []
                                                prev_fragment_length = end_list[junction_to_pos_index - 1] - start_list[junction_to_pos_index - 1] + 1
                                                for i in range(min(overlap, prev_fragment_length)):
                                                        pos = HTSeq.GenomicPosition(gene_chrom, end_list[junction_to_pos_index - 1] - i, gene_strand)
                                                        alt_prev_pos_list.append(pos)
                                        else:
                                                return
                                elif gene_strand == "-":
                                        if junction_to_pos_index != len(end_list) - 1:
                                                # alt_prev_pos = HTSeq.GenomicPosition(gene_chrom, start_list[junction_to_pos_index + 1], gene_strand)
                                                alt_prev_pos_list = []
                                                prev_fragment_length = end_list[junction_to_pos_index + 1] - start_list[junction_to_pos_index + 1] + 1
                                                for i in range(min(overlap, prev_fragment_length)):
                                                        pos = HTSeq.GenomicPosition(gene_chrom, start_list[junction_to_pos_index + 1] + i, gene_strand)
                                                        alt_prev_pos_list.append(pos)
                                        else:
                                                return
                                        
                                # If all the pos in alt_prev_pos_list are located in constituitive exonic region, we count this as a splicing junction read.
                                if self.is_spliced_read_entirely_in_CER(alt_prev_pos_list):
                                        self.CJ_counts[gene_id][CJ_number]["CJ_spliced_reads"] += 1 
                                        
        def assign_read_to_CIR(self, alt_iv_seq, CIR):
                gene_id, gene_chrom, gene_strand, CIR_number = CIR.attr["gene_id"], CIR.iv.chrom, CIR.iv.strand, CIR.attr["constitutive_intronic_region_number"]
                upstream_CJ = self.gene_CIR_associated_CJ_database[gene_id][CIR_number]['upstream_constitutive_junction']
                downstream_CJ = self.gene_CIR_associated_CJ_database[gene_id][CIR_number]['downstream_constitutive_junction']
                overlap = self.params['minoverlap']   
                
                start_list, end_list = self.generate_bisect_list(alt_iv_seq)
                
                if gene_strand == "+":
                        upstream_junction_from_pos = upstream_CJ.iv.start - 1
                        upstream_junction_from_pos_list = list(range(upstream_junction_from_pos - overlap + 1, upstream_junction_from_pos + 1))
                        
                        downstream_junction_to_pos = downstream_CJ.iv.start
                        downstream_junction_to_pos_list = list(range(downstream_junction_to_pos, downstream_junction_to_pos + overlap))
                        
                elif gene_strand == "-":
                        upstream_junction_from_pos = upstream_CJ.iv.start + 1
                        upstream_junction_from_pos_list = list(range(upstream_junction_from_pos, upstream_junction_from_pos + overlap))
                        
                        downstream_junction_to_pos = downstream_CJ.iv.start
                        downstream_junction_to_pos_list = list(range(downstream_junction_to_pos - overlap + 1, downstream_junction_to_pos + 1))   
                        
                upstream_junction_from_pos_index = self.find_pos_in_bisect_list(upstream_junction_from_pos_list, start_list, end_list)
                downstream_junction_to_pos_index = self.find_pos_in_bisect_list(downstream_junction_to_pos_list, start_list, end_list)
                        
                if upstream_junction_from_pos_index != -1 and downstream_junction_to_pos_index != -1:
                        if gene_strand == "+" and downstream_junction_to_pos_index - upstream_junction_from_pos_index == 1 and upstream_junction_from_pos == end_list[upstream_junction_from_pos_index] and downstream_junction_to_pos == start_list[downstream_junction_to_pos_index]:
                                self.CIR_counts[gene_id][CIR_number]['CIR_spliced_reads'] += 1
                        elif gene_strand ==  "-" and downstream_junction_to_pos_index - upstream_junction_from_pos_index == -1 and upstream_junction_from_pos == start_list[upstream_junction_from_pos_index] and downstream_junction_to_pos == end_list[downstream_junction_to_pos_index]:
                                self.CIR_counts[gene_id][CIR_number]['CIR_spliced_reads'] += 1 
                                
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
                
                logging.info("Counting number of retained reads and spliced reads for each constitutive junction (CJ) and constitutive intronic region (CIR)")                        
                
                bamfile = HTSeq.BAM_Reader(self.params['altfile'])
                if self.params['readtype'] == "single":
                        for alt in bamfile:
                                # Consider the alignments that are aligned and uniquely mapped.
                                if alt.aligned and self.unique_aligned(alt) and re.match('chr', alt.iv.chrom):                                              
                                        alt_iv_seq = self.get_alt_iv(alt)    
                                        if self.is_read_in_gene_region(alt_iv_seq) and self.is_read_in_CIR_or_CER(alt_iv_seq):                                             
                                                gene_id = self.read_associated_gene(alt_iv_seq)
                                                for CJ in list(self.gene_CJ_database[gene_id].values()):
                                                        self.assign_read_to_CJ(alt_iv_seq, CJ)
                                                for CIR in list(self.gene_CIR_database[gene_id].values()):
                                                        self.assign_read_to_CIR(alt_iv_seq, CIR)
                                                        
                elif self.params['readtype'] == "paired":
                        for alt_first, alt_second in HTSeq.pair_SAM_alignments(bamfile):
                                if alt_first == None or alt_second == None:
                                        continue
                                if alt_first.aligned and self.unique_aligned(alt_first) and alt_second.aligned and self.unique_aligned(alt_second) and alt_first.iv.chrom == alt_second.iv.chrom and re.match('chr', alt_first.iv.chrom) and re.match('chr', alt_second.iv.chrom):
                                        alt_first_iv_seq, alt_second_iv_seq = self.get_pair_alt_iv(alt_first, alt_second)
                                        alt_iv_seq = self.combine_pair_iv_seq(alt_first_iv_seq, alt_second_iv_seq)
                                        
                                        if self.is_read_in_gene_region(alt_iv_seq) and self.is_read_in_CIR_or_CER(alt_iv_seq):   
                                                gene_id = self.read_associated_gene(alt_iv_seq)
                                                for CJ in list(self.gene_CJ_database[gene_id].values()):
                                                        self.assign_read_to_CJ(alt_iv_seq, CJ)
                                                for CIR in list(self.gene_CIR_database[gene_id].values()):
                                                        self.assign_read_to_CIR(alt_iv_seq, CIR)
                                                        
        @staticmethod
        def get_quantile_index(read_count_qantile_list, read_count):
                return bisect.bisect_right(read_count_qantile_list, read_count) - 1
                                                        
        def empirical_filter(self, row, read_count_qantile_list, filter_cutoff_quantile_list):
                if row.CIR_retained_reads == 0:
                        return row.intron_IRC
                else:
                        if row.max_percentage <= filter_cutoff_quantile_list[self.get_quantile_index(read_count_qantile_list, row.CIR_retained_reads)]:
                                return row.intron_IRC
                        else:
                                if row.max_percentage == row["5_percentage"]:
                                        return "NA (5'AS)"
                                elif row.max_percentage == row["3_percentage"]:
                                        return "NA (3'AS)"
                                                        
        def apply_5_3_unbalanced_filter(self, df, outlier=0.01):
                logging.info("Appling filter to remove fake intron retention events")
                
                df["gene_id"] = df.CIR_id.str[:-4]
                df["CIR_number"] = df.CIR_id.str[-3:]
                
                df["CIR_retained_reads"] = df["CIR_5'retained_reads"] + df["CIR_3'retained_reads"]
                df["5_percentage"] = df["CIR_5'retained_reads"] / df["CIR_retained_reads"]
                df["3_percentage"] = df["CIR_3'retained_reads"] / df["CIR_retained_reads"]
                df["max_percentage"] = df.apply(lambda row: max(row["5_percentage"], row["3_percentage"]), axis=1)       
                
                try:
                        read_count_qantile_list = []
                        filter_cutoff_quantile_list = []
                        for q in [0, 25, 50, 75]:
                                lower, upper = np.percentile(df[df.CIR_retained_reads > 0].CIR_retained_reads, q), np.percentile(df[df.CIR_retained_reads > 0].CIR_retained_reads, q+25)
                                max_percentage_cutoff = np.percentile(df.loc[(df.CIR_retained_reads > lower) & (df.CIR_retained_reads <= upper), 'max_percentage'], (1 - outlier) * 100)
                                
                                read_count_qantile_list.append(lower)
                                filter_cutoff_quantile_list.append(max_percentage_cutoff)
                        
                        df['intron_IRC'] = df.apply(lambda row: self.empirical_filter(row, read_count_qantile_list, filter_cutoff_quantile_list), axis=1)
                        
                        self.filtered_CIR_id_list = list(df[df.intron_IRC.isin(["NA (5'AS)", "NA (3'AS)"])].CIR_id)  
                except IndexError:
                        self.filtered_CIR_id_list = []

                logging.info("{} constitutive intronic regions (CIR) are unlikely to be intron retention events and are filtered".format(len(self.filtered_CIR_id_list)))
                
                return df.drop(labels=['gene_id', 'CIR_number', 'CIR_retained_reads', '5_percentage', '3_percentage', 'max_percentage'], axis=1)
                                                        
        def output_IRC_junction_level(self):                 
                logging.info("Calculating IRC for each constitutive junction (CJ)")
                
                IRC_junction_level_data = []
                for gene_id in sorted(self.CJ_counts.keys()):
                        for CJ_number in sorted(self.CJ_counts[gene_id].keys()):
                                CJ_id = gene_id + ":" + CJ_number
                                CJ_type = self.gene_CJ_database[gene_id][CJ_number].attr["constitutive_junction_type"]
                                CJ_retained_reads = self.CJ_counts[gene_id][CJ_number]["CJ_retained_reads"]
                                CJ_spliced_reads = self.CJ_counts[gene_id][CJ_number]["CJ_spliced_reads"]                               
                                
                                junction_IRC = np.divide(CJ_retained_reads * 1.0, CJ_retained_reads + CJ_spliced_reads)
                                
                                IRC_junction_level_dict = {"CJ_id": CJ_id,
                                                           "CJ_iv": self.CJ_id2iv[CJ_id],
                                                           "CJ_type": CJ_type,
                                                           "CJ_retained_reads": CJ_retained_reads,
                                                           "CJ_spliced_reads": CJ_spliced_reads,
                                                           "junction_IRC": junction_IRC}
                                
                                IRC_junction_level_data.append(IRC_junction_level_dict)                                
                                                                   
                self.IRC_junction_level_df = pd.DataFrame(IRC_junction_level_data, columns=["CJ_id", "CJ_iv", "CJ_type", "CJ_retained_reads", "CJ_spliced_reads", "junction_IRC"])                
                
                outfile = self.params['name'] + ".quant.IRC.junctions.txt" 
                outfile_fullpath = os.path.join(self.params['outdir'], outfile)
                logging.info("Writing junction level result to file: {}".format(outfile_fullpath))
                self.IRC_junction_level_df.to_csv(outfile_fullpath, index=None, sep='\t', na_rep="NA")
                
        def output_IRC_intron_level(self):                 
                logging.info("Calculating IRC for each constitutive intronic region (CIR)")
                
                IRC_intron_level_data = []
                for gene_id in sorted(self.CIR_counts.keys()):
                        for CIR_number in sorted(self.CIR_counts[gene_id].keys()):
                                CIR_id = gene_id + ":" + CIR_number
                                upstream_CJ_number = self.gene_CIR_database[gene_id][CIR_number].attr["upstream_constitutive_junction_number"]
                                downstream_CJ_number = self.gene_CIR_database[gene_id][CIR_number].attr["downstream_constitutive_junction_number"]
                                CIR_five_retained_reads = self.CIR_counts[gene_id][CIR_number]["CIR_5'retained_reads"] = self.CJ_counts[gene_id][upstream_CJ_number]["CJ_retained_reads"]
                                CIR_three_retained_reads = self.CIR_counts[gene_id][CIR_number]["CIR_3'retained_reads"] = self.CJ_counts[gene_id][downstream_CJ_number]["CJ_retained_reads"]
                                
                                CIR_spliced_reads = self.CIR_counts[gene_id][CIR_number]["CIR_spliced_reads"]                                
                                intron_IRC = np.divide((CIR_five_retained_reads + CIR_three_retained_reads) / 2.0, (CIR_five_retained_reads + CIR_three_retained_reads) / 2.0 + CIR_spliced_reads)
                                                                
                                IRC_intron_level_dict = {"CIR_id": CIR_id,
                                                         "CIR_iv": self.CIR_id2iv[CIR_id],
                                                         "CIR_5'retained_reads": CIR_five_retained_reads,
                                                         "CIR_3'retained_reads": CIR_three_retained_reads,
                                                         "CIR_spliced_reads": CIR_spliced_reads,
                                                         "intron_IRC": intron_IRC}
                                
                                IRC_intron_level_data.append(IRC_intron_level_dict)                                
                                                                   
                self.IRC_intron_level_df = pd.DataFrame(IRC_intron_level_data, columns=["CIR_id", "CIR_iv", "CIR_5'retained_reads", "CIR_3'retained_reads", "CIR_spliced_reads", "intron_IRC"])                
                
                if self.filter:
                        self.IRC_intron_level_df = self.apply_5_3_unbalanced_filter(self.IRC_intron_level_df)
                
                outfile = self.params['name'] + ".quant.IRC.introns.txt" 
                outfile_fullpath = os.path.join(self.params['outdir'], outfile)
                logging.info("Writing intron level result to file: {}".format(outfile_fullpath))
                self.IRC_intron_level_df.to_csv(outfile_fullpath, index=None, sep='\t', na_rep="NA")
                
        def filter_CIR_id(self, filtered_CIR_id_list):
                for CIR_id in filtered_CIR_id_list:
                        gene_id, CIR_number = CIR_id[:-4], CIR_id[-3:]
                        del self.CIR_counts[gene_id][CIR_number]                 
                
        def output_IRC_gene_level(self, filtered_CIR_id_list=None):                 
                logging.info("Counting number of retained reads and spliced reads for each gene")        
                logging.info("Calculating IRC for each gene")
                
                if filtered_CIR_id_list is None:
                        filtered_CIR_id_list = self.filtered_CIR_id_list
                self.filter_CIR_id(filtered_CIR_id_list)                
                
                IRC_gene_level_data = []
                for gene_id in sorted(self.CIR_counts.keys()):
                        gene_retained_reads = self.gene_counts["gene_retained_reads"][gene_id] = sum([(self.CIR_counts[gene_id][CIR_number]["CIR_5'retained_reads"] + self.CIR_counts[gene_id][CIR_number]["CIR_3'retained_reads"]) / 2.0 for CIR_number in list(self.CIR_counts[gene_id].keys())])
                        gene_spliced_reads = self.gene_counts["gene_spliced_reads"][gene_id] = sum([self.CIR_counts[gene_id][CIR_number]["CIR_spliced_reads"] for CIR_number in list(self.CIR_counts[gene_id].keys())])
                        
                        gene_IRC = np.divide(gene_retained_reads * 1.0, gene_retained_reads + gene_spliced_reads)
                        
                        IRC_gene_level_dict = {"gene_id": gene_id,
                                               "gene_iv": self.gene_id2iv[gene_id],
                                               "gene_retained_reads": gene_retained_reads,
                                               "gene_spliced_reads": gene_spliced_reads,
                                               "gene_IRC": gene_IRC}
                
                        IRC_gene_level_data.append(IRC_gene_level_dict)                                    
                                                                          
                self.IRC_gene_level_df = pd.DataFrame(IRC_gene_level_data, columns=["gene_id", "gene_iv", "gene_retained_reads", "gene_spliced_reads", "gene_IRC"])                
                
                outfile = self.params['name'] + ".quant.IRC.genes.txt" 
                outfile_fullpath = os.path.join(self.params['outdir'], outfile)
                logging.info("Writing gene level result to file: {}".format(outfile_fullpath))
                self.IRC_gene_level_df.to_csv(outfile_fullpath, index=None, sep='\t', na_rep="NA")
                
        def output_IRC_genome_wide(self):
                total_retained_reads = sum(self.gene_counts["gene_retained_reads"].values())
                total_spliced_reads = sum(self.gene_counts["gene_spliced_reads"].values())
                total_reads = total_retained_reads + total_spliced_reads
                
                IRC_genome_wide = total_retained_reads * 1.0 / total_reads
                
                logging.info("Genome-wide intron retention statistics of the library:")
                print(("\tTotal retained reads: %i (%.2f%%)" % (total_retained_reads, total_retained_reads * 100.0 / total_reads)))
                print(("\tTotal spliced reads: %i (%.2f%%)" % (total_spliced_reads, total_spliced_reads * 100.0 / total_reads)))
                print(("\tTotal junction reads (retained reads + spliced reads): %i" % (total_retained_reads + total_spliced_reads)))
                print(("\tGenome-wide IRC: %f" % IRC_genome_wide))
