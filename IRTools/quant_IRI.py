import os
import logging
import re
import urllib.request, urllib.parse, urllib.error
import collections
import itertools
import networkx as nx
import numpy as np
import pandas as pd
import bisect
import HTSeq
import warnings
from functools import reduce

class IRI_quant(object):       
        def __init__(self, args):
                self.params = args.__dict__.copy()
                self.stranded = self.is_stranded(self.params['libtype'])
                
                self.gtffile = self.load_gtffile()
                
                self.gene_id2iv = self.get_gene_iv()
                self.valid_genes = set(self.gene_id2iv.keys())
                
                self.CIR_id2iv = self.get_CIR_iv()
                
                self.gene_map_score = self.init_mappability_GenomicArray(map_score_cutoff = 0.1)   
                
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
        
        def init_mappability_GenomicArray(self, map_score_cutoff):
                gene_map_score = HTSeq.GenomicArray("auto", stranded=self.stranded, typecode="i")              
                
                mapfile = self.params['mapfile']
                if mapfile:
                        from bx.bbi.bigwig_file import BigWigFile

                        if re.search('bigWig$', mapfile):
                                mapfile_data = BigWigFile(open(mapfile, "rb")) # [edit] added binary option
                        elif mapfile == "hg19" or mapfile == "mm9":
                                downloaded_mapfile = self.download_mappability_file_by_species_name(mapfile)
                                mapfile_data = BigWigFile(open(downloaded_mapfile, "rb")) # [edit] added binary option
                        else:
                                raise Exception("\"{}\" is neither a bigWig file nor a supported species name (hg19 or mm9).".format(mapfile))
                                                        
                        for feature in self.gtffile:
                                if feature.type == "constitutive_exonic_region":
                                        gene_map_score[feature.iv] = 1
        
                                elif feature.type == "constitutive_intronic_region":
                                        CIR_iv = feature.iv
                                        for start, end, score in mapfile_data.get(CIR_iv.chrom.encode('utf-8'), CIR_iv.start, CIR_iv.end): # [edit] made chrom binary
                                                score = 0 if score < map_score_cutoff else 1
                                                iv = HTSeq.GenomicInterval(CIR_iv.chrom, start, end, CIR_iv.strand)
                                                gene_map_score[iv] = score                                         
                        
                else:
                        for feature in self.gtffile:
                                if feature.type == "constitutive_exonic_region" or feature.type == "constitutive_intronic_region":
                                        gene_map_score[feature.iv] = 1                  
                return gene_map_score
        
        def get_CIR_effective_length(self):
                CIR_effective_length = collections.defaultdict( lambda: collections.Counter() )	
                for feature in self.gtffile:
                        gene_id = feature.attr["gene_id"]
                        if feature.type == "constitutive_intronic_region" and gene_id in self.valid_genes:
                                CIR_number = feature.attr["constitutive_intronic_region_number"]
                                CIR_effective_length[gene_id][CIR_number] = self.get_effective_length(feature.iv, self.gene_map_score)
                return CIR_effective_length
        
        def get_CER_length(self):
                CER_length = collections.defaultdict( lambda: collections.Counter() )
                for feature in self.gtffile:
                        gene_id = feature.attr["gene_id"]
                        if feature.type == "constitutive_exonic_region" and gene_id in self.valid_genes:
                                CER_number = feature.attr["constitutive_exonic_region_number"]
                                CER_length[gene_id][CER_number] = feature.iv.length
                return CER_length
        
        @staticmethod
        def get_effective_length(alt_iv_seq, gene_map_score):
                if type(alt_iv_seq) == list:
                        read_length = 0
                        for alt_iv in alt_iv_seq:
                                read_length += int(sum(gene_map_score[alt_iv]))
                        return read_length
                else:
                        return int(sum(gene_map_score[alt_iv_seq]))   
                
        def init_GenomicArrayOfSets_and_Counter_for_quant_IRI(self):
                genes = HTSeq.GenomicArrayOfSets("auto", stranded=self.stranded)
                gene_region = HTSeq.GenomicArrayOfSets("auto", stranded=self.stranded)                
                counts = collections.defaultdict( lambda:  collections.defaultdict( lambda:  collections.Counter()))  
                
                for feature in self.gtffile:
                        gene_id = feature.attr["gene_id"]
                        if feature.type == "gene_region":
                                gene_region[feature.iv] += gene_id
                        elif feature.type == "constitutive_exonic_region" and gene_id in self.valid_genes:
                                CER_number = feature.attr["constitutive_exonic_region_number"]
                                genes[feature.iv] += (gene_id, feature.type, CER_number)
                                counts[gene_id]["constitutive_exonic_region"][CER_number] = 0
                        elif feature.type == "constitutive_intronic_region" and gene_id in self.valid_genes:
                                CIR_number = feature.attr["constitutive_intronic_region_number"]
                                if self.CIR_effective_length[gene_id][CIR_number] > 0: 
                                        genes[feature.iv] += (gene_id, feature.type, CIR_number) 
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
        
        def init_GenomicArrayOfSets_and_Counter_for_bin_filter(self):
                bins = HTSeq.GenomicArrayOfSets("auto", stranded=self.stranded)                   
                bin_counts = collections.defaultdict( lambda:  collections.defaultdict( lambda:  collections.Counter()))
                        
                for feature in self.gtffile:
                        gene_id = feature.attr["gene_id"]
                        if feature.type == "constitutive_intronic_region" and gene_id in self.valid_genes:
                                CIR_number = feature.attr["constitutive_intronic_region_number"]
                                CIR_effective_length = self.CIR_effective_length[gene_id][CIR_number]
                                if CIR_effective_length >= self.num_bins:
                                        bins_iv_list = self.binnize(feature.iv, CIR_effective_length, self.gene_map_score, self.num_bins)
                                        for bin_number in range(self.num_bins):
                                                bins[bins_iv_list[bin_number]] += (gene_id, CIR_number, bin_number) 
                                                bin_counts[gene_id][CIR_number][bin_number] = 0 
                                else:
                                        for bin_number in range(self.num_bins):
                                                bin_counts[gene_id][CIR_number][bin_number] = np.nan                                         
                return bins, bin_counts
        
        def get_constitutive_junction_graph(self):
                G = collections.defaultdict( lambda: nx.Graph() )
                for feature in self.gtffile:
                        gene_id = feature.attr["gene_id"]
                        if feature.type == "constitutive_junction":
                                G[gene_id].add_edge(feature.attr['upstream'], feature.attr['downstream'], CJ_number=feature.attr['constitutive_junction_number'])                                   
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
        
        def assign_read_to_region(self, alt_iv_seq):
                alt_read_length = self.get_effective_length(alt_iv_seq, self.gene_map_score)
                if alt_read_length > 0:
                        for alt_iv in alt_iv_seq:
                                for iv, step_set in self.genes[alt_iv].steps():
                                        if step_set:
                                                gene_id, feature_type, region_number = list(step_set)[0]
                                                self.counts[gene_id][feature_type][region_number] += self.get_effective_length(iv, self.gene_map_score) * 1.0 / alt_read_length 
                                        
        def assign_read_to_bin_filter(self, alt_iv_seq):
                alt_read_length = self.get_effective_length(alt_iv_seq, self.gene_map_score)
                if alt_read_length > 0:
                        for alt_iv in alt_iv_seq:
                                for iv, step_set in self.bins[alt_iv].steps():
                                        if step_set:
                                                gene_id, CIR_number, bin_number = list(step_set)[0]
                                                self.bin_counts[gene_id][CIR_number][bin_number] += self.get_effective_length(iv, self.gene_map_score) * 1.0 / alt_read_length  
        
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
		bamfile = HTSeq.BAM_Reader(self.params['altfile'])
		# Single end
		if self.params['readtype'] == "single":
			for alt in bamfile:
				# Consider the alignments that are aligned and uniquely mapped.
				if alt.aligned and self.unique_aligned(alt) and re.match('chr', alt.iv.chrom):
					self.total_read_count += 1                                               
					alt_iv_seq = self.get_alt_iv(alt)
								
					# Eligible alignments are those mapped into one gene's either constitutive exonic region (CER) or constitutive intronic region (CIR).
					# For each eligible alignment, we count by fraction of length. i.e. If an alt has 50 bps, 30 bps in CER "001", 20 bps in CIR "001". Then, count in CER "001" is 0.6,
					# and count in CIR "001" is 0.4. (IRI is considered in intron level, so count is distributed in intron level)                                                
					if self.is_read_in_gene_region(alt_iv_seq) and self.is_read_in_CIR_or_CER(alt_iv_seq):    
						self.assign_read_to_region(alt_iv_seq)                                   
						if self.bin_filter:
							self.assign_read_to_bin_filter(alt_iv_seq)
							
		elif self.params['readtype'] == "paired":
			for alt_first, alt_second in HTSeq.pair_SAM_alignments(bamfile):
				if alt_first == None or alt_second == None:
					continue
				if alt_first.aligned and self.unique_aligned(alt_first) and alt_second.aligned and self.unique_aligned(alt_second) and alt_first.iv.chrom == alt_second.iv.chrom and re.match('chr', alt_first.iv.chrom) and re.match('chr', alt_second.iv.chrom):
					self.total_read_count += 1   
					alt_first_iv_seq, alt_second_iv_seq = self.get_pair_alt_iv(alt_first, alt_second)
					alt_iv_seq = self.combine_pair_iv_seq(alt_first_iv_seq, alt_second_iv_seq)
					
					if self.is_read_in_gene_region(alt_iv_seq) and self.is_read_in_CIR_or_CER(alt_iv_seq):    
						self.assign_read_to_region(alt_iv_seq)                                   
						if self.bin_filter:
							self.assign_read_to_bin_filter(alt_iv_seq)    
							
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
                
                df['gene_id'] = df.CIR_id.str[:-4]
                df['CIR_number'] = df.CIR_id.str[-3:]
                for i in range(self.num_bins):
                        bin_read_count_column = 'bin{}_read_count'.format(i)
                        df[bin_read_count_column] = df.apply(lambda row: self.bin_counts[row['gene_id']][row['CIR_number']][i], axis=1)
                        
                        bin_percentage_column = 'bin{}_percentage'.format(i)
                        df[bin_percentage_column] = df[bin_read_count_column] / df['CIR_read_count']
                
                bin_read_count_columns = ['bin{}_read_count'.format(i) for i in range(self.num_bins)]        
                bin_percentage_columns = ['bin{}_percentage'.format(i) for i in range(self.num_bins)]
                df['bin_max_percentage'] = df.apply(lambda row: max(row[bin_percentage_columns]), axis=1)                
                
                try: 
                        read_count_qantile_list = []
                        bin_filter_cutoff_quantile_list = []
                        bin_filter_df = df[pd.notnull(df.bin_max_percentage)]
                        for q in [0, 25, 50, 75]:
                                lower, upper = np.percentile(bin_filter_df[bin_filter_df.CIR_read_count > 0].CIR_read_count, q), np.percentile(bin_filter_df[bin_filter_df.CIR_read_count > 0].CIR_read_count, q+25)
                                bin_max_percentage_cutoff = np.percentile(bin_filter_df.loc[(bin_filter_df.CIR_read_count > lower) & (bin_filter_df.CIR_read_count <= upper), 'bin_max_percentage'], (1 - outlier) * 100)
                                
                                read_count_qantile_list.append(lower)
                                bin_filter_cutoff_quantile_list.append(bin_max_percentage_cutoff)
              
                        df['intron_IRI'] = df.apply(lambda row: self.empirical_bin_filter(row, read_count_qantile_list, bin_filter_cutoff_quantile_list), axis=1)   
                        
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
                                CIR_RPKM = CIR_read_count / (CIR_effective_length / 1000.0) / (self.total_read_count / 1000000.0)     
                                
                                adjacent_CER_list = consitutive_junction_graph.neighbors(CIR)
                                adjacent_CER_number_list = [item[-3:] for item in adjacent_CER_list]
                                adjacent_CER_read_count = 0
                                adjacent_CER_length = 0 
                                for adjacent_CER_number in adjacent_CER_number_list:
                                        adjacent_CER_read_count += self.counts[gene_id]["constitutive_exonic_region"][adjacent_CER_number]
                                        adjacent_CER_length += self.CER_length[gene_id][adjacent_CER_number]
                                adjacent_CER_RPKM = adjacent_CER_read_count / (adjacent_CER_length / 1000.0) / (self.total_read_count / 1000000.0)
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
                        
                        gene_CIR_RPKM = gene_CIR_read_count / (gene_CIR_effective_length / 1000.0) / (self.total_read_count / 1000000.0)
                        gene_CER_RPKM = gene_CER_read_count / (gene_CER_length / 1000.0) / (self.total_read_count / 1000000.0)
                        
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
                
                IRI_genome_wide = (total_CIR_read_count * 1.0 / total_CIR_effective_length) / (total_CER_read_count * 1.0 / total_CER_length)
                
                logging.info("Genome-wide intron retention statistics of the library:")
                print(("\tTotal CER read count: %f (%.2f%%)" % (total_CER_read_count, total_CER_read_count * 100.0 / self.total_read_count)))
                print(("\tTotal CIR read count: %f (%.2f%%)" % (total_CIR_read_count, total_CIR_read_count * 100.0 / self.total_read_count)))
                print(("\tTotal read count: %i" % self.total_read_count))
                print(("\tGenome-wide intron retention index: %f" % IRI_genome_wide))                
                
                
             
        
                        
                
                
                
                
