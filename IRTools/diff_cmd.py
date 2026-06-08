# Developer Note: IRTools diff works by first filtering out NA and inf
# values using a filtering algorithm dependent on whether the analysis
# is specified as paired or unpaired. In paired mode, rows missing from
# a replicate quant output are treated as 0 so replicate pairs stay
# aligned. In unpaired mode, only observed rows are collected.

import sys
import os
import logging
import re
import pandas as pd
import subprocess
import collections
import copy
import networkx as nx
import HTSeq
import math
import warnings
import scipy.stats
import statsmodels.stats.multitest
from statistics import mean
from IRTools.quant_IRI import IRI_quant
from IRTools.quant_IRC import IRC_quant

LOG2FC_PSEUDOCOUNT = 1e-6

def parse_replicate_names(value):
        return [item.strip() for item in value.split(',') if item.strip()]

def validate_diff_replicates(args):
        s1files = parse_replicate_names(args.s1files)
        s2files = parse_replicate_names(args.s2files)
        analysis_label = "RT" if args.quanttype == "RT" else "IR"

        if len(s1files) < 1 or len(s2files) < 1:
                logging.info("Run Aborted: Differential %s analysis requires at least one sample per group. Please check input.", analysis_label)
                sys.exit(1)

        if args.analysistype == "P" and len(s1files) != len(s2files):
                logging.info("Run Aborted: Samples must have the same number of replicates for paired analysis. Please check input.")
                sys.exit(1)

        args.s1files = ",".join(s1files)
        args.s2files = ",".join(s2files)

def read_diff_metric_file(indir, sample_name, suffix, value_index, id_header):
        file_path = os.path.join(indir, sample_name + suffix)
        values = {}
        with open(file_path, "r") as data:
                for line in data:
                        if line == "\n":
                                continue
                        fields = line.strip("\n").split()
                        if not fields or fields[0] == id_header:
                                continue
                        if len(fields) <= value_index:
                                raise ValueError("File {} does not have column index {} on line: {}".format(file_path, value_index, line.rstrip("\n")))
                        values[fields[0]] = fields[value_index]
        return values

def generate_diff_input_dict(params, suffix, value_index, id_header):
        s1files = parse_replicate_names(params['s1files'])
        s2files = parse_replicate_names(params['s2files'])
        s1_values_by_sample = {
                sample_name: read_diff_metric_file(params['indir'], sample_name, suffix, value_index, id_header)
                for sample_name in s1files
        }
        s2_values_by_sample = {
                sample_name: read_diff_metric_file(params['indir'], sample_name, suffix, value_index, id_header)
                for sample_name in s2files
        }

        feature_ids = set()
        for values_by_sample in (s1_values_by_sample, s2_values_by_sample):
                for sample_values in values_by_sample.values():
                        feature_ids.update(sample_values.keys())

        input_dict = {}
        paired_mode = params['analysistype'] == "P"
        descriptive_mode = len(s1files) < 2 or len(s2files) < 2
        for feature_id in feature_ids:
                if paired_mode or descriptive_mode:
                        s1_values = [s1_values_by_sample[sample_name].get(feature_id, "0") for sample_name in s1files]
                        s2_values = [s2_values_by_sample[sample_name].get(feature_id, "0") for sample_name in s2files]
                else:
                        s1_values = [s1_values_by_sample[sample_name][feature_id] for sample_name in s1files if feature_id in s1_values_by_sample[sample_name]]
                        s2_values = [s2_values_by_sample[sample_name][feature_id] for sample_name in s2files if feature_id in s2_values_by_sample[sample_name]]
                input_dict[feature_id] = (s1_values, s2_values)
        return input_dict

def diff_descriptive_mode(params):
        return len(parse_replicate_names(params['s1files'])) < 2 or len(parse_replicate_names(params['s2files'])) < 2

def numeric_diff_value(value):
        try:
                number = float(value)
        except (TypeError, ValueError):
                return None
        if math.isnan(number) or math.isinf(number):
                return None
        return number

def collect_numeric_diff_values(values_s1, values_s2, paired_mode):
        num_s1 = []
        num_s2 = []
        if paired_mode:
                for value_s1, value_s2 in zip(values_s1, values_s2):
                        number_s1 = numeric_diff_value(value_s1)
                        number_s2 = numeric_diff_value(value_s2)
                        if number_s1 is not None and number_s2 is not None:
                                num_s1.append(number_s1)
                                num_s2.append(number_s2)
        else:
                num_s1 = [number for number in (numeric_diff_value(value) for value in values_s1) if number is not None]
                num_s2 = [number for number in (numeric_diff_value(value) for value in values_s2) if number is not None]
        return num_s1, num_s2

def count_distinct_vals(num_s1, num_s2):
        distinct_vals = []
        for x in num_s1:
                if x not in distinct_vals:
                        distinct_vals.append(x)
        for x in num_s2:
                if x not in distinct_vals:
                        distinct_vals.append(x)
        return len(distinct_vals)

def calc_change_value(num_s1, num_s2, params):
        s1_mean = mean(num_s1)
        s2_mean = mean(num_s2)
        if params.get("method", "minus") == "log2fc":
                return math.log2((s2_mean + LOG2FC_PSEUDOCOUNT) / (s1_mean + LOG2FC_PSEUDOCOUNT))
        return s2_mean - s1_mean

def change_column_name(base_name, params):
        if params.get("method", "minus") == "log2fc":
                return base_name + "_log2FC"
        return base_name + "_difference"

def run_ttest(num_s1, num_s2, paired_mode):
        with warnings.catch_warnings():
                warnings.filterwarnings("ignore", category=RuntimeWarning)
                if paired_mode:
                        return scipy.stats.ttest_rel(num_s1, num_s2)[1]
                return scipy.stats.ttest_ind(num_s1, num_s2)[1]

def run_diff_analysis(input_dict, params):
        descriptive_mode = diff_descriptive_mode(params)
        paired_mode = params["analysistype"] == "P"
        pending_results = []
        pval_list = []

        for id in sorted(input_dict.keys()):
                values_s1 = input_dict[id][0]
                values_s2 = input_dict[id][1]
                num_s1, num_s2 = collect_numeric_diff_values(values_s1, values_s2, paired_mode)
                if len(num_s1) == 0 or len(num_s2) == 0:
                        continue

                diff = calc_change_value(num_s1, num_s2, params)
                if descriptive_mode:
                        pending_results.append([id, "NA", "NA", values_s1, values_s2, diff])
                        continue

                if count_distinct_vals(num_s1, num_s2) == 1:
                        pending_results.append([id, "NA", "NA", values_s1, values_s2, diff])
                        continue
                if len(num_s1) < 2 or len(num_s2) < 2:
                        pending_results.append([id, "NA", "NA", values_s1, values_s2, diff])
                        continue
                pval = run_ttest(num_s1, num_s2, paired_mode)
                if math.isnan(pval):
                        pending_results.append([id, "NA", "NA", values_s1, values_s2, diff])
                        continue
                pending_results.append([id, pval, None, values_s1, values_s2, diff])
                pval_list.append(pval)

        if not descriptive_mode and pval_list:
                fdr_bool_list, fdr_pval_list = statsmodels.stats.multitest.fdrcorrection(pval_list)
                fdr_index = 0
                for result in pending_results:
                        if result[1] != "NA":
                                result[2] = fdr_pval_list[fdr_index]
                                fdr_index += 1

        return pending_results

class IRI_diff(object):
        def __init__(self, args):
                self.params = args.__dict__.copy()
                
                self.temp_dir = self.check_temp_dir(self.params['outdir']) 
                print("\tNote: Running \"IRTools diff\" will produce some intermediate files saved in directory: {}/".format(self.temp_dir))
                sys.stdout.flush()   
                
                self.logger = logging.getLogger()
        
        @staticmethod                        
        def check_temp_dir(outdir):
                temp_dir = os.path.join(outdir, "temp")
                if not os.path.exists( temp_dir ):
                        try:
                                os.makedirs( temp_dir ) 
                        except:
                                temp_dir = outdir
                                             
                return temp_dir                  
        
        @staticmethod
        def count_distinct_vals(num_IRI_S1, num_IRI_S2):
                distinct_vals = []
                for x in num_IRI_S1:
                        if x not in distinct_vals:
                                distinct_vals.append(x)
                for x in num_IRI_S2:
                        if x not in distinct_vals:
                                distinct_vals.append(x)
                return len(distinct_vals)
                                            
        def generate_input_intron_level(self):
                logging.info("Generating inputs for analysis for differential IR in intron level")

                self.input_intron_dict = generate_diff_input_dict(self.params, ".quant.IRI.introns.txt", 8, "CIR_id")

                input_file_path = os.path.join(self.temp_dir, self.params['name'] + ".diff.input.IRI.introns.txt")
                input_file = open(input_file_path, "w")
                input_file.write("CIR_id\tintron_IRI_S1\tintron_IRI_S2\n")
                for id in sorted(self.input_intron_dict.keys()):
                        input_file.write(id + "\t" + ",".join(self.input_intron_dict[id][0]) + "\t" + ",".join(self.input_intron_dict[id][1]) + "\n")
                input_file.close()
                
                logging.info("Intron level analysis inputs can be found in " + input_file_path)

        def generate_input_gene_level(self):
                logging.info("Generating inputs for analysis for differential IR in gene level")

                self.input_gene_dict = generate_diff_input_dict(self.params, ".quant.IRI.genes.txt", 8, "gene_id")

                input_file_path = os.path.join(self.temp_dir, self.params['name'] + ".diff.input.IRI.genes.txt")
                input_file = open(input_file_path, "w")
                input_file.write("gene_id\tgene_IRI_S1\tgene_IRI_S2\n")
                for id in sorted(self.input_gene_dict.keys()):
                        input_file.write(id + "\t" + ",".join(self.input_gene_dict[id][0]) + "\t" + ",".join(self.input_gene_dict[id][1]) + "\n")
                input_file.close()

                logging.info("Gene level analysis inputs can be found at " + input_file_path)

        def run_analysis_intron_level(self):
                logging.info("Running analysis for differential IR in intron level")

                results = run_diff_analysis(self.input_intron_dict, self.params)
                
                results_file_path = os.path.join(self.params['outdir'], self.params['name'] + ".diff.IRI.introns.txt")
                results_file = open(results_file_path, "w")
                results_file.write("CIR_id\tPValue\tFDR\tintron_IRI_S1\tintron_IRI_S2\t{}\n".format(change_column_name("intron_IRI", self.params)))
                for id, pval, fdr, values_s1, values_s2, diff in results:
                        results_file.write(id + "\t" + str(pval) + "\t" + str(fdr) + "\t" + ",".join(values_s1) + "\t" + ",".join(values_s2) + "\t" + str(diff) + "\n")
                results_file.close()

                logging.info("Intron level differential IR results can be found in " + results_file_path)

        def run_analysis_gene_level(self):
                logging.info("Running analysis for differential IR in gene level")

                results = run_diff_analysis(self.input_gene_dict, self.params)

                results_file_path = os.path.join(self.params['outdir'], self.params['name'] + ".diff.IRI.genes.txt")
                results_file = open(results_file_path, "w")
                results_file.write("gene_id\tPValue\tFDR\tgene_IRI_S1\tgene_IRI_S2\t{}\n".format(change_column_name("gene_IRI", self.params)))
                for id, pval, fdr, values_s1, values_s2, diff in results:
                        results_file.write(id + "\t" + str(pval) + "\t" + str(fdr) + "\t" + ",".join(values_s1) + "\t" + ",".join(values_s2) + "\t" + str(diff) + "\n")
                results_file.close()

                logging.info("Gene level differential IR results can be found in " + results_file_path)

class IRC_diff(object):      
        def __init__(self, args):
                self.params = args.__dict__.copy()
                
                self.temp_dir = self.check_temp_dir(self.params['outdir']) 
                print("\tNote: Running \"IRTools diff\" will produce some intermediate files saved in directory: {}/".format(self.temp_dir))
                sys.stdout.flush()   
                
                self.logger = logging.getLogger()  
                
        @staticmethod                        
        def check_temp_dir(outdir):
                temp_dir = os.path.join(outdir, "temp")
                if not os.path.exists( temp_dir ):
                        try:
                                os.makedirs( temp_dir ) 
                        except:
                                temp_dir = outdir
                                             
                return temp_dir  
        
        @staticmethod
        def count_distinct_vals(num_IRC_S1, num_IRC_S2):
                distinct_vals = []
                for x in num_IRC_S1:
                        if x not in distinct_vals:
                                distinct_vals.append(x)
                for x in num_IRC_S2:
                        if x not in distinct_vals:
                                distinct_vals.append(x)
                return len(distinct_vals)
                                            
        def generate_input_intron_level(self):
                logging.info("Generating inputs for analysis for differential IR in intron level")

                self.input_intron_dict = generate_diff_input_dict(self.params, ".quant.IRC.introns.txt", 5, "CIR_id")

                input_file_path = os.path.join(self.temp_dir, self.params['name'] + ".diff.input.IRC.introns.txt")
                input_file = open(input_file_path, "w")
                input_file.write("CIR_id\tintron_IRC_S1\tintron_IRC_S2\n")
                for id in sorted(self.input_intron_dict.keys()):
                        input_file.write(id + "\t" + ",".join(self.input_intron_dict[id][0]) + "\t" + ",".join(self.input_intron_dict[id][1]) + "\n")
                input_file.close()

                logging.info("Intron level analysis inputs can be found in " + input_file_path)

        def generate_input_gene_level(self):
                logging.info("Generating inputs for analysis for differential IR in gene level")

                self.input_gene_dict = generate_diff_input_dict(self.params, ".quant.IRC.genes.txt", 4, "gene_id")

                input_file_path = os.path.join(self.temp_dir, self.params['name'] + ".diff.input.IRC.genes.txt")
                input_file = open(input_file_path, "w")
                input_file.write("gene_id\tgene_IRC_S1\tgene_IRC_S2\n")
                for id in sorted(self.input_gene_dict.keys()):
                        input_file.write(id + "\t" + ",".join(self.input_gene_dict[id][0]) + "\t" + ",".join(self.input_gene_dict[id][1]) + "\n")
                input_file.close()

                logging.info("Gene level analysis inputs can be found at " + input_file_path)

        def generate_input_junction_level(self):
                logging.info("Generating inputs for analysis for differential IR in junction level")

                self.input_junction_dict = generate_diff_input_dict(self.params, ".quant.IRC.junctions.txt", 5, "CJ_id")

                input_file_path = os.path.join(self.temp_dir, self.params['name'] + ".diff.input.IRC.junctions.txt")
                input_file = open(input_file_path, "w")
                input_file.write("CJ_id\tjunction_IRC_S1\tjunction_IRC_S2\n")
                for id in sorted(self.input_junction_dict.keys()):
                        input_file.write(id + "\t" + ",".join(self.input_junction_dict[id][0]) + "\t" + ",".join(self.input_junction_dict[id][1]) + "\n")
                input_file.close()

                logging.info("junction level analysis inputs can be found at " + input_file_path)

        def run_analysis_intron_level(self):
                logging.info("Running analysis for differential IR in intron level")

                results = run_diff_analysis(self.input_intron_dict, self.params)

                results_file_path = os.path.join(self.params['outdir'], self.params['name'] + ".diff.IRC.introns.txt")
                results_file = open(results_file_path, "w")
                results_file.write("CIR_id\tPValue\tFDR\tintron_IRC_S1\tintron_IRC_S2\t{}\n".format(change_column_name("intron_IRC", self.params)))
                for id, pval, fdr, values_s1, values_s2, diff in results:
                        results_file.write(id + "\t" + str(pval) + "\t" + str(fdr) + "\t" + ",".join(values_s1) + "\t" + ",".join(values_s2) + "\t" + str(diff) + "\n")
                results_file.close()

                logging.info("Intron level differential IR results can be found in " + results_file_path)

        def run_analysis_gene_level(self):
                logging.info("Running analysis for differential IR in gene level")

                results = run_diff_analysis(self.input_gene_dict, self.params)

                results_file_path = os.path.join(self.params['outdir'], self.params['name'] + ".diff.IRC.genes.txt")
                results_file = open(results_file_path, "w")
                results_file.write("gene_id\tPValue\tFDR\tgene_IRC_S1\tgene_IRC_S2\t{}\n".format(change_column_name("gene_IRC", self.params)))
                for id, pval, fdr, values_s1, values_s2, diff in results:
                        results_file.write(id + "\t" + str(pval) + "\t" + str(fdr) + "\t" + ",".join(values_s1) + "\t" + ",".join(values_s2) + "\t" + str(diff) + "\n")
                results_file.close()

                logging.info("Gene level differential IR results can be found in " + results_file_path)

        def run_analysis_junction_level(self):
                logging.info("Running analysis for differential IR in junction level")

                results = run_diff_analysis(self.input_junction_dict, self.params)

                results_file_path = os.path.join(self.params['outdir'], self.params['name'] + ".diff.IRC.junctions.txt")
                results_file = open(results_file_path, "w")
                results_file.write("CJ_id\tPValue\tFDR\tjunction_IRC_S1\tjunction_IRC_S2\t{}\n".format(change_column_name("junction_IRC", self.params)))
                for id, pval, fdr, values_s1, values_s2, diff in results:
                        results_file.write(id + "\t" + str(pval) + "\t" + str(fdr) + "\t" + ",".join(values_s1) + "\t" + ",".join(values_s2) + "\t" + str(diff) + "\n")
                results_file.close()

                logging.info("Junction level differential IR results can be found in " + results_file_path)

class RT_diff(object):
        def __init__(self, args):
                self.params = args.__dict__.copy()

                self.temp_dir = self.check_temp_dir(self.params['outdir'])
                print("\tNote: Running \"IRTools diff\" will produce some intermediate files saved in directory: {}/".format(self.temp_dir))
                sys.stdout.flush()

                self.logger = logging.getLogger()

        @staticmethod
        def check_temp_dir(outdir):
                temp_dir = os.path.join(outdir, "temp")
                if not os.path.exists(temp_dir):
                        try:
                                os.makedirs(temp_dir)
                        except:
                                temp_dir = outdir

                return temp_dir

        @staticmethod
        def count_distinct_vals(num_RT_S1, num_RT_S2):
                distinct_vals = []
                for x in num_RT_S1:
                        if x not in distinct_vals:
                                distinct_vals.append(x)
                for x in num_RT_S2:
                        if x not in distinct_vals:
                                distinct_vals.append(x)
                return len(distinct_vals)

        def generate_input_gene_level(self):
                logging.info("Generating inputs for analysis for differential RT in gene level")

                self.input_gene_dict = generate_diff_input_dict(self.params, ".quant.RT.genes.txt", 10, "gene_id")

                input_file_path = os.path.join(self.temp_dir, self.params['name'] + ".diff.input.RT.genes.txt")
                input_file = open(input_file_path, "w")
                input_file.write("gene_id\tgene_RT_S1\tgene_RT_S2\n")
                for id in sorted(self.input_gene_dict.keys()):
                        input_file.write(id + "\t" + ",".join(self.input_gene_dict[id][0]) + "\t" + ",".join(self.input_gene_dict[id][1]) + "\n")
                input_file.close()

                logging.info("Gene level analysis inputs can be found at " + input_file_path)

        def run_analysis_gene_level(self):
                logging.info("Running analysis for differential RT in gene level")

                results = run_diff_analysis(self.input_gene_dict, self.params)

                results_file_path = os.path.join(self.params['outdir'], self.params['name'] + ".diff.RT.genes.txt")
                results_file = open(results_file_path, "w")
                results_file.write("gene_id\tPValue\tFDR\tgene_RT_S1\tgene_RT_S2\t{}\n".format(change_column_name("gene_RT", self.params)))
                for id, pval, fdr, values_s1, values_s2, diff in results:
                        results_file.write(id + "\t" + str(pval) + "\t" + str(fdr) + "\t" + ",".join(values_s1) + "\t" + ",".join(values_s2) + "\t" + str(diff) + "\n")
                results_file.close()

                logging.info("Gene level differential RT results can be found in " + results_file_path)

def run(args):
        validate_diff_replicates(args)

        if args.quanttype == "IRI":
                IRI_differ = IRI_diff(args)
                IRI_differ.generate_input_intron_level()
                IRI_differ.generate_input_gene_level()
                IRI_differ.run_analysis_intron_level()
                IRI_differ.run_analysis_gene_level()

        elif args.quanttype == "IRC":
                IRC_differ = IRC_diff(args)
                IRC_differ.generate_input_intron_level()
                IRC_differ.generate_input_gene_level()
                IRC_differ.generate_input_junction_level()
                IRC_differ.run_analysis_intron_level()
                IRC_differ.run_analysis_gene_level()
                IRC_differ.run_analysis_junction_level()

        elif args.quanttype == "RT":
                RT_differ = RT_diff(args)
                RT_differ.generate_input_gene_level()
                RT_differ.run_analysis_gene_level()
