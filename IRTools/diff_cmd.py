import sys
import os
import logging
import pkg_resources
import re
import pandas as pd
import subprocess
import collections
import copy
import networkx as nx
import HTSeq
import math
import scipy.stats
import statsmodels.stats.multitest
from statistics import mean
from IRTools.quant_IRI import IRI_quant
from IRTools.quant_IRC import IRC_quant

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
        
        def run_IRI_quant_for_all_samples(self, args):
                logging.info("Perform \"IRTools quant\" for all replciates from both samples")
                self.logger.disabled = True
                
                filtered_CIR_id_list = []
                
                s1files_list = self.params['s1files'].split(',')
                s2files_list = self.params['s2files'].split(',')                

                IRI_intron_level_df_s1_list = []
                IRI_intron_level_df_s2_list = []
                
                counts_s1_list = []
                counts_s2_list = []
                
                total_read_count_s1_list = []
                total_read_count_s2_list = []                
                
                IRI_quanter = IRI_quant(args)
                IRI_quanter.params['outdir'] = self.temp_dir 
                
                for i, s1file in enumerate(s1files_list):
                        IRI_quanter.params['name'] = self.params['name'] + "_S1_R%d" % (i+1)
                        IRI_quanter.params['altfile'] = s1file
                        
                        IRI_quanter.quant()
                        IRI_quanter.output_IRI_intron_level()
                        IRI_intron_level_df_s1_list.append(IRI_quanter.IRI_intron_level_df)
                        
                        counts_s1_list.append(copy.deepcopy(IRI_quanter.counts))
                        total_read_count_s1_list.append(IRI_quanter.total_read_count)
                        filtered_CIR_id_list.extend(IRI_quanter.filtered_CIR_id_list)
                                               
                for i, s2file in enumerate(s2files_list):
                        IRI_quanter.params['name'] = self.params['name'] + "_S2_R%d" % (i+1)
                        IRI_quanter.params['altfile'] = s2file
                        
                        IRI_quanter.quant()
                        IRI_quanter.output_IRI_intron_level()
                        IRI_intron_level_df_s2_list.append(IRI_quanter.IRI_intron_level_df)
                        
                        counts_s2_list.append(copy.deepcopy(IRI_quanter.counts))
                        total_read_count_s2_list.append(IRI_quanter.total_read_count)
                        filtered_CIR_id_list.extend(IRI_quanter.filtered_CIR_id_list)
                        
                filtered_CIR_id_list = list(set(filtered_CIR_id_list))
                
                for i, IRI_intron_level_df in enumerate(IRI_intron_level_df_s1_list):
                        IRI_intron_level_df_s1_list[i] = IRI_intron_level_df[IRI_intron_level_df.CIR_id.isin(filtered_CIR_id_list) == False]
                
                for i, IRI_intron_level_df in enumerate(IRI_intron_level_df_s2_list):
                        IRI_intron_level_df_s2_list[i] = IRI_intron_level_df[IRI_intron_level_df.CIR_id.isin(filtered_CIR_id_list) == False]                     
                        
                self.IRI_intron_level_data = {'s1_data': IRI_intron_level_df_s1_list,
                                              's2_data': IRI_intron_level_df_s2_list} 
                
                IRI_gene_level_df_s1_list = []
                IRI_gene_level_df_s2_list = []                
                
                for i, s1file in enumerate(s1files_list):
                        IRI_quanter.params['name'] = self.params['name'] + "_S1_R%d" % (i+1)
                        
                        IRI_quanter.counts = counts_s1_list[i]   
                        IRI_quanter.total_read_count = total_read_count_s1_list[i]
                        IRI_quanter.output_IRI_gene_level(filtered_CIR_id_list)
                        IRI_gene_level_df_s1_list.append(IRI_quanter.IRI_gene_level_df)
                
                        self.logger.disabled = False
                        logging.info("Sample: " + IRI_quanter.params['name'])                   
                        IRI_quanter.output_IRI_genome_wide()
                        self.logger.disabled = True
                        
                for i, s2file in enumerate(s2files_list):
                        IRI_quanter.params['name'] = self.params['name'] + "_S2_R%d" % (i+1)
                        
                        IRI_quanter.counts = counts_s2_list[i]
                        IRI_quanter.total_read_count = total_read_count_s2_list[i]
                        IRI_quanter.output_IRI_gene_level(filtered_CIR_id_list)
                        IRI_gene_level_df_s2_list.append(IRI_quanter.IRI_gene_level_df)
                
                        self.logger.disabled = False
                        logging.info("Sample: " + IRI_quanter.params['name'])  
                        sys.stdout.flush()
                        IRI_quanter.output_IRI_genome_wide()  
                        self.logger.disabled = True
                        
                self.IRI_gene_level_data = {'s1_data': IRI_gene_level_df_s1_list,
                                            's2_data': IRI_gene_level_df_s2_list}
                
                self.logger.disabled = False
        
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

                temp_dict = {}

                for i in range(len(self.params['s1files'].split(','))):
                        file_path = os.path.join(self.temp_dir, self.params['name'] + "_S1_R%d.quant.IRI.introns.txt" % (i + 1))
                        data = open(file_path, "r")
                        lines = [x.strip("\n") for x in data if x != "\n"]
                        if not temp_dict:
                                for l in lines:
                                        temp_dict[l.split()[0]] = ([l.split()[8]], [])
                        else:
                                for l in lines:
                                        temp_dict[l.split()[0]][0].append(l.split()[8])
                        data.close()

                for i in range(len(self.params['s2files'].split(','))):
                        file_path = os.path.join(self.temp_dir, self.params['name'] + "_S2_R%d.quant.IRI.introns.txt" % (i + 1))
                        data = open(file_path, "r")
                        lines = [x.strip("\n") for x in data if x != "\n"]
                        for l in lines:
                                temp_dict[l.split()[0]][1].append(l.split()[8])
                        data.close()

                del temp_dict["CIR_id"]
                self.input_intron_dict = temp_dict

                input_file_path = os.path.join(self.temp_dir, self.params['name'] + ".diff.input.IRI.introns.txt")
                input_file = open(input_file_path, "w")
                input_file.write("CIR_id\tintron_IRI_S1\tintron_IRI_S2\n")
                for id in sorted(self.input_intron_dict.keys()):
                        input_file.write(id + "\t" + ",".join(self.input_intron_dict[id][0]) + "\t" + ",".join(self.input_intron_dict[id][1]) + "\n")
                input_file.close()
                
                logging.info("Intron level analysis inputs can be found in " + input_file_path)

        def generate_input_gene_level(self):
                logging.info("Generating inputs for analysis for differential IR in gene level")

                temp_dict = {}

                for i in range(len(self.params['s1files'].split(','))):
                        file_path = os.path.join(self.temp_dir, self.params['name'] + "_S1_R%d.quant.IRI.genes.txt" % (i + 1))
                        data = open(file_path, "r")
                        lines = [x.strip("\n") for x in data if x != "\n"]
                        if not temp_dict:
                                for l in lines:
                                        temp_dict[l.split()[0]] = ([l.split()[8]], [])
                        else:
                                for l in lines:
                                        temp_dict[l.split()[0]][0].append(l.split()[8])
                        data.close()

                for i in range(len(self.params['s2files'].split(','))):
                        file_path = os.path.join(self.temp_dir, self.params['name'] + "_S2_R%d.quant.IRI.genes.txt" % (i + 1))
                        data = open(file_path, "r")
                        lines = [x.strip("\n") for x in data if x != "\n"]
                        for l in lines:
                                temp_dict[l.split()[0]][1].append(l.split()[8])
                        data.close()

                del temp_dict["gene_id"]
                self.input_gene_dict = temp_dict

                input_file_path = os.path.join(self.temp_dir, self.params['name'] + ".diff.input.IRI.genes.txt")
                input_file = open(input_file_path, "w")
                input_file.write("gene_id\tgene_IRI_S1\tgene_IRI_S2\n")
                for id in sorted(self.input_gene_dict.keys()):
                        input_file.write(id + "\t" + ",".join(self.input_gene_dict[id][0]) + "\t" + ",".join(self.input_gene_dict[id][1]) + "\n")
                input_file.close()

                logging.info("Gene level analysis inputs can be found at " + input_file_path)

        def run_analysis_intron_level(self):
                logging.info("Running analysis for differential IR in intron level")

                filtered_introns = {}
                pval_list = []
                IRI_diff_list = []

                for id in sorted(self.input_intron_dict.keys()):
                        intron_IRI_S1 = self.input_intron_dict[id][0]
                        intron_IRI_S2 = self.input_intron_dict[id][1]
                        num_intron_IRI_S1 = []
                        num_intron_IRI_S2 = []
                        if self.params["analysistype"] == "P":
                                for i, val in enumerate(intron_IRI_S1):
                                        if val != "NA" and val != 'inf' and intron_IRI_S2[i] != "NA" and intron_IRI_S2[i] != 'inf':
                                                num_intron_IRI_S1.append(float(val))
                                                num_intron_IRI_S2.append(float(intron_IRI_S2[i]))
                        else:
                                num_intron_IRI_S1 = [float(x) for x in intron_IRI_S1 if x != "NA" and x != 'inf']
                                num_intron_IRI_S2 = [float(x) for x in intron_IRI_S2 if x != "NA" and x != 'inf']
                        if self.count_distinct_vals(num_intron_IRI_S1, num_intron_IRI_S2) == 1:
                                continue
                        if len(num_intron_IRI_S1) < 2 or len(num_intron_IRI_S2) < 2:
                                continue
                        if self.params["analysistype"] == "P":
                                pval = scipy.stats.ttest_rel(num_intron_IRI_S1, num_intron_IRI_S2)[1]
                        else:
                                pval = scipy.stats.ttest_ind(num_intron_IRI_S1, num_intron_IRI_S2)[1]
                        if math.isnan(pval):
                                continue
                        filtered_introns[id] = (intron_IRI_S1, intron_IRI_S2)
                        pval_list.append(pval)
                        diff = mean(num_intron_IRI_S2) - mean(num_intron_IRI_S1)
                        IRI_diff_list.append(diff)
                
                fdr_bool_list, fdr_pval_list = statsmodels.stats.multitest.fdrcorrection(pval_list)
                
                results_file_path = os.path.join(self.params['outdir'], self.params['name'] + ".diff.IRI.introns.txt")
                results_file = open(results_file_path, "w")
                results_file.write("CIR_id\tPValue\tFDR\tintron_IRI_S1\tintron_IRI_S2\tintron_IRI_difference\n")
                for i, id in enumerate(sorted(filtered_introns.keys())):
                        results_file.write(id + "\t" + str(pval_list[i]) + "\t" + str(fdr_pval_list[i]) + "\t" + ",".join(filtered_introns[id][0]) + "\t" + ",".join(filtered_introns[id][1]) + "\t" + str(IRI_diff_list[i]) + "\n")
                results_file.close()

                logging.info("Intron level differential IR results can be found in " + results_file_path)

        def run_analysis_gene_level(self):
                logging.info("Running analysis for differential IR in gene level")

                filtered_genes = {}
                pval_list = []
                IRI_diff_list = []

                for id in sorted(self.input_gene_dict.keys()):
                        gene_IRI_S1 = self.input_gene_dict[id][0]
                        gene_IRI_S2 = self.input_gene_dict[id][1]
                        num_gene_IRI_S1 = []
                        num_gene_IRI_S2 = []
                        if self.params["analysistype"] == "P":
                                for i, val in enumerate(gene_IRI_S1):
                                        if val != "NA" and val != 'inf' and gene_IRI_S2[i] != "NA" and gene_IRI_S2[i] != 'inf':
                                                num_gene_IRI_S1.append(float(val))
                                                num_gene_IRI_S2.append(float(gene_IRI_S2[i]))
                        else:
                                num_gene_IRI_S1 = [float(x) for x in gene_IRI_S1 if x != "NA" and x != 'inf']
                                num_gene_IRI_S2 = [float(x) for x in gene_IRI_S2 if x != "NA" and x != 'inf']
                        if self.count_distinct_vals(num_gene_IRI_S1, num_gene_IRI_S2) == 1:
                                continue
                        if len(num_gene_IRI_S1) < 2 or len(num_gene_IRI_S2) < 2:
                                continue
                        if self.params["analysistype"] == "P":
                                pval = scipy.stats.ttest_rel(num_gene_IRI_S1, num_gene_IRI_S2)[1]
                        else:
                                pval = scipy.stats.ttest_ind(num_gene_IRI_S1, num_gene_IRI_S2)[1]
                        if math.isnan(pval):
                                continue
                        filtered_genes[id] = (gene_IRI_S1, gene_IRI_S2)
                        pval_list.append(pval)
                        diff = mean(num_gene_IRI_S2) - mean(num_gene_IRI_S1)
                        IRI_diff_list.append(diff)

                fdr_bool_list, fdr_pval_list = statsmodels.stats.multitest.fdrcorrection(pval_list)

                results_file_path = os.path.join(self.params['outdir'], self.params['name'] + ".diff.IRI.genes.txt")
                results_file = open(results_file_path, "w")
                results_file.write("gene_id\tPValue\tFDR\tgene_IRI_S1\tgene_IRI_S2\tgene_IRI_difference\n")
                for i, id in enumerate(sorted(filtered_genes.keys())):
                        results_file.write(id + "\t" + str(pval_list[i]) + "\t" + str(fdr_pval_list[i]) + "\t" + ",".join(filtered_genes[id][0]) + "\t" + ",".join(filtered_genes[id][1]) + "\t" + str(IRI_diff_list[i]) + "\n")
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
        
        def run_IRC_quant_for_all_samples(self, args):
                logging.info("Perform \"IRTools quant\" for all replciates from both samples")
                self.logger.disabled = True
                
                filtered_CIR_id_list = []
                
                s1files_list = self.params['s1files'].split(',')
                s2files_list = self.params['s2files'].split(',')                

                IRC_junction_level_df_s1_list = []
                IRC_junction_level_df_s2_list = []
                
                IRC_intron_level_df_s1_list = []
                IRC_intron_level_df_s2_list = []                
                
                CIR_counts_s1_list = []
                CIR_counts_s2_list = []
                
                IRC_quanter = IRC_quant(args)
                IRC_quanter.params['outdir'] = self.temp_dir 
                
                for i, s1file in enumerate(s1files_list):
                        IRC_quanter.params['name'] = self.params['name'] + "_S1_R%d" % (i+1)
                        IRC_quanter.params['altfile'] = s1file
                        
                        IRC_quanter.quant()
                        IRC_quanter.output_IRC_junction_level()
                        IRC_junction_level_df_s1_list.append(IRC_quanter.IRC_junction_level_df)
                        
                        IRC_quanter.output_IRC_intron_level()
                        IRC_intron_level_df_s1_list.append(IRC_quanter.IRC_intron_level_df)
                        
                        CIR_counts_s1_list.append(copy.deepcopy(IRC_quanter.CIR_counts))
                        filtered_CIR_id_list.extend(IRC_quanter.filtered_CIR_id_list)
                                               
                for i, s2file in enumerate(s2files_list):
                        IRC_quanter.params['name'] = self.params['name'] + "_S2_R%d" % (i+1)
                        IRC_quanter.params['altfile'] = s2file
                        
                        IRC_quanter.quant()
                        IRC_quanter.output_IRC_junction_level()
                        IRC_junction_level_df_s2_list.append(IRC_quanter.IRC_junction_level_df)
                        
                        IRC_quanter.output_IRC_intron_level()
                        IRC_intron_level_df_s2_list.append(IRC_quanter.IRC_intron_level_df)
                        
                        CIR_counts_s2_list.append(copy.deepcopy(IRC_quanter.CIR_counts))
                        filtered_CIR_id_list.extend(IRC_quanter.filtered_CIR_id_list)
                        
                filtered_CIR_id_list = list(set(filtered_CIR_id_list))
                
                for i, IRC_intron_level_df in enumerate(IRC_intron_level_df_s1_list):
                        IRC_intron_level_df_s1_list[i] = IRC_intron_level_df[IRC_intron_level_df.CIR_id.isin(filtered_CIR_id_list) == False]
                
                for i, IRC_intron_level_df in enumerate(IRC_intron_level_df_s2_list):
                        IRC_intron_level_df_s2_list[i] = IRC_intron_level_df[IRC_intron_level_df.CIR_id.isin(filtered_CIR_id_list) == False]                     
                
                self.IRC_junction_level_data = {'s1_data': IRC_junction_level_df_s1_list,
                                                's2_data': IRC_junction_level_df_s2_list}  
                
                self.IRC_intron_level_data = {'s1_data': IRC_intron_level_df_s1_list,
                                              's2_data': IRC_intron_level_df_s2_list} 
                
                IRC_gene_level_df_s1_list = []
                IRC_gene_level_df_s2_list = []                
                
                for i, s1file in enumerate(s1files_list):
                        IRC_quanter.params['name'] = self.params['name'] + "_S1_R%d" % (i+1)
                        
                        IRC_quanter.CIR_counts = CIR_counts_s1_list[i]
                        IRC_quanter.output_IRC_gene_level(filtered_CIR_id_list)
                        IRC_gene_level_df_s1_list.append(IRC_quanter.IRC_gene_level_df)
                
                        self.logger.disabled = False
                        logging.info("Sample: " + IRC_quanter.params['name'])                   
                        IRC_quanter.output_IRC_genome_wide()
                        self.logger.disabled = True
                        
                for i, s2file in enumerate(s2files_list):
                        IRC_quanter.params['name'] = self.params['name'] + "_S2_R%d" % (i+1)
                        
                        IRC_quanter.CIR_counts = CIR_counts_s2_list[i]
                        IRC_quanter.output_IRC_gene_level(filtered_CIR_id_list)
                        IRC_gene_level_df_s2_list.append(IRC_quanter.IRC_gene_level_df)
                
                        self.logger.disabled = False
                        logging.info("Sample: " + IRC_quanter.params['name'])                   
                        IRC_quanter.output_IRC_genome_wide()
                        self.logger.disabled = True
                        
                self.IRC_gene_level_data = {'s1_data': IRC_gene_level_df_s1_list,
                                            's2_data': IRC_gene_level_df_s2_list}
                
                self.logger.disabled = False

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

                temp_dict = {}

                for i in range(len(self.params['s1files'].split(','))):
                        file_path = os.path.join(self.temp_dir, self.params['name'] + "_S1_R%d.quant.IRC.introns.txt" % (i + 1))
                        data = open(file_path, "r")
                        lines = [x.strip("\n") for x in data if x != "\n"]
                        if not temp_dict:
                                for l in lines:
                                        temp_dict[l.split()[0]] = ([l.split()[5]], [])
                        else:
                                for l in lines:
                                        temp_dict[l.split()[0]][0].append(l.split()[5])
                        data.close()

                for i in range(len(self.params['s2files'].split(','))):
                        file_path = os.path.join(self.temp_dir, self.params['name'] + "_S2_R%d.quant.IRC.introns.txt" % (i + 1))
                        data = open(file_path, "r")
                        lines = [x.strip("\n") for x in data if x != "\n"]
                        for l in lines:
                                temp_dict[l.split()[0]][1].append(l.split()[5])
                        data.close()

                del temp_dict["CIR_id"]
                self.input_intron_dict = temp_dict

                input_file_path = os.path.join(self.temp_dir, self.params['name'] + ".diff.input.IRC.introns.txt")
                input_file = open(input_file_path, "w")
                input_file.write("CIR_id\tintron_IRC_S1\tintron_IRC_S2\n")
                for id in sorted(self.input_intron_dict.keys()):
                        input_file.write(id + "\t" + ",".join(self.input_intron_dict[id][0]) + "\t" + ",".join(self.input_intron_dict[id][1]) + "\n")
                input_file.close()

                logging.info("Intron level analysis inputs can be found in " + input_file_path)

        def generate_input_gene_level(self):
                logging.info("Generating inputs for analysis for differential IR in gene level")

                temp_dict = {}

                for i in range(len(self.params['s1files'].split(','))):
                        file_path = os.path.join(self.temp_dir, self.params['name'] + "_S1_R%d.quant.IRC.genes.txt" % (i + 1))
                        data = open(file_path, "r")
                        lines = [x.strip("\n") for x in data if x != "\n"]
                        if not temp_dict:
                                for l in lines:
                                        temp_dict[l.split()[0]] = ([l.split()[4]], [])
                        else:
                                for l in lines:
                                        temp_dict[l.split()[0]][0].append(l.split()[4])
                        data.close()

                for i in range(len(self.params['s2files'].split(','))):
                        file_path = os.path.join(self.temp_dir, self.params['name'] + "_S2_R%d.quant.IRC.genes.txt" % (i + 1))
                        data = open(file_path, "r")
                        lines = [x.strip("\n") for x in data if x != "\n"]
                        for l in lines:
                                temp_dict[l.split()[0]][1].append(l.split()[4])
                        data.close()

                del temp_dict["gene_id"]
                self.input_gene_dict = temp_dict

                input_file_path = os.path.join(self.temp_dir, self.params['name'] + ".diff.input.IRC.genes.txt")
                input_file = open(input_file_path, "w")
                input_file.write("gene_id\tgene_IRC_S1\tgene_IRC_S2\n")
                for id in sorted(self.input_gene_dict.keys()):
                        input_file.write(id + "\t" + ",".join(self.input_gene_dict[id][0]) + "\t" + ",".join(self.input_gene_dict[id][1]) + "\n")
                input_file.close()

                logging.info("Gene level analysis inputs can be found at " + input_file_path)

        def generate_input_junction_level(self):
                logging.info("Generating inputs for analysis for differential IR in junction level")

                temp_dict = {}

                for i in range(len(self.params['s1files'].split(','))):
                        file_path = os.path.join(self.temp_dir, self.params['name'] + "_S1_R%d.quant.IRC.junctions.txt" % (i + 1))
                        data = open(file_path, "r")
                        lines = [x.strip("\n") for x in data if x != "\n"]
                        if not temp_dict:
                                for l in lines:
                                        temp_dict[l.split()[0]] = ([l.split()[5]], [])
                        else:
                                for l in lines:
                                        temp_dict[l.split()[0]][0].append(l.split()[5])
                        data.close()

                for i in range(len(self.params['s2files'].split(','))):
                        file_path = os.path.join(self.temp_dir, self.params['name'] + "_S2_R%d.quant.IRC.junctions.txt" % (i + 1))
                        data = open(file_path, "r")
                        lines = [x.strip("\n") for x in data if x != "\n"]
                        for l in lines:
                                temp_dict[l.split()[0]][1].append(l.split()[5])
                        data.close()

                del temp_dict["CJ_id"]
                self.input_junction_dict = temp_dict

                input_file_path = os.path.join(self.temp_dir, self.params['name'] + ".diff.input.IRC.junctions.txt")
                input_file = open(input_file_path, "w")
                input_file.write("CJ_id\tjunction_IRC_S1\tjunction_IRC_S2\n")
                for id in sorted(self.input_junction_dict.keys()):
                        input_file.write(id + "\t" + ",".join(self.input_junction_dict[id][0]) + "\t" + ",".join(self.input_junction_dict[id][1]) + "\n")
                input_file.close()

                logging.info("junction level analysis inputs can be found at " + input_file_path)

        def run_analysis_intron_level(self):
                logging.info("Running analysis for differential IR in intron level")

                filtered_introns = {}
                pval_list = []
                IRC_diff_list = []

                for id in sorted(self.input_intron_dict.keys()):
                        intron_IRC_S1 = self.input_intron_dict[id][0]
                        intron_IRC_S2 = self.input_intron_dict[id][1]
                        num_intron_IRC_S1 = []
                        num_intron_IRC_S2 = []
                        if self.params["analysistype"] == "P":
                                for i, val in enumerate(intron_IRC_S1):
                                        if val != "NA" and val != 'inf' and intron_IRC_S2[i] != "NA" and intron_IRC_S2[i] != 'inf':
                                                num_intron_IRC_S1.append(float(val))
                                                num_intron_IRC_S2.append(float(intron_IRC_S2[i]))
                        else:
                                num_intron_IRC_S1 = [float(x) for x in intron_IRC_S1 if x != "NA" and x != 'inf']
                                num_intron_IRC_S2 = [float(x) for x in intron_IRC_S2 if x != "NA" and x != 'inf']
                        if self.count_distinct_vals(num_intron_IRC_S1, num_intron_IRC_S2) == 1:
                                continue
                        if len(num_intron_IRC_S1) < 2 or len(num_intron_IRC_S2) < 2:
                                continue
                        if self.params["analysistype"] == "P":
                                pval = scipy.stats.ttest_rel(num_intron_IRC_S1, num_intron_IRC_S2)[1]
                        else:
                                pval = scipy.stats.ttest_ind(num_intron_IRC_S1, num_intron_IRC_S2)[1]
                        if math.isnan(pval):
                                continue
                        filtered_introns[id] = (intron_IRC_S1, intron_IRC_S2)
                        pval_list.append(pval)
                        diff = mean(num_intron_IRC_S2) - mean(num_intron_IRC_S1)
                        IRC_diff_list.append(diff)
                
                fdr_bool_list, fdr_pval_list = statsmodels.stats.multitest.fdrcorrection(pval_list)

                results_file_path = os.path.join(self.params['outdir'], self.params['name'] + ".diff.IRC.introns.txt")
                results_file = open(results_file_path, "w")
                results_file.write("CIR_id\tPValue\tFDR\tintron_IRC_S1\tintron_IRC_S2\tintron_IRC_difference\n")
                for i, id in enumerate(sorted(filtered_introns.keys())):
                        results_file.write(id + "\t" + str(pval_list[i]) + "\t" + str(fdr_pval_list[i]) + "\t" + ",".join(filtered_introns[id][0]) + "\t" + ",".join(filtered_introns[id][1]) + "\t" + str(IRC_diff_list[i]) + "\n")
                results_file.close()

                logging.info("Intron level differential IR results can be found in " + results_file_path)

        def run_analysis_gene_level(self):
                logging.info("Running analysis for differential IR in gene level")

                filtered_genes = {}
                pval_list = []
                IRC_diff_list = []

                for id in sorted(self.input_gene_dict.keys()):
                        gene_IRC_S1 = self.input_gene_dict[id][0]
                        gene_IRC_S2 = self.input_gene_dict[id][1]
                        num_gene_IRC_S1 = []
                        num_gene_IRC_S2 = []
                        if self.params["analysistype"] == "P":
                                for i, val in enumerate(gene_IRC_S1):
                                        if val != "NA" and val != 'inf' and gene_IRC_S2[i] != "NA" and gene_IRC_S2[i] != 'inf':
                                                num_gene_IRC_S1.append(float(val))
                                                num_gene_IRC_S2.append(float(gene_IRC_S2[i]))
                        else:
                                num_gene_IRC_S1 = [float(x) for x in gene_IRC_S1 if x != "NA" and x != 'inf']
                                num_gene_IRC_S2 = [float(x) for x in gene_IRC_S2 if x != "NA" and x != 'inf']
                        if self.count_distinct_vals(num_gene_IRC_S1, num_gene_IRC_S2) == 1:
                                continue
                        if len(num_gene_IRC_S1) < 2 or len(num_gene_IRC_S2) < 2:
                                continue
                        if self.params["analysistype"] == "P":
                                pval = scipy.stats.ttest_rel(num_gene_IRC_S1, num_gene_IRC_S2)[1]
                        else:
                                pval = scipy.stats.ttest_ind(num_gene_IRC_S1, num_gene_IRC_S2)[1]
                        if math.isnan(pval):
                                continue
                        filtered_genes[id] = (gene_IRC_S1, gene_IRC_S2)
                        pval_list.append(pval)
                        diff = mean(num_gene_IRC_S2) - mean(num_gene_IRC_S1)
                        IRC_diff_list.append(diff)

                fdr_bool_list, fdr_pval_list = statsmodels.stats.multitest.fdrcorrection(pval_list)

                results_file_path = os.path.join(self.params['outdir'], self.params['name'] + ".diff.IRC.genes.txt")
                results_file = open(results_file_path, "w")
                results_file.write("gene_id\tPValue\tFDR\tgene_IRC_S1\tgene_IRC_S2\tgene_IRC_difference\n")
                for i, id in enumerate(sorted(filtered_genes.keys())):
                        results_file.write(id + "\t" + str(pval_list[i]) + "\t" + str(fdr_pval_list[i]) + "\t" + ",".join(filtered_genes[id][0]) + "\t" + ",".join(filtered_genes[id][1]) + "\t" + str(IRC_diff_list[i]) + "\n")
                results_file.close()

                logging.info("Gene level differential IR results can be found in " + results_file_path)

        def run_analysis_junction_level(self):
                logging.info("Running analysis for differential IR in junction level")

                filtered_junctions = {}
                pval_list = []
                IRC_diff_list = []

                for id in sorted(self.input_junction_dict.keys()):
                        junction_IRC_S1 = self.input_junction_dict[id][0]
                        junction_IRC_S2 = self.input_junction_dict[id][1]
                        num_junction_IRC_S1 = []
                        num_junction_IRC_S2 = []
                        if self.params["analysistype"] == "P":
                                for i, val in enumerate(junction_IRC_S1):
                                        if val != "NA" and val != 'inf' and junction_IRC_S2[i] != "NA" and junction_IRC_S2[i] != 'inf':
                                                num_junction_IRC_S1.append(float(val))
                                                num_junction_IRC_S2.append(float(junction_IRC_S2[i]))
                        else:
                                num_junction_IRC_S1 = [float(x) for x in junction_IRC_S1 if x != "NA" and x != 'inf']
                                num_junction_IRC_S2 = [float(x) for x in junction_IRC_S2 if x != "NA" and x != 'inf']
                        if self.count_distinct_vals(num_junction_IRC_S1, num_junction_IRC_S2) == 1:
                                continue
                        if len(num_junction_IRC_S1) < 2 or len(num_junction_IRC_S2) < 2:
                                continue
                        if self.params["analysistype"] == "P":
                                pval = scipy.stats.ttest_rel(num_junction_IRC_S1, num_junction_IRC_S2)[1]
                        else:
                                pval = scipy.stats.ttest_ind(num_junction_IRC_S1, num_junction_IRC_S2)[1]
                        if math.isnan(pval):
                                continue
                        filtered_junctions[id] = (junction_IRC_S1, junction_IRC_S2)
                        pval_list.append(pval)
                        diff = mean(num_junction_IRC_S2) - mean(num_junction_IRC_S1)
                        IRC_diff_list.append(diff)

                fdr_bool_list, fdr_pval_list = statsmodels.stats.multitest.fdrcorrection(pval_list)

                results_file_path = os.path.join(self.params['outdir'], self.params['name'] + ".diff.IRC.junctions.txt")
                results_file = open(results_file_path, "w")
                results_file.write("CJ_id\tPValue\tFDR\tjunction_IRC_S1\tjunction_IRC_S2\tjunction_IRC_difference\n")
                for i, id in enumerate(sorted(filtered_junctions.keys())):
                        results_file.write(id + "\t" + str(pval_list[i]) + "\t" + str(fdr_pval_list[i]) + "\t" + ",".join(filtered_junctions[id][0]) + "\t" + ",".join(filtered_junctions[id][1]) + "\t" + str(IRC_diff_list[i]) + "\n")
                results_file.close()

                logging.info("Junction level differential IR results can be found in " + results_file_path)

def run(args):
        if args.quanttype == "IRI":
                IRI_differ = IRI_diff(args)
                if ',' not in IRI_differ.params['s1files'] or ',' not in IRI_differ.params['s2files']:
                        logging.info("Run Aborted: Differential IR analysis requires at least two replicates per sample. Please check input.")
                        exit()
                if IRI_differ.params['analysistype'] == "P" and IRI_differ.params['s1files'].count(',') != IRI_differ.params['s2files'].count(','):
                        logging.info("Run Aborted: Samples must have the same number of replicates for paired analysis. Please check input.")
                        exit()
                IRI_differ.run_IRI_quant_for_all_samples(args)
                IRI_differ.generate_input_intron_level()
                IRI_differ.generate_input_gene_level()
                IRI_differ.run_analysis_intron_level()
                IRI_differ.run_analysis_gene_level()

        elif args.quanttype == "IRC":
                IRC_differ = IRC_diff(args)
                if ',' not in IRC_differ.params['s1files'] or ',' not in IRC_differ.params['s2files']:
                        logging.info("Run Aborted: Differential IR analysis requires at least two replicates per sample. Please check input.")
                        exit()
                if IRC_differ.params['analysistype'] == "P" and IRC_differ.params['s1files'].count(',') != IRC_differ.params['s2files'].count(','):
                        logging.info("Run Aborted: Samples must have the same number of replicates for paired analysis. Please check input.")
                        exit()
                IRC_differ.run_IRC_quant_for_all_samples(args)
                IRI_differ.generate_input_intron_level()
                IRI_differ.generate_input_gene_level()
                IRI_differ.generate_input_junction_level()
                IRI_differ.run_analysis_intron_level()
                IRI_differ.run_analysis_gene_level()
                IRI_differ.run_analysis_junction_level()



