import logging
import os

from IRTools.quant_RT import RT_quant


class ONT_RT_quant(RT_quant):
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

        def quant(self):
                logging.info("Accumulating aligned ONT bases in exonic and read-through regions")
                self.total_read_count = 0
                self.assigned_read_count = 0
                self.total_aligned_bases = 0
                self.aligned_base_counts = {
                        gene_id: {
                                "constitutive_exonic_region": 0,
                                "readthrough_region": 0,
                        }
                        for gene_id in self.valid_genes
                }

                for alt in self.iter_single_reads():
                        if not self.usable_alignment(alt):
                                continue

                        alt_iv_seq = self.get_alt_iv(alt)
                        if not alt_iv_seq:
                                continue

                        self.total_read_count += 1
                        self.total_aligned_bases += sum(iv.length for iv in alt_iv_seq)
                        read_assignments = self.collect_read_assignments(alt_iv_seq)
                        if read_assignments:
                                self.assigned_read_count += 1
                                _, region_hits = read_assignments
                                for iv, feature_assignment in region_hits:
                                        gene_id, feature_type = feature_assignment
                                        self.aligned_base_counts[gene_id][feature_type] += iv.length

                if self.total_aligned_bases == 0:
                        return

                mean_aligned_length = self.total_aligned_bases * 1.0 / self.total_read_count
                for gene_id, feature_counts in self.aligned_base_counts.items():
                        for feature_type, aligned_bases in feature_counts.items():
                                self.counts[gene_id][feature_type] = aligned_bases * 1.0 / mean_aligned_length

        def output_RT_gene_level(self):
                super().output_RT_gene_level()
                self.RT_gene_level_df["gene_CER_aligned_bases"] = [
                        self.aligned_base_counts[gene_id]["constitutive_exonic_region"]
                        for gene_id in self.RT_gene_level_df["gene_id"]
                ]
                self.RT_gene_level_df["RT_aligned_bases"] = [
                        self.aligned_base_counts[gene_id]["readthrough_region"]
                        for gene_id in self.RT_gene_level_df["gene_id"]
                ]
                outfile_fullpath = os.path.join(self.params["outdir"], self.RT_gene_output_filename())
                self.RT_gene_level_df.to_csv(outfile_fullpath, index=None, sep="\t", na_rep="NA")

        def output_RT_genome_wide(self):
                stats = super().output_RT_genome_wide()
                stats["Assigned alignment count"] = self.assigned_read_count
                stats["Total aligned bases"] = self.total_aligned_bases
                return stats

        def RT_gene_output_filename(self):
                return self.params["name"] + ".ont.RT.genes.txt"
