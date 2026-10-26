import csv
import os
import tempfile
import unittest
from types import SimpleNamespace

from IRTools.diff_cmd import generate_diff_input_dict, run


class DiffInputFormatTests(unittest.TestCase):
        def setUp(self):
                self.temp_dir = tempfile.TemporaryDirectory()

        def tearDown(self):
                self.temp_dir.cleanup()

        def write_result(self, sample_name, suffix, headers, values):
                file_path = os.path.join(self.temp_dir.name, sample_name + suffix)
                with open(file_path, "w", newline="") as output:
                        writer = csv.writer(output, delimiter="\t")
                        writer.writerow(headers)
                        writer.writerows(values)

        def params(self, input_format="auto"):
                return {
                        "indir": self.temp_dir.name,
                        "s1files": "control",
                        "s2files": "treated",
                        "analysistype": "U",
                        "inputformat": input_format,
                }

        def test_auto_detects_ont_iri_and_reads_named_metric_columns(self):
                headers = ["CIR_id", "CIR_iv", "intron_IRI", "CIR_RPKM"]
                self.write_result("control", ".ont.IRI.introns.txt", headers, [["gene1:001", "chr1:1-10", "0.2", "5"]])
                self.write_result("treated", ".ont.IRI.introns.txt", headers, [["gene1:001", "chr1:1-10", "0.6", "8"]])

                values = generate_diff_input_dict(
                        self.params(), "IRI", "introns", "intron_IRI", "CIR_id"
                )

                self.assertEqual(values["gene1:001"], (["0.2"], ["0.6"]))

        def test_auto_detects_ont_rt(self):
                headers = ["gene_id", "gene_RT_index", "gene_CER_aligned_bases", "RT_aligned_bases"]
                self.write_result("control", ".ont.RT.genes.txt", headers, [["gene1", "0.1", "1000", "100"]])
                self.write_result("treated", ".ont.RT.genes.txt", headers, [["gene1", "0.4", "1000", "400"]])

                values = generate_diff_input_dict(
                        self.params(), "RT", "genes", "gene_RT_index", "gene_id"
                )

                self.assertEqual(values["gene1"], (["0.1"], ["0.4"]))

        def test_existing_ngs_format_still_works(self):
                headers = ["gene_id", "gene_iv", "gene_IRI"]
                self.write_result("control", ".quant.IRI.genes.txt", headers, [["gene1", "chr1:1-10", "0.3"]])
                self.write_result("treated", ".quant.IRI.genes.txt", headers, [["gene1", "chr1:1-10", "0.5"]])

                values = generate_diff_input_dict(
                        self.params(), "IRI", "genes", "gene_IRI", "gene_id"
                )

                self.assertEqual(values["gene1"], (["0.3"], ["0.5"]))

        def test_auto_rejects_mixed_ngs_and_ont_samples(self):
                headers = ["gene_id", "gene_RT_index"]
                self.write_result("control", ".quant.RT.genes.txt", headers, [["gene1", "0.1"]])
                self.write_result("treated", ".ont.RT.genes.txt", headers, [["gene1", "0.4"]])

                with self.assertRaisesRegex(ValueError, "cannot be mixed"):
                        generate_diff_input_dict(
                                self.params(), "RT", "genes", "gene_RT_index", "gene_id"
                        )

        def test_explicit_format_resolves_ambiguous_sample_names(self):
                headers = ["gene_id", "gene_RT_index"]
                for sample_name in ("control", "treated"):
                        self.write_result(sample_name, ".quant.RT.genes.txt", headers, [["gene1", "9.9"]])
                        self.write_result(sample_name, ".ont.RT.genes.txt", headers, [["gene1", "0.2"]])

                values = generate_diff_input_dict(
                        self.params("ONT"), "RT", "genes", "gene_RT_index", "gene_id"
                )

                self.assertEqual(values["gene1"], (["0.2"], ["0.2"]))

        def test_ont_iri_runs_end_to_end(self):
                intron_headers = ["CIR_id", "intron_IRI"]
                gene_headers = ["gene_id", "gene_IRI"]
                for sample_name, intron_iri, gene_iri in (
                        ("control", "0.2", "0.3"),
                        ("treated", "0.6", "0.8"),
                ):
                        self.write_result(
                                sample_name,
                                ".ont.IRI.introns.txt",
                                intron_headers,
                                [["gene1:001", intron_iri]],
                        )
                        self.write_result(
                                sample_name,
                                ".ont.IRI.genes.txt",
                                gene_headers,
                                [["gene1", gene_iri]],
                        )

                args = SimpleNamespace(
                        quanttype="IRI",
                        indir=self.temp_dir.name,
                        inputformat="ONT",
                        s1files="control",
                        s2files="treated",
                        analysistype="U",
                        method="minus",
                        name="control_vs_treated",
                        outdir=self.temp_dir.name,
                )
                run(args)

                output_path = os.path.join(
                        self.temp_dir.name,
                        "control_vs_treated.diff.IRI.introns.txt",
                )
                with open(output_path, "r", newline="") as output:
                        result = next(csv.DictReader(output, delimiter="\t"))

                self.assertEqual(result["PValue"], "NA")
                self.assertEqual(result["FDR"], "NA")
                self.assertAlmostEqual(float(result["intron_IRI_difference"]), 0.4)


if __name__ == "__main__":
        unittest.main()
