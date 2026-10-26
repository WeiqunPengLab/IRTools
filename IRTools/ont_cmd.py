import csv
import logging
from copy import copy

from IRTools.quant_cmd import matched_altfiles, open_unique_summary_file, sample_name_from_altfile


def write_ont_iri_summary(outdir, summary_rows):
        outfile, handle = open_unique_summary_file(outdir, "ONT_IRI_ratio.txt")
        with handle:
                writer = csv.writer(handle, delimiter="\t")
                columns = [
                        "sample_name", "Total CIR", "Total CER", "CIR ratio", "CER ratio",
                        "Eligible alignment count", "Assigned alignment count", "Total aligned bases",
                        "Genome-wide ONT IRI", "Total retained reads", "Total spliced reads",
                        "Retained ratio", "Spliced ratio", "Total IR reads",
                        "Genome-wide ONT IR event ratio",
                ]
                writer.writerow(columns)
                for row in summary_rows:
                        writer.writerow([row[column] for column in columns])
        logging.info("Wrote ONT IRI summary file: %s", outfile)


def write_ont_rt_summary(outdir, summary_rows):
        outfile, handle = open_unique_summary_file(outdir, "ONT_RT_ratio.txt")
        with handle:
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow(["sample_name", "Total RT", "Total CER", "RT ratio", "CER ratio", "Eligible alignment count", "Assigned alignment count", "Total aligned bases", "Genome-wide ONT RT"])
                for row in summary_rows:
                        writer.writerow([row["sample_name"], row["Total RT"], row["Total CER"], row["RT ratio"], row["CER ratio"], row["Eligible read count"], row["Assigned alignment count"], row["Total aligned bases"], row["Genome-wide RT"]])
        logging.info("Wrote ONT RT summary file: %s", outfile)


def run_iri(args):
        from IRTools.ont_IRI import ONT_IRI_quant

        summary_rows = []
        altfiles = matched_altfiles(args.altfile)
        multifile_mode = len(altfiles) > 1 or any(token in args.altfile for token in ["*", "?", "["])

        for altfile in altfiles:
                run_args = copy(args)
                run_args.altfile = altfile
                if multifile_mode:
                        run_args.name = sample_name_from_altfile(altfile)

                ONT_quanter = ONT_IRI_quant(run_args)
                ONT_quanter.quant()
                ONT_quanter.output_ONT_intron_level()
                ONT_quanter.output_ONT_gene_level()
                summary_rows.append(ONT_quanter.output_ONT_genome_wide())

        write_ont_iri_summary(args.outdir, summary_rows)


def run_rt(args):
        from IRTools.ont_RT import ONT_RT_quant

        summary_rows = []
        altfiles = matched_altfiles(args.altfile)
        multifile_mode = len(altfiles) > 1 or any(token in args.altfile for token in ["*", "?", "["])

        for altfile in altfiles:
                run_args = copy(args)
                run_args.altfile = altfile
                if multifile_mode:
                        run_args.name = sample_name_from_altfile(altfile)

                ONT_quanter = ONT_RT_quant(run_args)
                ONT_quanter.quant()
                ONT_quanter.output_RT_gene_level()
                summary_rows.append(ONT_quanter.output_RT_genome_wide())

        write_ont_rt_summary(args.outdir, summary_rows)


def run(args):
        if args.quanttype == "IRI":
                run_iri(args)
        elif args.quanttype == "RT":
                run_rt(args)
