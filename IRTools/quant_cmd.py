import csv
import glob
import logging
import os
from copy import copy


def matched_altfiles(altfile_pattern):
        matches = sorted(glob.glob(altfile_pattern))
        if matches:
                return matches
        return [altfile_pattern]


def sample_name_from_altfile(altfile):
        return os.path.splitext(os.path.basename(altfile))[0]


def open_unique_summary_file(outdir, basename):
        stem, ext = os.path.splitext(basename)
        candidate_names = [basename]
        candidate_names.extend("{}_{}{}".format(stem, i, ext) for i in range(1, 1000000))

        for candidate_name in candidate_names:
                outfile = os.path.join(outdir, candidate_name)
                try:
                        handle = open(outfile, "x", newline="", encoding="utf-8")
                        return outfile, handle
                except FileExistsError:
                        continue

        raise RuntimeError("Could not create a unique summary file for {}".format(basename))


def write_iri_summary(outdir, summary_rows):
        outfile, handle = open_unique_summary_file(outdir, "CIR_ratio.txt")
        with handle:
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow(["sample_name", "Total CIR", "Total CER", "CIR ratio", "CER ratio", "Eligible read count"])
                for row in summary_rows:
                        writer.writerow([row["sample_name"], row["Total CIR"], row["Total CER"], row["CIR ratio"], row["CER ratio"], row["Eligible read count"]])
        logging.info("Wrote IRI summary file: %s", outfile)


def write_rt_summary(outdir, summary_rows):
        outfile, handle = open_unique_summary_file(outdir, "RT_ratio.txt")
        with handle:
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow(["sample_name", "Total RT", "Total CER", "RT ratio", "CER ratio", "Eligible read count"])
                for row in summary_rows:
                        writer.writerow([row["sample_name"], row["Total RT"], row["Total CER"], row["RT ratio"], row["CER ratio"], row["Eligible read count"]])
        logging.info("Wrote RT summary file: %s", outfile)


def write_irc_summary(outdir, summary_rows):
        outfile, handle = open_unique_summary_file(outdir, "IRC_ratio.txt")
        with handle:
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow(["sample_name", "Total retained reads", "Total spliced reads", "Retained ratio", "Spliced ratio", "Total junction reads", "Genome-wide IRC"])
                for row in summary_rows:
                        writer.writerow([row["sample_name"], row["Total retained reads"], row["Total spliced reads"], row["Retained ratio"], row["Spliced ratio"], row["Total junction reads"], row["Genome-wide IRC"]])
        logging.info("Wrote IRC summary file: %s", outfile)


def run_iri_quant(args):
        from IRTools.quant_IRI import IRI_quant
        summary_rows = []
        altfiles = matched_altfiles(args.altfile)
        multifile_mode = len(altfiles) > 1 or any(token in args.altfile for token in ["*", "?", "["])

        for altfile in altfiles:
                run_args = copy(args)
                run_args.altfile = altfile
                if multifile_mode:
                        run_args.name = sample_name_from_altfile(altfile)

                IRI_quanter = IRI_quant(run_args)
                IRI_quanter.quant()
                IRI_quanter.output_IRI_intron_level()
                IRI_quanter.output_IRI_gene_level()
                summary_rows.append(IRI_quanter.output_IRI_genome_wide())

        write_iri_summary(args.outdir, summary_rows)


def run_irc_quant(args):
        from IRTools.quant_IRC import IRC_quant
        summary_rows = []
        altfiles = matched_altfiles(args.altfile)
        multifile_mode = len(altfiles) > 1 or any(token in args.altfile for token in ["*", "?", "["])

        for altfile in altfiles:
                run_args = copy(args)
                run_args.altfile = altfile
                if multifile_mode:
                        run_args.name = sample_name_from_altfile(altfile)

                IRC_quanter = IRC_quant(run_args)
                IRC_quanter.quant()
                IRC_quanter.output_IRC_junction_level()
                IRC_quanter.output_IRC_intron_level()
                IRC_quanter.output_IRC_gene_level()
                summary_rows.append(IRC_quanter.output_IRC_genome_wide())

        write_irc_summary(args.outdir, summary_rows)


def run_rt_quant(args):
        from IRTools.quant_RT import RT_quant
        summary_rows = []
        altfiles = matched_altfiles(args.altfile)
        multifile_mode = len(altfiles) > 1 or any(token in args.altfile for token in ["*", "?", "["])

        for altfile in altfiles:
                run_args = copy(args)
                run_args.altfile = altfile
                if multifile_mode:
                        run_args.name = sample_name_from_altfile(altfile)

                RT_quanter = RT_quant(run_args)
                RT_quanter.quant()
                RT_quanter.output_RT_gene_level()
                summary_rows.append(RT_quanter.output_RT_genome_wide())

        write_rt_summary(args.outdir, summary_rows)

def run(args):
        if args.quanttype == 'IRI':
                run_iri_quant(args)

        elif args.quanttype == 'IRC':
                run_irc_quant(args)

        elif args.quanttype == 'RT':
                run_rt_quant(args)
                
                




