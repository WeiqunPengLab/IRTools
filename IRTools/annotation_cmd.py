import sys
import os
import logging
import re
import collections
import copy
import gzip
import hashlib
import shutil
import tempfile
import time
import urllib.request, urllib.parse, urllib.error
import subprocess
import HTSeq

INPUT_FEATURE_TYPES = {"exon"}
TRANSCRIPT_FEATURE_TYPES = {"transcript", "mRNA", "lnc_RNA", "ncRNA", "rRNA", "tRNA", "miRNA", "snRNA", "snoRNA"}
MAPFILE_SPECIES_ALIASES = {"mm11": "mm39"}
SPECIES_MAPPABILITY_SOURCES = {
        # hg19/hg38/mm9/mm10 follow the Hoffman Lab raw layout visible on the site.
        "hg19": ("bedgraph", "https://bismap.hoffmanlab.org/raw/hg19/k{read_length}.umap.bedgraph.gz"),
        "hg38": ("bedgraph", "https://bismap.hoffmanlab.org/raw/hg38/k{read_length}.umap.bedgraph.gz"),
        "mm9": ("bedgraph", "https://bismap.hoffmanlab.org/raw/mm9/k{read_length}.umap.bedgraph.gz"),
        "mm10": ("bedgraph", "https://bismap.hoffmanlab.org/raw/mm10/k{read_length}.umap.bedgraph.gz"),
        "mm39": ("bigwig", "https://bismap.hoffmanlab.org/raw/mm39/umap/bigWigFiles/k{read_length}.merged.GRCm39_umap.MultiReadMappability.bigWig"),
}

def normalize_parent_ids(value):
        if value is None:
                return []
        if isinstance(value, list):
                parent_values = value
        else:
                parent_values = [value]
        normalized = []
        for item in parent_values:
                normalized.extend([part for part in str(item).split(",") if part])
        return normalized

def first_attr(feature, keys):
        for key in keys:
                value = feature.attr.get(key)
                if value:
                        if isinstance(value, list):
                                if len(value) > 0:
                                        return value[0]
                        else:
                                return value
        return None

def extract_refseq_accession_from_text(text):
        if not text:
                return None
        match = re.search(r'([NX][MR]_\d+(?:\.\d+)?)', str(text))
        if match:
                return match.group(1)
        return None

def normalize_map_species(species_name):
        return MAPFILE_SPECIES_ALIASES.get(species_name, species_name)

def default_mappability_cache_dir():
        return os.path.join(os.path.dirname(__file__), "data", "mappability_cache")

def mappability_cache_dir(cache_dir=None):
        if cache_dir:
                return os.path.abspath(os.path.expanduser(cache_dir))
        return default_mappability_cache_dir()

def cache_dir_is_writable(cache_dir):
        try:
                os.makedirs(cache_dir, exist_ok=True)
                with tempfile.NamedTemporaryFile(prefix=".irtools_write_test_", dir=cache_dir, delete=True):
                        pass
                return True
        except OSError:
                return False

def require_mappability_cache_dir(cache_dir):
        if cache_dir_is_writable(cache_dir):
                return

        default_cache_dir = default_mappability_cache_dir()
        if os.path.abspath(cache_dir) == os.path.abspath(default_cache_dir):
                raise Exception(
                        "IRTools cannot write the mappability cache inside the installed package:\n"
                        "  {}\n"
                        "Please rerun annotation with a writable cache directory, for example:\n"
                        "  IRTools annotation ... --map-cache-dir /path/to/project/irtools_mappability_cache".format(cache_dir)
                )

        raise Exception(
                "IRTools cannot write the mappability cache directory passed to --map-cache-dir:\n"
                "  {}\n"
                "Please choose a writable directory and rerun annotation.".format(cache_dir)
        )

def annotation_cache_token(inputfile):
        inputfile = os.path.abspath(inputfile)
        stat = os.stat(inputfile)
        token_input = "{}|{}|{}".format(inputfile, stat.st_size, stat.st_mtime_ns)
        return hashlib.md5(token_input.encode("utf-8")).hexdigest()[:12]

def cached_blacklist_path(cache_stem, cache_dir=None):
        safe_stem = re.sub(r'[^A-Za-z0-9._-]+', "_", cache_stem)
        return os.path.join(mappability_cache_dir(cache_dir), safe_stem + ".bed.gz")

def report_cache_status(cache_path, exists):
        if exists:
                logging.info("Using cached low-mappability BED: %s", cache_path)
                print("Using cached low-mappability BED: {}".format(cache_path))
        else:
                logging.info("Low-mappability BED cache is missing: %s", cache_path)
                print("Low-mappability BED cache not found: {}".format(cache_path))
                print("IRTools needs to generate it now. This may take some time on the first run.")

def validate_download(temp_path, suffix):
        if suffix.endswith(".gz"):
                with gzip.open(temp_path, "rt", encoding="utf-8") as handle:
                        handle.read(1)

def download_with_wget(url, suffix, retries):
        wget_path = shutil.which("wget")
        if not wget_path:
                return None

        temp_dir = tempfile.mkdtemp(prefix="irtools_map_")
        temp_path = os.path.join(temp_dir, "download" + suffix)
        for attempt in range(1, retries + 1):
                logging.info("Downloading mappability resource with wget from %s (attempt %i/%i)", url, attempt, retries)
                print("Downloading mappability resource with wget from {} (attempt {}/{})".format(url, attempt, retries))
                try:
                        subprocess.run(
                                [wget_path, "-c", "-O", temp_path, url],
                                check=True,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                text=True,
                        )
                        validate_download(temp_path, suffix)
                        return temp_path
                except Exception as exc:
                        if attempt < retries:
                                wait_seconds = min(5 * attempt, 15)
                                logging.warning("wget download attempt %i failed for %s: %s. Retrying in %i seconds.", attempt, url, exc, wait_seconds)
                                print("wget download attempt {} failed. Retrying in {} seconds.".format(attempt, wait_seconds))
                                time.sleep(wait_seconds)
                        else:
                                if os.path.exists(temp_path):
                                        os.unlink(temp_path)
                                try:
                                        os.rmdir(temp_dir)
                                except OSError:
                                        pass
                                raise
        return None

def download_to_tempfile(url, suffix, retries=3):
        wget_download = download_with_wget(url, suffix, retries)
        if wget_download:
                return wget_download

        last_error = None
        for attempt in range(1, retries + 1):
                temp_handle = tempfile.NamedTemporaryFile(delete=False, suffix=suffix)
                temp_path = temp_handle.name
                temp_handle.close()
                logging.info("Downloading mappability resource from %s (attempt %i/%i)", url, attempt, retries)
                print("Downloading mappability resource from {} (attempt {}/{})".format(url, attempt, retries))
                try:
                        with urllib.request.urlopen(url, timeout=120) as response, open(temp_path, "wb") as out_handle:
                                expected_size = response.headers.get("Content-Length")
                                while True:
                                        chunk = response.read(10 * 1024 * 1024)
                                        if not chunk:
                                                break
                                        out_handle.write(chunk)

                        if expected_size is not None and os.path.getsize(temp_path) != int(expected_size):
                                raise IOError(
                                        "Downloaded file size {} does not match expected size {}.".format(
                                                os.path.getsize(temp_path),
                                                expected_size,
                                        )
                                )

                        validate_download(temp_path, suffix)
                        return temp_path
                except Exception as exc:
                        last_error = exc
                        if os.path.exists(temp_path):
                                os.unlink(temp_path)
                        if attempt < retries:
                                wait_seconds = min(5 * attempt, 15)
                                logging.warning("Download attempt %i failed for %s: %s. Retrying in %i seconds.", attempt, url, exc, wait_seconds)
                                print("Download attempt {} failed. Retrying in {} seconds.".format(attempt, wait_seconds))
                                time.sleep(wait_seconds)

        raise IOError("Failed to download a valid mappability resource from {} after {} attempts: {}".format(url, retries, last_error))

def iter_low_mappability_segments(feature, entries, map_score_cutoff):
        current = feature.iv.start
        for entry_start, entry_end, score in sorted(entries, key=lambda item: (item[0], item[1])):
                start = max(feature.iv.start, int(entry_start))
                end = min(feature.iv.end, int(entry_end))
                if end <= start:
                        continue

                if start > current:
                        yield current, start

                if score < map_score_cutoff:
                        yield start, end

                if end > current:
                        current = end

        if current < feature.iv.end:
                yield current, feature.iv.end

def build_blacklist_from_bigwig(bigwig_path, blacklist_path, features, map_score_cutoff):
        from bx.bbi.bigwig_file import BigWigFile

        os.makedirs(os.path.dirname(blacklist_path), exist_ok=True)
        with open(bigwig_path, "rb") as bigwig_handle:
                bigwig = BigWigFile(bigwig_handle)
                with gzip.open(blacklist_path, "wt", encoding="utf-8") as out_handle:
                        for feature in features:
                                chrom = feature.iv.chrom.encode("utf-8")
                                entries = bigwig.get(chrom, feature.iv.start, feature.iv.end) or []
                                for start, end in iter_low_mappability_segments(feature, entries, map_score_cutoff):
                                        out_handle.write(
                                                "{}\t{}\t{}\t.\t0\t{}\n".format(
                                                        feature.iv.chrom,
                                                        start,
                                                        end,
                                                        feature.iv.strand,
                                                )
                                        )

def build_blacklist_from_bedgraph(bedgraph_path, blacklist_path, features, map_score_cutoff):
        feature_entries = collections.defaultdict(list)
        feature_state = {}
        for feature in features:
                feature_entries[(feature.iv.chrom, feature.iv.strand)].append(feature)

        for key in feature_entries:
                feature_entries[key].sort(key=lambda feature: feature.iv.start)
                feature_state[key] = {"index": 0, "active": []}

        opener = gzip.open if bedgraph_path.endswith(".gz") else open
        with opener(bedgraph_path, "rt", encoding="utf-8") as in_handle:
                for line in in_handle:
                        if not line.strip() or line.startswith("#") or line.startswith("track") or line.startswith("browser"):
                                continue
                        fields = line.rstrip("\n").split("\t")
                        if len(fields) < 4:
                                continue
                        chrom = fields[0]
                        start = int(fields[1])
                        end = int(fields[2])
                        score = float(fields[3])
                        for strand in ("+", "-"):
                                key = (chrom, strand)
                                chrom_features = feature_entries.get(key, [])
                                if not chrom_features:
                                        continue
                                state = feature_state[key]
                                while state["index"] < len(chrom_features) and chrom_features[state["index"]].iv.start < end:
                                        state["active"].append(chrom_features[state["index"]])
                                        state["index"] += 1
                                state["active"] = [feature for feature in state["active"] if feature.iv.end > start]
                                for feature in state["active"]:
                                        feature.attr.setdefault("_map_entries", []).append((start, end, score))

        os.makedirs(os.path.dirname(blacklist_path), exist_ok=True)
        with gzip.open(blacklist_path, "wt", encoding="utf-8") as out_handle:
                for feature in features:
                        entries = feature.attr.pop("_map_entries", [])
                        for start, end in iter_low_mappability_segments(feature, entries, map_score_cutoff):
                                out_handle.write(
                                        "{}\t{}\t{}\t.\t0\t{}\n".format(
                                                feature.iv.chrom,
                                                start,
                                                end,
                                                feature.iv.strand,
                                        )
                                )

def resolve_blacklist_bed(mapfile, inputfile, features, read_length, map_score_cutoff, cache_dir=None):
        if not mapfile:
                return None

        normalized_species = normalize_map_species(mapfile)
        source_token = annotation_cache_token(inputfile)
        resolved_cache_dir = mappability_cache_dir(cache_dir)

        if normalized_species in SPECIES_MAPPABILITY_SOURCES:
                cache_stem = "{}_k{}_cutoff{}_{}".format(normalized_species, read_length, str(map_score_cutoff).replace(".", "p"), source_token)
                blacklist_path = cached_blacklist_path(cache_stem, resolved_cache_dir)
                report_cache_status(blacklist_path, os.path.exists(blacklist_path))
                if os.path.exists(blacklist_path):
                        return blacklist_path
                require_mappability_cache_dir(resolved_cache_dir)

                source_type, source_url_template = SPECIES_MAPPABILITY_SOURCES[normalized_species]
                source_url = source_url_template.format(read_length=read_length)
                temp_suffix = ".bw" if source_type == "bigwig" else ".bedgraph.gz"
                temp_resource = download_to_tempfile(source_url, temp_suffix)
                try:
                        if source_type == "bigwig":
                                build_blacklist_from_bigwig(temp_resource, blacklist_path, features, map_score_cutoff)
                        else:
                                build_blacklist_from_bedgraph(temp_resource, blacklist_path, features, map_score_cutoff)
                finally:
                        if os.path.exists(temp_resource):
                                os.unlink(temp_resource)
                logging.info("Generated low-mappability BED cache: %s", blacklist_path)
                print("Generated low-mappability BED cache: {}".format(blacklist_path))
                return blacklist_path

        if re.search(r'\.bed(\.gz)?$', mapfile, re.IGNORECASE):
                return mapfile

        if re.search(r'\.bedgraph(\.gz)?$', mapfile, re.IGNORECASE):
                bedgraph_path = os.path.abspath(mapfile)
                stat = os.stat(bedgraph_path)
                map_token_input = "{}|{}|{}".format(bedgraph_path, stat.st_size, stat.st_mtime_ns)
                map_token = hashlib.md5(map_token_input.encode("utf-8")).hexdigest()[:12]
                cache_stem = "{}_k{}_cutoff{}_{}_{}".format(
                        os.path.splitext(os.path.basename(bedgraph_path))[0],
                        read_length,
                        str(map_score_cutoff).replace(".", "p"),
                        source_token,
                        map_token,
                )
                blacklist_path = cached_blacklist_path(cache_stem, resolved_cache_dir)
                report_cache_status(blacklist_path, os.path.exists(blacklist_path))
                if os.path.exists(blacklist_path):
                        return blacklist_path
                require_mappability_cache_dir(resolved_cache_dir)
                build_blacklist_from_bedgraph(bedgraph_path, blacklist_path, features, map_score_cutoff)
                logging.info("Generated low-mappability BED cache: %s", blacklist_path)
                print("Generated low-mappability BED cache: {}".format(blacklist_path))
                return blacklist_path

        if re.search(r'(bigWig|\.bw)$', mapfile, re.IGNORECASE):
                bigwig_path = os.path.abspath(mapfile)
                stat = os.stat(bigwig_path)
                map_token_input = "{}|{}|{}".format(bigwig_path, stat.st_size, stat.st_mtime_ns)
                map_token = hashlib.md5(map_token_input.encode("utf-8")).hexdigest()[:12]
                cache_stem = "{}_k{}_cutoff{}_{}_{}".format(
                        os.path.splitext(os.path.basename(bigwig_path))[0],
                        read_length,
                        str(map_score_cutoff).replace(".", "p"),
                        source_token,
                        map_token,
                )
                blacklist_path = cached_blacklist_path(cache_stem, resolved_cache_dir)
                report_cache_status(blacklist_path, os.path.exists(blacklist_path))
                if os.path.exists(blacklist_path):
                        return blacklist_path
                require_mappability_cache_dir(resolved_cache_dir)
                build_blacklist_from_bigwig(bigwig_path, blacklist_path, features, map_score_cutoff)
                logging.info("Generated low-mappability BED cache: %s", blacklist_path)
                print("Generated low-mappability BED cache: {}".format(blacklist_path))
                return blacklist_path

        raise Exception(
                "\"{}\" is neither a supported species shortcut, nor a local BED file, nor a local bedGraph file, nor a local bigWig file.".format(mapfile)
        )

def effective_lengths_from_blacklist(features, blacklist_bed):
        gene_map_score = HTSeq.GenomicArray("auto", stranded=True, typecode="i")
        for feature in features:
                gene_map_score[feature.iv] = 1

        opener = gzip.open if blacklist_bed.endswith(".gz") else open
        with opener(blacklist_bed, "rt", encoding="utf-8") as bed_handle:
                for line in bed_handle:
                        if not line.strip() or line.startswith("#"):
                                continue
                        fields = line.rstrip("\n").split("\t")
                        chrom = fields[0]
                        start = int(fields[1])
                        end = int(fields[2])
                        strand = fields[5] if len(fields) >= 6 else "."
                        gene_map_score[HTSeq.GenomicInterval(chrom, start, end, strand)] = 0

        effective_lengths = {}
        for feature in features:
            effective_lengths[(feature.attr["gene_id"], feature.attr["constitutive_intronic_region_number"])] = int(sum(gene_map_score[feature.iv]))
        return effective_lengths

def transcript_accession(feature):
        attr_candidates = []
        for key in ["transcript_id", "transcript_name", "Name", "ID", "Dbxref", "transcript"]:
                value = feature.attr.get(key)
                if isinstance(value, list):
                        attr_candidates.extend(value)
                elif value:
                        attr_candidates.append(value)
        for candidate in attr_candidates:
                accession = extract_refseq_accession_from_text(candidate)
                if accession:
                        return accession
        return None

def normalize_input_features(raw_features):
        transcript_info = {}

        for feature in raw_features:
                if feature.type not in TRANSCRIPT_FEATURE_TYPES:
                        continue
                raw_id = first_attr(feature, ["ID", "transcript_id", "transcript_name", "Name"])
                parent_gene_id = first_attr(feature, ["gene_id", "gene", "gene_name"])
                if not parent_gene_id:
                        parent_ids = normalize_parent_ids(feature.attr.get("Parent"))
                        if parent_ids:
                                parent_gene_id = parent_ids[0]
                transcript_id = first_attr(feature, ["transcript_id", "transcript_name", "Name", "ID"]) or raw_id
                accession = transcript_accession(feature)
                info = {
                        "gene_id": parent_gene_id,
                        "transcript_id": transcript_id,
                        "accession": accession,
                }
                for key in {raw_id, transcript_id}:
                        if key:
                                transcript_info[key] = info

        normalized = []
        for feature in raw_features:
                if feature.type not in INPUT_FEATURE_TYPES:
                        continue

                direct_gene_id = first_attr(feature, ["gene_id", "gene", "gene_name"])
                direct_transcript_id = first_attr(feature, ["transcript_id", "transcript_name"])
                parent_ids = normalize_parent_ids(feature.attr.get("Parent"))

                transcript_keys = []
                if direct_transcript_id:
                        transcript_keys.append(direct_transcript_id)
                transcript_keys.extend(parent_ids)
                if not transcript_keys:
                        transcript_keys = [None]

                for transcript_key in transcript_keys:
                        info = transcript_info.get(transcript_key, {})
                        transcript_id = direct_transcript_id or info.get("transcript_id") or transcript_key
                        gene_id = direct_gene_id or info.get("gene_id")

                        if not transcript_id or not gene_id:
                                continue

                        normalized_feature = copy.deepcopy(feature)
                        normalized_feature.attr = dict(normalized_feature.attr)
                        normalized_feature.attr["gene_id"] = gene_id
                        normalized_feature.attr["transcript_id"] = transcript_id
                        normalized.append(normalized_feature)

        return normalized

# bad genes are those which have exons in different chromosomes or strands.
def find_bad_genes(gtffile):
        exon_dict = {}
        bad_gene_list = []
        for feature in gtffile:
                if feature.type == "exon" and feature.attr["gene_id"] not in bad_gene_list:
                        gene_id = feature.attr["gene_id"]
                        gene_chrom = feature.iv.chrom
                        gene_strand = feature.iv.strand
                        if gene_id not in exon_dict.keys():
                                exon_dict[gene_id] = (gene_chrom, gene_strand)
                        else:
                                if gene_chrom != exon_dict[gene_id][0] or gene_strand != exon_dict[gene_id][1]:
                                        bad_gene_list.append(gene_id)    
        return bad_gene_list

# Feature is an exon, the transcript region will extend to include the exon region.
def extend_transcript_region(feature, transcript_region):
        gene_id = feature.attr["gene_id"]
        transcript_id = feature.attr["transcript_id"]
        if gene_id not in list(transcript_region.keys()) or transcript_id not in list(transcript_region[gene_id].keys()):
                transcript_region[gene_id][transcript_id] = feature.iv.copy()
        else:
                transcript_iv = transcript_region[gene_id][transcript_id]
                assert transcript_iv.chrom == feature.iv.chrom
                assert transcript_iv.strand == feature.iv.strand
                transcript_region[gene_id][transcript_id].start = min(transcript_region[gene_id][transcript_id].start, feature.iv.start)
                transcript_region[gene_id][transcript_id].end = max(transcript_region[gene_id][transcript_id].end, feature.iv.end)

# Gene region will cover all transcripts region that belong to this gene.                
def find_gene_region(transcript_region):
        gene_region = {}
        for gene_id in list(transcript_region.keys()):
                for transcript_id in list(transcript_region[gene_id].keys()):
                        if gene_id not in list(gene_region.keys()):
                                gene_region[gene_id] = transcript_region[gene_id][transcript_id].copy()
                        else:
                                gene_region[gene_id].start = min(gene_region[gene_id].start, transcript_region[gene_id][transcript_id].start)
                                gene_region[gene_id].end = max(gene_region[gene_id].end, transcript_region[gene_id][transcript_id].end)
        return gene_region

# Gene region length: including all transcripts region length, the overlapping region in different transcripts will only be counted once.
def find_gene_region_length(gene_region, transcript_region):
        gene_region_length = collections.Counter()
        for gene_id in list(transcript_region.keys()):
                transcripts = HTSeq.GenomicArrayOfSets("auto", stranded=True)
                for transcript_id in list(transcript_region[gene_id].keys()):
                        transcripts[transcript_region[gene_id][transcript_id]] += gene_id
                for iv, step_set in transcripts[gene_region[gene_id]].steps():
                        if len(step_set) != 0:
                                gene_region_length[gene_id] += iv.length
        return gene_region_length

# If iv (intron iv) overlaps any of the exon iv in exon_iv_list, return True; Otherwise return False.
def overlap(iv, exon_iv_list):
        for exon_iv in exon_iv_list:
                if iv.overlaps(exon_iv):
                        return True
        return False

def find_region_number(pos, feature_list, start_d_or_end_d_as_pos):
        for feature in feature_list:
                if pos == getattr(feature.iv, start_d_or_end_d_as_pos):
                        region_number_attr_name = feature.type + "_number"
                        return feature.attr[region_number_attr_name]
                
def find_region_number_and_index(pos, feature_list, start_d_or_end_d_as_pos):
        for index in range(len(feature_list)):
                feature = feature_list[index]
                if pos == getattr(feature.iv, start_d_or_end_d_as_pos):
                        region_number_attr_name = feature.type + "_number"
                        return feature.attr[region_number_attr_name], index               

def run(args):
        exons = collections.defaultdict( lambda: HTSeq.GenomicArrayOfSets("auto", stranded=True) )      
        gene_region = {}
        gene_region_length = collections.Counter()
        transcript_region = collections.defaultdict( lambda: dict() )
               
        
        # Read features from the input annotation file (GTF or GFF3) and normalize the attributes used downstream.
        gtffile = HTSeq.GFF_Reader(args.inputfile, end_included=True)
        gtffile = normalize_input_features(list(gtffile))
        bad_gene_list = find_bad_genes(gtffile)
        logging.info("Removing genes with exons in different chromosomes or strands (%i discarded)" % len(bad_gene_list))
        gtffile = filter(lambda feature: feature.attr['gene_id'] not in bad_gene_list, gtffile)
        for feature in gtffile:
                if feature.type == "exon":
                        gene_id = feature.attr["gene_id"]
                        exons[gene_id][feature.iv] += feature.attr["transcript_id"]
                        extend_transcript_region(feature, transcript_region)
                        
        gene_region = find_gene_region(transcript_region)
        gene_region_length = find_gene_region_length(gene_region, transcript_region)
                 
        introns = collections.defaultdict( lambda: HTSeq.GenomicArrayOfSets("auto", stranded=True) )               
        for gene_id in list(transcript_region.keys()):
                for transcript_id in list(transcript_region[gene_id].keys()):
                        transcript_iv = transcript_region[gene_id][transcript_id]
                        for iv, step_set in exons[gene_id][transcript_iv].steps():
                                if transcript_id not in step_set:
                                        introns[gene_id][iv] += transcript_id
        
        # gene_exons_bins redefines the exons in one gene. All exons in the gene region are split into exon bins. 
        # Each exon bin is a feature (feature type is "exonic_region"), which has attributes: "gene_id" and "transcripts". 
        # One exon bin is possibly shared by multiple transcripts.                                
        gene_exons_bins = collections.defaultdict( lambda: list() )

        for gene_id in list(gene_region.keys()):
                gene_iv = gene_region[gene_id]
                gene_iv = gene_region[gene_id]
                for iv, step_set in exons[gene_id][gene_iv].steps():
                        transcript_list = sorted(list(step_set))
                        if len(transcript_list) != 0:
                                feature = HTSeq.GenomicFeature(gene_id, "exonic_region", iv)
                                feature.source = "IR_annotation"
                                feature.attr = {}
                                feature.attr["gene_id"] = gene_id
                                feature.attr["transcripts"] = "+".join(transcript_list)
                                gene_exons_bins[gene_id].append(feature)
                if gene_iv.strand == "-":                
                        gene_exons_bins[gene_id] = gene_exons_bins[gene_id][::-1]
        
        # Number the exon bins with attrubute "exonic_region_number" starting from "001".                           
        for exons_bins_list in list(gene_exons_bins.values()):
                for i in range(len(exons_bins_list)):
                        exons_bins_list[i].attr["exonic_region_number"] = "%03d" % ( i+1 )
        
        # gene_introns_bins redefines the introns in one gene. All introns in the gene region are split into intron bins. 
        # Each intron bin is a feature (feature type is "intronic_region"), which has attributes: "gene_id" and "transcripts". 
        # One intron bin is possibly shared by multiple transcripts. If it isn't shared by one transcript, this intron bin
        # either overlaps with the exonic region of that transcript or lies outside of the whole region of that transcript.
        gene_introns_bins = collections.defaultdict( lambda: list() )
        
        for gene_id in list(gene_region.keys()):
                gene_iv = gene_region[gene_id]
                for iv, step_set in introns[gene_id][gene_iv].steps():
                        transcript_list = sorted(list(step_set))
                        if len(transcript_list) != 0:
                                feature = HTSeq.GenomicFeature(gene_id, "intronic_region", iv)
                                feature.source = "IR_annotation"
                                feature.attr = {}
                                feature.attr["gene_id"] = gene_id
                                feature.attr["transcripts"] = "+".join(transcript_list)
                                gene_introns_bins[gene_id].append(feature)
                if gene_iv.strand == "-":                
                        gene_introns_bins[gene_id] = gene_introns_bins[gene_id][::-1]
        
        # Number the intron bins with attrubute "intronic_region_number" starting from "001".        
        for introns_bins_list in list(gene_introns_bins.values()):
                for i in range(len(introns_bins_list)):
                        introns_bins_list[i].attr["intronic_region_number"] = "%03d" % ( i+1 )                          
        
        
        # gene_constitutive_exons_bins defines that kind of exon bins that shared by all the transcripts in one gene.
        # Each constitutive exon bin is a feature (feature type is "constitutive_exonic_region"), which has attribute: "gene_id". 
        logging.info("Generating constitutive exonic region (CER) annotation")
        
        gene_constitutive_exons_bins = collections.defaultdict( lambda: list() )
        gene_constitutive_exons_start_d = collections.defaultdict( lambda: set() )
        gene_constitutive_exons_end_d = collections.defaultdict( lambda: set() )
        gene_constitutive_exons_length = collections.Counter()
        gene_constitutive_exons_number = collections.Counter()
        
        for gene_id in list(gene_region.keys()):
                transcripts_in_gene = len(transcript_region[gene_id])
                gene_iv = gene_region[gene_id]
                for iv, step_set in exons[gene_id][gene_iv].steps():
                        transcript_list = list(step_set)
                        if len(transcript_list) == transcripts_in_gene:
                                feature = HTSeq.GenomicFeature(gene_id, "constitutive_exonic_region", iv)
                                feature.source = "IR_annotation"
                                feature.attr = {}
                                feature.attr["gene_id"] = gene_id
                                gene_constitutive_exons_bins[gene_id].append(feature)
                                gene_constitutive_exons_start_d[gene_id].add(feature.iv.start_d_as_pos)
                                gene_constitutive_exons_end_d[gene_id].add(feature.iv.end_d_as_pos)
                                gene_constitutive_exons_length[gene_id] += feature.iv.length
                                gene_constitutive_exons_number[gene_id] += 1
                if gene_iv.strand == "-":
                        gene_constitutive_exons_bins[gene_id] = gene_constitutive_exons_bins[gene_id][::-1]
                                
        for constitutive_exons_bins_list in list(gene_constitutive_exons_bins.values()):
                for i in range(len(constitutive_exons_bins_list)):
                        constitutive_exons_bins_list[i].attr["constitutive_exonic_region_number"] = "%03d" % ( i+1 )
        
        
        # gene_constitutive_introns_bins defines that kind of intron bins that shared by all the transcripts in one gene.
        # Each constitutive intron bin is a feature (feature type is "constitutive_intronic_region"), which has attribute: "gene_id".
        logging.info("Generating constitutive intronic region (CIR) annotation")
        
        gene_constitutive_introns_bins = collections.defaultdict( lambda: list() )
        gene_constitutive_introns_start_d = collections.defaultdict( lambda: set() )
        gene_constitutive_introns_end_d = collections.defaultdict( lambda: set() )
        gene_constitutive_introns_length = collections.Counter()
        gene_constitutive_introns_number = collections.Counter()
        
        for gene_id in list(gene_region.keys()):
                transcripts_in_gene = len(transcript_region[gene_id])
                gene_iv = gene_region[gene_id]
                for iv, step_set in introns[gene_id][gene_iv].steps():
                        transcript_list = list(step_set)
                        if len(transcript_list) == transcripts_in_gene:
                                feature = HTSeq.GenomicFeature(gene_id, "constitutive_intronic_region", iv)
                                feature.source = "IR_annotation"
                                feature.attr = {}
                                feature.attr["gene_id"] = gene_id
                                gene_constitutive_introns_bins[gene_id].append(feature)
                                gene_constitutive_introns_start_d[gene_id].add(feature.iv.start_d_as_pos)
                                gene_constitutive_introns_end_d[gene_id].add(feature.iv.end_d_as_pos)
                                gene_constitutive_introns_length[gene_id] += feature.iv.length
                                gene_constitutive_introns_number[gene_id] += 1
                if gene_iv.strand == "-":
                        gene_constitutive_introns_bins[gene_id] = gene_constitutive_introns_bins[gene_id][::-1]
        
        for constitutive_introns_bins_list in list(gene_constitutive_introns_bins.values()):
                for i in range(len(constitutive_introns_bins_list)):
                        constitutive_introns_bins_list[i].attr["constitutive_intronic_region_number"] = "%03d" % ( i+1 )

        all_constitutive_introns = []
        for constitutive_introns_bins_list in list(gene_constitutive_introns_bins.values()):
                all_constitutive_introns.extend(constitutive_introns_bins_list)

        gene_constitutive_introns_effective_length = collections.Counter()
        if args.mapfile:
                blacklist_bed = resolve_blacklist_bed(
                        args.mapfile,
                        args.inputfile,
                        all_constitutive_introns,
                        args.mapreadlength,
                        args.mapcutoff,
                        getattr(args, "mapcachedir", None),
                )
                effective_length_map = effective_lengths_from_blacklist(all_constitutive_introns, blacklist_bed)
        else:
                effective_length_map = {
                        (feature.attr["gene_id"], feature.attr["constitutive_intronic_region_number"]): feature.iv.length
                        for feature in all_constitutive_introns
                }

        for feature in all_constitutive_introns:
                key = (feature.attr["gene_id"], feature.attr["constitutive_intronic_region_number"])
                effective_length = effective_length_map[key]
                feature.attr["constitutive_intronic_region_effective_length"] = effective_length
                gene_constitutive_introns_effective_length[feature.attr["gene_id"]] += effective_length
                                
                                                        
        # gene_constitutive_junction defines that kind of junction positions that join constitutive exon bin and constitutive intron bin in one gene.   
        # Each constitutive junction is a feature (feature type is "constitutive_junction"), which has attributes: "gene_id", "constitutive_junction_type", "upstream", "downstream".
        # attr["constitutive_junction_type"] can be the value: "5'_splice_junction", which means the upstream of the junction position is a constitutive exon bin,
        # and the downstream of the junction position is a constitutive intron bin. In this case, attr["upstream"] will be like "constitutive_exonic_region_number 002" (shows exactly which
        # constitutive exon bin in the upstream), and similarly for attr["downstream"].
        # On the other hand, attr["constitutive_junction_type"] can be the value: "3'_splice_junction"
        logging.info("Generating constitutive junctions (CJ) annotation")
        
        gene_constitutive_junction = collections.defaultdict( lambda: list() )
        for gene_id in list(gene_constitutive_exons_start_d.keys()):
                if gene_id in list(gene_constitutive_introns_start_d.keys()):
                        gene_constitutive_junction_from_exon_to_intron_set = gene_constitutive_exons_end_d[gene_id] & gene_constitutive_introns_start_d[gene_id]
                        for gene_constitutive_junction_pos in gene_constitutive_junction_from_exon_to_intron_set:
                                feature = HTSeq.GenomicFeature(gene_id, "constitutive_junction", gene_constitutive_junction_pos)
                                feature.source = "IR_annotation"
                                feature.attr = {}
                                feature.attr["gene_id"] = gene_id
                                from_region_number = find_region_number(gene_constitutive_junction_pos, gene_constitutive_exons_bins[gene_id], "end_d_as_pos")
                                feature.attr["upstream"] =  "constitutive_exonic_region_number " + from_region_number
                                to_region_number, index = find_region_number_and_index(gene_constitutive_junction_pos, gene_constitutive_introns_bins[gene_id], "start_d_as_pos")
                                feature.attr["downstream"] =  "constitutive_intronic_region_number " + to_region_number
                                feature.attr["constitutive_junction_type"] = "5'_splice_junction"
                                gene_constitutive_junction[gene_id].append((feature, index))
                                
                                
                        gene_constitutive_junction_from_intron_to_exon_set = gene_constitutive_exons_start_d[gene_id] & gene_constitutive_introns_end_d[gene_id]
                        for gene_constitutive_junction_pos in gene_constitutive_junction_from_intron_to_exon_set:
                                feature = HTSeq.GenomicFeature(gene_id, "constitutive_junction", gene_constitutive_junction_pos)
                                feature.source = "IR_annotation"
                                feature.attr = {}
                                feature.attr["gene_id"] = gene_id
                                from_region_number, index = find_region_number_and_index(gene_constitutive_junction_pos, gene_constitutive_introns_bins[gene_id], "end_d_as_pos")
                                feature.attr["upstream"] =  "constitutive_intronic_region_number " + from_region_number
                                to_region_number = find_region_number(gene_constitutive_junction_pos, gene_constitutive_exons_bins[gene_id], "start_d_as_pos")
                                feature.attr["downstream"] =  "constitutive_exonic_region_number " + to_region_number
                                feature.attr["constitutive_junction_type"] = "3'_splice_junction"
                                gene_constitutive_junction[gene_id].append((feature, index))
                                
                                
                if len(gene_constitutive_junction[gene_id]) > 0:
                        gene_strand = gene_constitutive_junction[gene_id][0][0].iv.strand
                        if gene_strand == "+":
                                gene_constitutive_junction[gene_id].sort(key = lambda f: ( f[0].iv.chrom, f[0].iv.start ))
                        else:
                                gene_constitutive_junction[gene_id].sort(key = lambda f: ( f[0].iv.chrom, -f[0].iv.start ))
        
        for gene_constitutive_junction_list in list(gene_constitutive_junction.values()):
                for i in range(len(gene_constitutive_junction_list)):
                        gene_constitutive_junction_list[i][0].attr["constitutive_junction_number"] = "%03d" % ( i+1 ) 
                        
                        feature = gene_constitutive_junction_list[i][0]
                        gene_id = feature.attr["gene_id"]
                        constitutive_junction_type = feature.attr["constitutive_junction_type"]
                        if constitutive_junction_type == "5'_splice_junction":
                                index = gene_constitutive_junction_list[i][1]
                                gene_constitutive_introns_bins[gene_id][index].attr["upstream_constitutive_junction_number"] =  "%03d" % ( i+1 ) 
                        elif constitutive_junction_type == "3'_splice_junction":
                                index = gene_constitutive_junction_list[i][1]
                                gene_constitutive_introns_bins[gene_id][index].attr["downstream_constitutive_junction_number"] =  "%03d" % ( i+1 ) 
                                
        
        # gene_region_features defines a feature for each gene that summarize some key info. feature iv is the gene region.
        # Feature type is "gene_region". Each feature has attributes: "gene_id", "transcripts_in_gene" (count how many transcripts in this gene), 
        # "gene_region_length", "constitutive_exonic_region_length", "constitutive_intronic_region_length". 
        gene_region_features = []
        for gene_id in list(gene_region.keys()):
                iv = gene_region[gene_id]
                feature = HTSeq.GenomicFeature(gene_id, "gene_region", iv)
                feature.source = "IR_annotation"
                feature.attr = {}
                feature.attr["gene_id"] = gene_id
                feature.attr["transcripts_in_gene"] = len(transcript_region[gene_id])
                feature.attr["gene_region_length"] = gene_region_length[gene_id]
                feature.attr["constitutive_exonic_region_length"] = gene_constitutive_exons_length[gene_id]
                feature.attr["constitutive_intronic_region_length"] = gene_constitutive_introns_length[gene_id]
                feature.attr["constitutive_intronic_region_effective_length"] = gene_constitutive_introns_effective_length[gene_id]
                feature.attr["constitutive_exonic_region_number"] = gene_constitutive_exons_number[gene_id]
                feature.attr["constitutive_intronic_region_number"] = gene_constitutive_introns_number[gene_id]
                gene_region_features.append(feature)
                
        gene_region_features.sort( key = lambda f: ( f.iv.chrom, f.iv.start ) )
        
        # transcript_region_features defines a feature for each transcript. feature iv is the transcript region.
        # Feature type is "transcript_region". Each feature has attributes: "gene_id", "transcript_id".
        transcript_region_features = collections.defaultdict( lambda: list() )
        for gene_id in list(transcript_region.keys()):
                for transcript_id in list(transcript_region[gene_id].keys()):
                        iv = transcript_region[gene_id][transcript_id]
                        feature = HTSeq.GenomicFeature(gene_id, "transcript_region", iv)
                        feature.source = "IR_annotation"
                        feature.attr = {}
                        feature.attr["gene_id"] = gene_id
                        feature.attr["transcript_id"] = transcript_id
                        transcript_region_features[gene_id].append(feature)
        
        # Write all newly defined features into new gtf annotation file.
        logging.info("Writing annotation to file: %s" % os.path.join(args.outdir, args.annofile))
        
        f = open(os.path.join(args.outdir, args.annofile), "w")
        for gene_region_feature in gene_region_features:
                f.write( gene_region_feature.get_gff_line() )
                gene_id = gene_region_feature.attr["gene_id"]
                
                for feature in transcript_region_features[gene_id]:
                        f.write( feature.get_gff_line() )
                
                for feature in gene_exons_bins[gene_id]:
                        f.write( feature.get_gff_line() )
                        
                for feature in gene_introns_bins[gene_id]:
                        f.write( feature.get_gff_line() )
                        
                for feature in gene_constitutive_exons_bins[gene_id]:
                        f.write( feature.get_gff_line() )
                        
                for feature in gene_constitutive_introns_bins[gene_id]:
                        f.write( feature.get_gff_line() )  
                        
                for feature in [item[0] for item in gene_constitutive_junction[gene_id]]:
                        f.write( feature.get_gff_line() ) 
                        
        f.close()
        
                
                                
                        
                
                                        
        
                                        
        
