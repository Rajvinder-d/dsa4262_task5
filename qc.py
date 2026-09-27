#!/usr/bin/env python3
"""Read-only FASTQ/BAM QC and an offline report; Python standard library only."""
import argparse
import bisect
import csv
import gzip
import html
import itertools
import json
import math
import re
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
import subprocess
import sys
import zlib


def write_tsv(path, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def fastq_metrics(path, short_length):
    lengths = Counter()
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="ascii") as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            seq, plus, qual = (handle.readline() for _ in range(3))
            if not header.startswith("@") or not plus.startswith("+") or not qual:
                raise ValueError("Malformed or incomplete four-line FASTQ record")
            seq, qual = seq.rstrip("\r\n"), qual.rstrip("\r\n")
            if not seq or len(seq) != len(qual):
                raise ValueError("Empty sequence or sequence/quality length mismatch")
            lengths[len(seq)] += 1
    if not lengths:
        raise ValueError("FASTQ contains no reads")
    ordered = sorted(lengths)
    cumulative = list(itertools.accumulate(lengths[n] for n in ordered))
    count = cumulative[-1]

    def quantile(p):
        rank = (count - 1) * p
        low, high = math.floor(rank), math.ceil(rank)
        a = ordered[bisect.bisect_right(cumulative, low)]
        b = ordered[bisect.bisect_right(cumulative, high)]
        return a + (b - a) * (rank - low)

    bases = sum(n * count_at_n for n, count_at_n in lengths.items())
    return {
        "fastq_reads": count, "fastq_bases": bases,
        "min_length": ordered[0], "q1_length": quantile(0.25),
        "median_length": quantile(0.5), "mean_length": round(bases / count, 3),
        "q3_length": quantile(0.75), "max_length": ordered[-1],
        "short_read_pct": round(100 * sum(v for n, v in lengths.items()
                                          if n < short_length) / count, 3)
    }, lengths


def alignment_metrics(path, raw_dir, sample):
    check = subprocess.run(["samtools", "quickcheck", "-v", str(path)],
                           capture_output=True, text=True)
    (raw_dir / f"{sample}.quickcheck.txt").write_text(check.stdout + check.stderr)
    stats = subprocess.run(["samtools", "flagstat", "-@", "1", "-O", "json", str(path)],
                           capture_output=True, text=True)
    (raw_dir / f"{sample}.flagstat.stderr.txt").write_text(stats.stderr)
    if stats.returncode:
        raise ValueError("samtools flagstat could not read the BAM; see its stderr log")
    data = json.loads(stats.stdout)
    (raw_dir / f"{sample}.flagstat.json").write_text(stats.stdout)

    def total(key):
        return sum(data[group][key] for group in ("QC-passed reads", "QC-failed reads"))

    primary, mapped = total("primary"), total("primary mapped")
    return {
        "bam_quickcheck": "PASS" if check.returncode == 0 else "FAIL",
        "alignment_records": total("total"), "primary_reads": primary,
        "primary_mapped": mapped,
        "primary_mapped_pct": round(100 * mapped / primary, 3) if primary else None,
        "secondary_alignments": total("secondary"),
        "supplementary_alignments": total("supplementary")
    }


def sample_flags(row, args):
    flags = []
    if row["fastq_status"] != "PASS":
        flags.append("FASTQ_READ_ERROR")
    if row["bam_quickcheck"] != "PASS":
        flags.append("BAM_CHECK_ERROR")
    if row["fastq_reads"] is not None and row["primary_reads"] is not None:
        if row["fastq_reads"] != row["primary_reads"]:
            flags.append("FASTQ_BAM_COUNT_MISMATCH")
    if row["primary_reads"] is not None and row["primary_reads"] < args.min_primary_reads:
        flags.append("LOW_PRIMARY_READ_COUNT")
    if row["primary_mapped_pct"] is None:
        flags.append("MAPPING_RATE_UNAVAILABLE")
    elif row["primary_mapped_pct"] < args.min_mapping_pct:
        flags.append("LOW_PRIMARY_MAPPING_RATE")
    if row["short_read_pct"] is not None and row["short_read_pct"] > args.max_short_pct:
        flags.append("MOST_READS_SHORT")
    if float(row["transcripts_ge10"]) == 0:
        flags.append("NO_TRANSCRIPTS_GE10")
    return flags


def table(rows, fields=None):
    fields = fields or list(rows[0])
    esc = lambda x: html.escape("NA" if x is None else str(x))
    head = "".join(f"<th>{esc(x)}</th>" for x in fields)
    body = "".join("<tr>" + "".join(f"<td>{esc(row[x])}</td>" for x in fields)
                   + "</tr>" for row in rows)
    return f'<div class="scroll"><table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table></div>'


def make_report(out, rows, args):
    totals = read_tsv(out / "transcript_totals.tsv")
    classes = read_tsv(out / "transcripts_by_class.tsv")
    warnings = html.escape((out / "bambu_warnings.txt").read_text())
    now = datetime.now(timezone.utc).isoformat(timespec="seconds")
    data = {
        "generated_utc": now, "reporting_only": True,
        "annotation_mode": args.annotation_mode,
        "policy": {"min_primary_reads": args.min_primary_reads,
                   "min_primary_mapping_pct": args.min_mapping_pct,
                   "short_read_length_lt": args.short_length,
                   "max_short_read_pct": args.max_short_pct},
        "samples": rows, "transcript_totals": totals, "transcripts_by_class": classes
    }
    (out / "qc.json").write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
    summary_fields = ["sample", "qc_status", "flags", "fastq_reads", "mean_length",
                      "median_length", "short_read_pct", "primary_mapped_pct",
                      "estimated_counts", "transcripts_ge10"]
    errors = "\n".join(f'{row["sample"]}: {row["detail"]}' for row in rows if row["detail"])
    page = f'''<!doctype html>
<html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>SG-NEx QC report</title><style>
body{{font:15px/1.5 system-ui,sans-serif;color:#172b3a;margin:32px auto;padding:0 24px;max-width:1500px}}
h1,h2{{color:#123e51}}h2{{margin-top:32px}}.scroll{{overflow-x:auto}}table{{border-collapse:collapse;width:100%}}
th,td{{padding:9px 12px;text-align:left;border-bottom:1px solid #dce3e8;white-space:nowrap}}
th{{background:#edf3f5}}.note{{padding:16px;background:#fff4d8;border-left:4px solid #bd8112}}
pre{{white-space:pre-wrap;overflow-wrap:anywhere;background:#f4f6f8;padding:16px}}a{{color:#075c85}}
</style><h1>SG-NEx QC report</h1><p>Generated {now}. Metrics use the supplied final Bambu object.</p>
<p>Annotation mode: <strong>{args.annotation_mode}</strong>. In unannotated mode, novel means
discovered without a supplied reference transcript annotation. Those totals are not novelty
relative to the reference GTF and should not be compared directly with annotated-mode novelty.</p>
<p class="note">Reporting only: all samples and transcripts are retained. REVIEW means a reporting rule
was triggered. ERROR means a file could not be checked or parsed successfully. NO_FLAGS means none of
these rules triggered; it is not a comprehensive biological quality assessment.</p>
<h2>Sample overview</h2>{table(rows, summary_fields)}
<h2>Transcript totals before additional QC filtering</h2>{table(totals)}
<p>The novel_single_exon_transcripts value is the total requested for the assignment. It has no
additional estimated-count, unique-count or full-length-count threshold.</p>
<h2>Transcript support by class</h2>{table(classes)}
<p>detected_ge1 and expressed_ge10 require the stated estimated count in any sample.
full_length_ge1 and unique_ge1 use pooled support across all samples. These columns are independent
descriptions, not simultaneous filters. Single-exon models lack splice-junction evidence; counts
alone do not establish a complete biological transcript.</p>
<h2>Reporting rules</h2><ul>
<li>Primary reads below {args.min_primary_reads:,}.</li>
<li>Primary mapping below {args.min_mapping_pct:g}%.</li>
<li>More than {args.max_short_pct:g}% of reads shorter than {args.short_length} nt.</li>
<li>FASTQ count differs from primary BAM count, a file check fails, or no transcript reaches 10 counts.</li>
</ul><p>These configurable rules are workflow heuristics, not validated universal cutoffs.
Primary mapping combines SAM QC-pass and QC-fail records. Secondary/supplementary values count
alignment records, not independent reads. quickcheck tests headers and file termination; flagstat
then reads through the BAM. Sample names must match across FASTQ, BAM and se.</p>
<h2>Bambu warnings preserved from se</h2><pre>{warnings}</pre>
<h2>File-check details</h2><pre>{html.escape(errors or "No file-check errors recorded.")}</pre>
<h2>Data files</h2><ul>
<li><a href="sample_qc.tsv">Complete sample metrics</a></li>
<li><a href="transcript_totals.tsv">Transcript totals</a></li>
<li><a href="transcripts_by_class.tsv">Transcript support by class</a></li>
<li><a href="bambu_samples.tsv">Bambu sample metrics</a></li>
<li><a href="qc.json">Machine-readable report and rules</a></li>
</ul><p>read_lengths/ contains exact per-sample read-length frequency tables.
alignment/ contains raw samtools results. No transcripts or samples were removed.</p></html>'''
    (out / "qc_report.html").write_text(page)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reads-dir", type=Path, required=True)
    parser.add_argument("--bam-dir", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--annotation-mode", choices=["annotated", "unannotated"], required=True)
    parser.add_argument("--min-primary-reads", type=int, default=10000)
    parser.add_argument("--min-mapping-pct", type=float, default=50)
    parser.add_argument("--short-length", type=int, default=200)
    parser.add_argument("--max-short-pct", type=float, default=90)
    args = parser.parse_args()
    if (args.min_primary_reads < 0 or args.short_length < 1
            or not 0 <= args.min_mapping_pct <= 100 or not 0 <= args.max_short_pct <= 100):
        parser.error("Invalid reporting threshold")
    out = args.outdir
    raw_dir, length_dir = out / "alignment", out / "read_lengths"
    raw_dir.mkdir(parents=True, exist_ok=True)
    length_dir.mkdir(parents=True, exist_ok=True)
    bambu_rows = read_tsv(out / "bambu_samples.tsv")
    names = {r["sample"] for r in bambu_rows}
    reads = {re.sub(r"\.(fastq|fq)(\.gz)?$", "", p.name): p
             for p in args.reads_dir.iterdir()
             if re.search(r"\.(fastq|fq)(\.gz)?$", p.name)}
    bams = {p.stem: p for p in args.bam_dir.glob("*.bam")}
    if len(names) != len(bambu_rows) or names != set(reads) or names != set(bams):
        raise ValueError("Sample names must match exactly across se, FASTQs and BAMs")
    rows = []
    for bambu in bambu_rows:
        sample = bambu["sample"]
        row = {"sample": sample, "qc_status": "", "flags": "", "fastq_status": "PASS"}
        row.update(dict.fromkeys(["fastq_reads", "fastq_bases", "min_length", "q1_length",
                                  "median_length", "mean_length", "q3_length", "max_length",
                                  "short_read_pct"]))
        row.update({"bam_quickcheck": "ERROR", "alignment_records": None,
                    "primary_reads": None, "primary_mapped": None, "primary_mapped_pct": None,
                    "secondary_alignments": None, "supplementary_alignments": None})
        details = []
        try:
            metrics, lengths = fastq_metrics(reads[sample], args.short_length)
            row.update(metrics)
            write_tsv(length_dir / f"{sample}.read_lengths.tsv",
                      [{"length_bp": n, "reads": lengths[n]} for n in sorted(lengths)])
        except (ValueError, OSError, EOFError, UnicodeError, zlib.error) as error:
            row["fastq_status"] = "ERROR"
            details.append(str(error))
        try:
            row.update(alignment_metrics(bams[sample], raw_dir, sample))
        except (ValueError, KeyError) as error:
            details.append(str(error))
        for key, value in bambu.items():
            if key != "sample":
                row[key] = float(value) if "." in value or "e" in value.lower() else int(value)
        row["detail"] = "; ".join(details)
        flags = sample_flags(row, args)
        row["flags"] = ";".join(flags) if flags else "none"
        file_error = row["fastq_status"] != "PASS" or row["bam_quickcheck"] != "PASS"
        row["qc_status"] = "ERROR" if file_error else ("REVIEW" if flags else "NO_FLAGS")
        rows.append(row)
        print(f'{sample}: {row["qc_status"]} ({row["flags"]})', flush=True)
    write_tsv(out / "sample_qc.tsv", rows)
    make_report(out, rows, args)
    versions = subprocess.run(["samtools", "--version"], capture_output=True, text=True, encoding="utf-8", errors="replace", check=True)
    (out / "software_versions.txt").write_text(f"Python {sys.version}\n{versions.stdout}")
    print(f"Report written to {out / 'qc_report.html'}", flush=True)


if __name__ == "__main__":
    main()
