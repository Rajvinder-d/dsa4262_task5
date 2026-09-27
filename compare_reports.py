#!/usr/bin/env python3
"""Compare Task 5 runtimes from two native Nextflow report.html files."""
import argparse
import html
import json
from pathlib import Path
import re


def duration_seconds(value):
    scales = {"d": 86400, "h": 3600, "m": 60, "s": 1, "ms": 0.001}
    pattern = r"(\d+(?:\.\d+)?)\s*(ms|d|h|m|s)"
    matches = list(re.finditer(pattern, value))
    if not matches or re.sub(pattern, "", value).strip():
        raise ValueError(f"Unrecognized Nextflow duration: {value!r}")
    return sum(float(m[1]) * scales[m[2]] for m in matches)


def read_report(path, mode):
    source = path.read_text()
    if "Workflow execution completed successfully!" not in source:
        raise ValueError(f"{path}: report does not show successful completion")
    match = re.search(r"duration:\s*<strong>([^<]+)</strong>", source)
    if not match:
        raise ValueError(f"{path}: cannot locate Summary duration; read it in your browser")
    duration = html.unescape(match[1]).strip()
    marker = source.rfind("// Nextflow report data")
    payload_match = re.search(r"window\.data\s*=\s*", source[marker:]) if marker >= 0 else None
    if not payload_match:
        raise ValueError(f"{path}: cannot locate Nextflow task data; use the Tasks table")
    payload = source[marker + payload_match.end():]
    # Nextflow escapes apostrophes for JavaScript; JSON does not allow \\'.
    # Consume escape pairs so genuine escaped backslashes remain intact.
    payload = re.sub(r"\\(.)", lambda m: "'" if m[1] == "'" else m[0], payload)
    data, _ = json.JSONDecoder().raw_decode(payload)
    tasks = data.get("trace") or []
    bambu = [t for t in tasks if "BAMBU_DISCOVERY" in t.get("name", "")]
    if len(bambu) != 1 or f"({mode})" not in bambu[0]["name"]:
        raise ValueError(f"{path}: expected one BAMBU_DISCOVERY ({mode}) task")
    if bambu[0]["status"] != "COMPLETED":
        raise ValueError(f"{path}: Bambu was cached or incomplete; use the original scenario report")
    return {"duration": duration, "seconds": duration_seconds(duration),
            "bambu_seconds": float(bambu[0]["realtime"]) / 1000,
            "cached": sum(t["status"] == "CACHED" for t in tasks),
            "executed": sum(t["status"] == "COMPLETED" for t in tasks)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("annotated_report", type=Path)
    parser.add_argument("unannotated_report", type=Path)
    args = parser.parse_args()
    a = read_report(args.annotated_report, "annotated")
    b = read_report(args.unannotated_report, "unannotated")
    print("| Measurement from report.html | Scenario 1: annotated | Scenario 2: unannotated + resume |")
    print("| --- | ---: | ---: |")
    print(f'| Workflow duration (Summary) | {a["duration"]} | {b["duration"]} |')
    print(f'| Bambu task realtime (seconds) | {a["bambu_seconds"]:.3f} | {b["bambu_seconds"]:.3f} |')
    print(f'| Completed tasks | {a["executed"]} | {b["executed"]} |')
    print(f'| Cached tasks | {a["cached"]} | {b["cached"]} |')
    delta = b["seconds"] - a["seconds"]
    print(f"\nWorkflow runtime difference (scenario 2 minus scenario 1): {delta:+.3f} seconds.")
    if a["seconds"]:
        print(f'Workflow runtime reduction relative to scenario 1: {-100 * delta / a["seconds"]:.2f}% (negative means slower).')
    print(f'Bambu task realtime difference (scenario 2 minus scenario 1): {b["bambu_seconds"] - a["bambu_seconds"]:+.3f} seconds.')
    print("\nWorkflow durations have the precision displayed in report.html. "
          "Scenario 2 reuses reference preparation, alignment and BAM conversion; "
          "the whole-workflow difference includes caching as well as the annotation change.")
    if (a["cached"], a["executed"], b["cached"], b["executed"]) != (0, 11, 9, 2):
        print("\nCheck the Tasks tables: expected 11 completed tasks in scenario 1, "
              "then 9 cached and 2 completed tasks in scenario 2.")


if __name__ == "__main__":
    main()
