"""Restore the phase (frame) column of coding features in a Liftoff GTF.

Liftoff writes "." as the phase of every feature, and tools that build CDS frames
(e.g. gtfToGenePred -genePredExt, used by SQANTI3) reject that. The phase of each
CDS / start_codon / stop_codon segment follows from the length of the preceding
segments of the same feature type in transcription order:
    next_phase = (phase - length) mod 3
The first segment's phase is taken from the source (GENCODE) annotation, which is
non-zero for transcripts with an incomplete 5' CDS (cds_start_NF); otherwise 0.

Start/stop codons of transcripts left without any CDS by a partial lift-over are
dropped: gtfToGenePred counts them as coding sequence and then finds no CDS frames.
"""
import re
import sys
from collections import defaultdict

CODING = {"CDS", "start_codon", "stop_codon"}
TRANSCRIPT_ID = re.compile(r'transcript_id "([^"]+)"')


def coding_segments(path):
    """Yield (line_number, fields, transcript_id) for coding features of a GTF."""
    with open(path) as f:
        for n, line in enumerate(f):
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) == 9 and fields[2] in CODING:
                yield n, fields, TRANSCRIPT_ID.search(fields[8]).group(1)


def group(path):
    """(transcript_id, feature) -> [(line_number, start, end, strand, phase)]."""
    groups = defaultdict(list)
    for n, f, tid in coding_segments(path):
        groups[(tid, f[2])].append((n, int(f[3]), int(f[4]), f[6], f[7]))
    return groups


def transcription_order(segments):
    return sorted(segments, key=lambda s: s[1] if s[3] == "+" else -s[2])


def first_phases(path):
    """Phase of the first segment (in transcription order) of each coding feature."""
    return {key: transcription_order(segs)[0][4] for key, segs in group(path).items()}


def compute_phases(groups, initial):
    phases = {}
    for key, segments in groups.items():
        start_phase = initial.get(key, "0")
        phase = int(start_phase) if start_phase in ("0", "1", "2") else 0
        for n, start, end, _, _ in transcription_order(segments):
            phases[n] = phase
            phase = (phase - (end - start + 1)) % 3
    return phases


def main(source_gtf, lifted_gtf, out_gtf):
    groups = group(lifted_gtf)
    with_cds = {tid for tid, feature in groups if feature == "CDS"}
    orphan_codons = {n for (tid, _), segments in groups.items() if tid not in with_cds
                     for n, *_ in segments}
    phases = compute_phases(groups, first_phases(source_gtf))
    with open(lifted_gtf) as f, open(out_gtf, "w") as out:
        for n, line in enumerate(f):
            if n in orphan_codons:
                continue
            if n in phases:
                fields = line.rstrip("\n").split("\t")
                fields[7] = str(phases[n])
                line = "\t".join(fields) + "\n"
            out.write(line)
    print(f"Set the phase of {len(phases) - len(orphan_codons)} coding segments; dropped "
          f"{len(orphan_codons)} start/stop codons of transcripts without CDS", file=sys.stderr)


if __name__ == "__main__":
    if "snakemake" in globals():
        main(snakemake.input.source_gtf, snakemake.input.lifted_gtf, snakemake.output[0])
    else:
        main(*sys.argv[1:4])
