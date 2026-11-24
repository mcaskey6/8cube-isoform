"""Keep barcodes with at least MIN_UMIS UMIs in a kb cells x TCC MatrixMarket matrix.

Usage: filter_barcodes.py MIN_UMIS IN.mtx IN.barcodes.txt OUT.mtx OUT.barcodes.txt

Rows are barcodes (1-based, in the order of the barcodes file); values are UMI counts per
equivalence class. Kept rows are renumbered in their original order, so OUT.barcodes.txt
lines up with OUT.mtx. Columns (equivalence classes) are unchanged, so the original
cells_x_tcc.ec.txt still applies.
"""
import sys
from collections import Counter


def data_lines(path):
    """Yield the MatrixMarket size line, then each entry line (comments skipped)."""
    with open(path) as f:
        for line in f:
            if not line.startswith("%"):
                yield line


def main(min_umis, mtx_in, barcodes_in, mtx_out, barcodes_out):
    min_umis = float(min_umis)

    lines = data_lines(mtx_in)
    _, n_cols, _ = next(lines).split()
    umis = Counter()
    for line in lines:
        row, _, value = line.split()
        umis[int(row)] += float(value)

    with open(barcodes_in) as f:
        barcodes = f.read().split()
    kept = [row for row in range(1, len(barcodes) + 1) if umis[row] >= min_umis]
    new_row = {row: i for i, row in enumerate(kept, start=1)}

    with open(barcodes_out, "w") as out:
        out.writelines(barcodes[row - 1] + "\n" for row in kept)

    lines = data_lines(mtx_in)
    next(lines)
    entries = [line.split() for line in lines]
    entries = [(new_row[int(r)], c, v) for r, c, v in entries if int(r) in new_row]
    with open(mtx_out, "w") as out:
        out.write("%%MatrixMarket matrix coordinate real general\n")
        out.write(f"{len(kept)} {n_cols} {len(entries)}\n")
        out.writelines(f"{r} {c} {v}\n" for r, c, v in entries)

    total = sum(umis.values())
    kept_umis = sum(umis[row] for row in kept)
    print(f"Kept {len(kept):,} of {len(barcodes):,} barcodes with >= {min_umis:g} UMIs "
          f"({kept_umis / total:.1%} of UMIs)", file=sys.stderr)


if __name__ == "__main__":
    main(*sys.argv[1:6])
