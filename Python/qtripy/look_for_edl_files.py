#!/usr/bin/env python3
"""
Scan a directory for .vedl/.aedl files, print size, header hex (first 64 bytes)
and basic heuristics (first 2 uint32 dims, possible float32 payload shape,
ASCII numeric file probe, MATLAB v5 header).
Usage:
    python3 check_vedl_aedl.py [path]
If path omitted, current directory is used.
"""
import argparse
import os
import sys
from pathlib import Path
import textwrap

data_path = os.environ.get('env_data_path')

def hexdump(b: bytes, width=16):
    hexs = b.hex()
    out = []
    for i in range(0, len(b), width):
        chunk = b[i:i+width]
        hex_bytes = ' '.join(f"{x:02x}" for x in chunk)
        out.append(f"{i:04x}: {hex_bytes}")
    return "\n".join(out)

def probe_file(fn: Path):
    info = {}
    info['path'] = str(fn)
    info['size'] = fn.stat().st_size
    with fn.open('rb') as f:
        hdr = f.read(128)
    info['hdr64'] = hdr[:64]
    # first two uint32 little-endian
    if len(hdr) >= 8:
        r0 = int.from_bytes(hdr[0:4], 'little', signed=False)
        r1 = int.from_bytes(hdr[4:8], 'little', signed=False)
        info['uint32_0_1'] = (r0, r1, r0 * r1)
    else:
        info['uint32_0_1'] = None
    # ASCII numeric heuristic
    try:
        txt = info['hdr64'].decode('ascii', errors='ignore')
        info['ascii_like'] = ('\n' in txt or '\r' in txt) and any(ch.isdigit() for ch in txt)
    except Exception:
        info['ascii_like'] = False
    info['matlab5'] = info['hdr64'].startswith(b'MATLAB 5.0 MAT-file')
    # quick float32 payload detection (8-byte dims header + float32 body)
    if info['uint32_0_1'] is not None:
        r0, r1, prod = info['uint32_0_1']
        if info['size'] == 8 + prod * 4 and prod > 0:
            info['likely_float32_matrix'] = (r0, r1)
        else:
            info['likely_float32_matrix'] = None
    else:
        info['likely_float32_matrix'] = None
    return info

def try_load_shapes(fn: Path):
    """Try to infer shapes by attempting to read as several dtypes (small probe)."""
    import numpy as np
    out = []
    size = fn.stat().st_size
    with fn.open('rb') as f:
        data = f.read()
    # try uint8 length -> trivial
    out.append(("bytes", (size,)))
    # try interpreting as float32 array (no header)
    if size % 4 == 0 and size >= 4:
        arr32_count = size // 4
        k = int((arr32_count)**0.5)
        if k*k == arr32_count:
            out.append(("float32 square", (k, k)))
        # common shaped guesses
        for ncols in (1544, 3088, 551, 550):
            if arr32_count % ncols == 0:
                out.append((f"float32 -> (rows,cols={ncols})", (arr32_count//ncols, ncols)))
    # try float64
    if size % 8 == 0 and size >= 8:
        arr64_count = size // 8
        k = int((arr64_count)**0.5)
        if k*k == arr64_count:
            out.append(("float64 square", (k, k)))
    return out

def main():
    p = argparse.ArgumentParser(description="Scan .vedl/.aedl files and print header info")
    p.add_argument("path", nargs="?", default=".", help="directory to scan")
    p.add_argument("--exts", nargs="+", default=[".vedl", ".aedl"], help="file extensions to search")
    args = p.parse_args()

    base = Path(args.path)
    if not base.exists():
        print("Path not found:", base, file=sys.stderr)
        sys.exit(2)

    files = []
    for ext in args.exts:
        files.extend(sorted(base.glob(f"*{ext}")))
        files.extend(sorted(base.glob(f"*{ext.upper()}")))

    if not files:
        print("No .vedl/.aedl files found in", base)
        return

    for fn in files:
        info = probe_file(fn)
        print("="*80)
        print(fn, "-", info['size'], "bytes")
        if info['uint32_0_1']:
            r0, r1, prod = info['uint32_0_1']
            print(f"first8 uint32 (little-endian): {r0}, {r1} -> product {prod}")
        print("MATLAB v5 header:", info['matlab5'])
        print("ASCII-like header:", info['ascii_like'])
        if info['likely_float32_matrix']:
            print("Likely 8-byte header + float32 payload matrix shape:", info['likely_float32_matrix'])
        print("\nFirst 64 bytes (hex):")
        print(hexdump(info['hdr64'], width=16))
        # quick shape guesses (may be noisy)
        shapes = try_load_shapes(fn)
        if shapes:
            print("\nQuick shape guesses from file size:")
            for sdesc, shp in shapes:
                print("  -", sdesc, "->", shp)
        print()

if __name__ == "__main__":
    main()