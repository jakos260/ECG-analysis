#!/usr/bin/env python3
"""Rename files in a directory by stripping the directory-name prefix.

Example:
  Directory: A/B/C contains files `C_file1.txt`, `C_file2.txt`.
  Run: python rename_strip_dir_prefix.py A/B/C
  Result: files renamed to `file1.txt`, `file2.txt`.

Options:
  --prefix PREFIX   Use explicit prefix instead of directory name.
  --dry-run         Print what would be done without renaming.
  --overwrite       Overwrite existing target files if present.
"""
from pathlib import Path
import argparse
import sys


def rename_strip_dir_prefix(target_dir: Path, prefix: str, dry_run: bool, overwrite: bool) -> int:
    prefix_with_sep = prefix + "_"
    renamed = 0
    skipped = 0
    errors = 0

    for p in sorted(target_dir.iterdir()):
        if not p.is_file():
            continue
        name = p.name
        if not name.startswith(prefix_with_sep):
            continue

        new_name = name[len(prefix_with_sep):]
        dest = target_dir / new_name

        try:
            if dest.exists():
                if overwrite:
                    if dry_run:
                        print(f"DRY RUN: overwrite {p.name} -> {new_name}")
                    else:
                        dest.unlink()
                        p.rename(dest)
                        print(f"OVERWRITTEN: {p.name} -> {new_name}")
                    renamed += 1
                else:
                    print(f"SKIP (exists): {p.name} -> {new_name}")
                    skipped += 1
            else:
                if dry_run:
                    print(f"DRY RUN: {p.name} -> {new_name}")
                else:
                    p.rename(dest)
                    print(f"RENAMED: {p.name} -> {new_name}")
                renamed += 1
        except Exception as e:
            print(f"ERROR renaming {p} -> {dest}: {e}", file=sys.stderr)
            errors += 1

    print(f"Done. Renamed={renamed}, Skipped={skipped}, Errors={errors}")
    return 0 if errors == 0 else 2


def main(argv=None):
    ap = argparse.ArgumentParser(description="Strip directory-name prefix from files in that directory")
    ap.add_argument("dir", help="Target directory (e.g. A/B/C)")
    ap.add_argument("--prefix", help="Prefix to remove (default: directory name)")
    ap.add_argument("--dry-run", action="store_true", help="Don't actually rename, just show actions")
    ap.add_argument("--overwrite", action="store_true", help="Overwrite existing target files")
    args = ap.parse_args(argv)

    target = Path(args.dir)
    if not target.exists() or not target.is_dir():
        print(f"Error: '{target}' is not a directory", file=sys.stderr)
        return 2

    prefix = args.prefix if args.prefix is not None else target.name
    return rename_strip_dir_prefix(target, prefix, args.dry_run, args.overwrite)


if __name__ == "__main__":
    sys.exit(main())
