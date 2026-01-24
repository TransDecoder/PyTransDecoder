#!/usr/bin/env bash
# Wrapper for gff3_file_to_bed.py
exec python3 "$(dirname "$0")/gff3_file_to_bed.py" "$@"
