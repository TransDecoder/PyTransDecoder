#!/usr/bin/env bash
# Wrapper for fasta_prot_checker.py
exec python3 "$(dirname "$0")/fasta_prot_checker.py" "$@"
