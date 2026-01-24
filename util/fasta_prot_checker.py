#!/usr/bin/env python3
"""
Validate protein FASTA file

Checks for:
- Starts with M (methionine)
- Ends with * (stop codon)
- No internal stop codons

Python port of fasta_prot_checker.pl
"""

import sys
from Bio import SeqIO


def check_protein_fasta(fasta_file):
    """
    Validate protein sequences in FASTA file
    
    Args:
        fasta_file: Path to protein FASTA file
        
    Returns:
        Number of sequences with intervening stops
    """
    seen_intervening_stop = 0
    
    with open(fasta_file) as f:
        for record in SeqIO.parse(f, "fasta"):
            header = record.id
            sequence = str(record.seq).upper()
            
            errors = []
            
            # Check for trailing stop codon
            if not sequence.endswith('*'):
                errors.append("No stop codon")
            
            # Count internal stops (excluding trailing)
            seq_no_trailing_stop = sequence.rstrip('*')
            num_stops = seq_no_trailing_stop.count('*')
            
            if num_stops > 0 or len(seq_no_trailing_stop) == 0:
                errors.append(f"*{num_stops}")
                seen_intervening_stop += 1
            
            # Check for start codon
            if not sequence.startswith('M'):
                errors.append("Doesn't start with M")
            
            # Check for trailing stop codon
            if not sequence.endswith('*'):
                errors.append("No stop codon")
            
            # Report errors (only if there are internal stops)
            if errors and num_stops > 0:
                print(f"{header}\tERRORS: {', '.join(errors)}", file=sys.stderr)
    
    if seen_intervening_stop > 0:
        print(f"\nError, found {seen_intervening_stop} proteins with intervening stop codons.",
              file=sys.stderr)
        sys.exit(1)
    
    return seen_intervening_stop


def main():
    if len(sys.argv) < 2:
        print("usage: fasta_prot_checker.py proteins.fa", file=sys.stderr)
        sys.exit(1)
    
    fasta_file = sys.argv[1]
    check_protein_fasta(fasta_file)
    print(f"OK - {fasta_file} passed validation")


if __name__ == '__main__':
    main()
