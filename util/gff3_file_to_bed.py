#!/usr/bin/env python3
"""
Convert GFF3 file to BED format

Python port of gff3_file_to_bed.pl
"""

import sys


def gff3_to_bed(gff3_file, output_file=None):
    """
    Convert GFF3 to BED format
    
    Args:
        gff3_file: Input GFF3 file
        output_file: Output BED file (default: stdout)
    """
    out = open(output_file, 'w') if output_file else sys.stdout
    
    try:
        with open(gff3_file) as f:
            for line in f:
                line = line.strip()
                
                # Skip comments and empty lines
                if not line or line.startswith('#'):
                    continue
                
                parts = line.split('\t')
                if len(parts) < 9:
                    continue
                
                # GFF3 format: seqid, source, type, start, end, score, strand, phase, attributes
                seqid = parts[0]
                start = int(parts[3]) - 1  # Convert to 0-based
                end = int(parts[4])  # Already 1-based end (inclusive)
                score = parts[5] if parts[5] != '.' else '0'
                strand = parts[6]
                attributes = parts[8]
                
                # Extract ID from attributes
                name = seqid
                for attr in attributes.split(';'):
                    if '=' in attr:
                        key, value = attr.split('=', 1)
                        if key == 'ID':
                            name = value
                            break
                
                # BED format: chrom, chromStart, chromEnd, name, score, strand
                out.write(f"{seqid}\t{start}\t{end}\t{name}\t{score}\t{strand}\n")
    
    finally:
        if output_file:
            out.close()


def main():
    if len(sys.argv) < 2:
        print("usage: gff3_file_to_bed.py file.gff3 [> output.bed]", file=sys.stderr)
        sys.exit(1)
    
    gff3_file = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else None
    
    gff3_to_bed(gff3_file, output_file)


if __name__ == '__main__':
    main()
