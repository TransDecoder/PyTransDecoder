# Porting notes (historical)

These notes were written while porting TransDecoder from Perl to Python (January to March 2026). They describe
the original port, which reproduced TransDecoder v6.0.0, including its 100-aa minimum, position-weight-matrix
start refinement and multiple ORFs per transcript. They do not describe PyTransDecoder v0.2.0 or later. For
current behaviour see the main [README](../../README.md), [CHANGELOG](../../CHANGELOG.md) and
[start codon refinement](../start_codon_refinement.md).
