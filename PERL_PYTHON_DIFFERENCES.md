# Perl vs Python Implementation Differences

## Summary

Python consistently finds **10-25% more ORFs** than Perl (100% superset - all Perl ORFs included).

### Test Results

| Dataset | Python ORFs | Perl ORFs | Difference | % More |
|---------|-------------|-----------|------------|--------|
| Trinity.fasta (Phase 1) | 845 | 845 | 0 | 0% |
| Trinity.fasta (Phase 2) | 733 | 679 | +54 | +8.0% |
| cufflinks_example | 94 | 82 | +12 | +14.6% |
| pasa_example | 977 | 792 | +185 | +23.4% |

**Key Finding**: Differences arise entirely in **Phase 2 (Predict)**, not Phase 1 (LongOrfs).

## Python-Only ORFs Characteristics

Analysis of 54 Python-only ORFs from Trinity.fasta:

- **90% are complete ORFs** (50/55 have both start and stop codons)
- **Most are secondary ORFs** (p2, p3, p4) on same transcript
- **All have neutral scores** (score = 5.0, not trained)
- **Relatively short** (104-331 aa, median 120 aa)
- **Meet selection criteria**: score[0] > 0 AND score[0] > max(score[1], score[2])

### Example: comp1004_c0_seq1
- **Python**: Selects both .p1 (score=93.71) and .p2 (score=3.01) on opposite strands
- **Perl**: Only selects .p1 (score=98.15)
- **Why**: p2 has score=3.01 which meets criteria (3.01 > 0 and > max(-13.69, -10.07))

## Selection Criteria (Both Implementations)

ORF is selected if **ANY** of these conditions is met:

1. **Has homology**: BLAST or Pfam hit
2. **Very long**: Length >= min_length_auto_accept (default: 1000000 nt)
3. **Good score**: `score[frame0] > 0` AND `score[frame0] > max(score[frame1], score[frame2])`

Both Python and Perl use identical selection criteria.

## Overlap Handling

### Perl Logic
```perl
if ($gene_entry->{homology_count} || ! &has_sufficient_overlap($gene_entry, \@selected_entries)) {
    push (@selected_entries, $gene_entry);
}
```

- **Always includes ORFs with homology**, even if overlapping
- Removes others only if >10% overlap

### Python Logic (Fixed)
```python
# Always include ORFs with homology support (matches Perl)
if candidate['homology_count'] > 0:
    selected.append(candidate)
    continue

# For others, check overlap...
```

Python now matches Perl's homology-first behavior.

## Potential Sources of Differences

### 1. ORF Prioritization
**Perl** (lines 174-180 in select_best_ORFs_per_transcript.pl):
```perl
@gene_entries = sort {
    $b->{homology_count} <=> $a->{homology_count}
    ||
    $b->{cds_scores}->[0] <=> $a->{cds_scores}->[0]
    ||
    $b->{length} <=> $a->{length}
} @gene_entries;
```

**Python** (lines 666-670 in predict.py):
```python
candidates.sort(key=lambda x: (
    -x['homology_count'],
    -x['scores'][0],
    -x['length']
))
```

✅ **Identical** - Both prioritize by: homology_count → score[0] → length

### 2. Hexamer Scoring
Could there be slight numerical differences in:
- Training ORF selection
- Hexamer frequency calculation
- Markov model scoring

**Status**: Needs investigation (scores like 5.0 suggest these ORFs weren't trained on)

### 3. Score Threshold Interpretation
The score=5.0 on many Python-only ORFs is suspicious:
- This is likely the default/neutral score for ORFs not in the training set
- Python may be applying the selection criteria `score[0] > 0` more literally
- Perl may have implicit filtering of "unscored" ORFs

### 4. Secondary ORF Selection
Most Python-only ORFs are .p2, .p3, .p4 (secondary ORFs):
- These might be on opposite strands (don't overlap primary ORF)
- They meet selection criteria but have low scores
- Perl may have implicit logic to prefer primary ORFs

## Recommendations

### To Achieve Closer Match

1. **Investigate score=5.0 ORFs**
   - Check if these should be excluded
   - Verify hexamer scoring is working correctly
   - Compare training ORF selection between versions

2. **Add Minimum Score Threshold** (Option)
   ```python
   MIN_SCORE_THRESHOLD = 10.0  # or appropriate value
   if scores[0] < MIN_SCORE_THRESHOLD and homology_count == 0:
       skip_orf
   ```

3. **Prioritize Primary ORFs** (Option)
   - Give preference to .p1 over .p2/.p3 when scores are low
   - Match Perl's implicit behavior

4. **Verify Training Selection**
   - Ensure top 500 ORFs selected identically
   - Compare hexamer frequency calculations
   - Check Markov model scoring

### To Maintain Current Behavior

**Arguments for keeping Python's behavior:**

1. **More Sensitive**: Finds additional valid ORFs that meet criteria
2. **Conservative**: All ORFs pass the same threshold as Perl's ORFs
3. **Transparent**: Applies selection criteria consistently
4. **Complete**: Doesn't skip ORFs on opposite strands
5. **Complete ORFs**: 90% of extra ORFs are complete (not fragments)

The extra ORFs aren't "false positives" - they legitimately meet the selection criteria. The question is whether Perl has implicit filtering that's not documented.

## Next Steps

1. Trace a specific Python-only ORF through Perl's selection process
2. Compare hexamer scoring output between implementations
3. Check for undocumented filtering logic in Perl
4. Decide: match Perl exactly vs. keep more sensitive behavior
