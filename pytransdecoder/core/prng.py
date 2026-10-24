"""
Portable pseudo-random numbers for decoy shuffles, identical in PyTransDecoder and TransDecoder (Perl).

xoshiro128** (Blackman & Vigna 2018) on 32-bit words, seeded from FNV-1a hashes of a seed string.
All arithmetic stays below 2^53, so the Perl implementation (PerlLib/Decoys.pm) gives the same stream.
"""
from typing import List

M32 = 0xFFFFFFFF


def fnv1a32(text: str) -> int:
    h = 0x811C9DC5
    for b in text.encode():
        h = ((h ^ b) * 0x01000193) & M32
    return h


def _rotl(x: int, k: int) -> int:
    return ((x << k) | (x >> (32 - k))) & M32


class Xoshiro128ss:
    def __init__(self, seed: str):
        s = [fnv1a32(f"{seed}#{i}") for i in range(4)]
        if not any(s):
            s[0] = 1
        self.s: List[int] = s

    def next32(self) -> int:
        s = self.s
        result = (_rotl((s[1] * 5) & M32, 7) * 9) & M32
        t = (s[1] << 9) & M32
        s[2] ^= s[0]
        s[3] ^= s[1]
        s[1] ^= s[2]
        s[0] ^= s[3]
        s[2] ^= t
        s[3] = _rotl(s[3], 11)
        return result

    def below(self, n: int) -> int:
        """Uniform integer in [0, n) by rejection (n >= 1)."""
        limit = (2 ** 32 // n) * n
        while True:
            r = self.next32()
            if r < limit:
                return r % n

    def shuffle(self, items: list) -> None:
        """Fisher-Yates, in place."""
        for i in range(len(items) - 1, 0, -1):
            j = self.below(i + 1)
            items[i], items[j] = items[j], items[i]
