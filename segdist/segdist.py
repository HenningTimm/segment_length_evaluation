#!/usr/bin/env python3
"""Compute the expected segment length distribution of compressed winnowing windows for
a given hash function codomain size C and winnowing windows size w for segments up to length k.

This code was developed in cooperation with Dr. Denis Kurz (https://ls11-www.cs.tu-dortmund.de/staff/kurz).
"""
import argparse
import sys


class UnboundedSegmentLengthDistribution:
    """Compute the probabilities \overline{p}_{k}:
    What is the probability for a segment starting at or before X_1, and ending at or after X_k?
    """

    def __init__(self, C, w):
        """Initialize an empty list of probabilities.
        """
        self._C = C
        self._w = w
        self._c_k_m = []

    def _grow_c_k_m(self):
        """Add another row to the _c_k_m array:
           m =  0 .. C
          k=0
           ..  p_{m,k}
        max_k

        In the thesis, this is the array P.
        Each entry \overline{p}_{m, k} contains the probability to
        observe segments of length k with a specific hash value m.
        """
        c_k = []
        k = len(self._c_k_m) + 1

        for m in range(1, self._C + 1):
            N = self._C - m
            p_first_min_before_xk = sum([(N ** (j - 1)) * self.count_length_at_least_k_with_m(k - j, m)
                            for j in range(1, min(self._w, k - 1) + 1)])
            p_first_min_in_overlap = 0
            if k <= self._w:
                p_first_min_in_overlap = (N ** (k - 1)) * ((N + 1) ** self._w)
                p_first_min_in_overlap -= (N ** self._w) * ((N + 1) ** (k - 1))
            c_k.append(p_first_min_before_xk + p_first_min_in_overlap)
        self._c_k_m.append(c_k)

    def count_length_at_least_k_with_m(self, k, m):
        """Offset for one-based indexing.
        E.g. the values for `k` are stored at `k-1`.
        """
        return self._c_k_m[k - 1][m - 1]

    def length_at_least(self, k):
        """Provide the probability \overline{p}_{k}.
        Used to compute the left bounded segment length distribution.
        """
        # Dynamically extend the _c_k_m array as required
        # for this values of k.
        while len(self._c_k_m) < k:
            self._grow_c_k_m()
        return sum([self.count_length_at_least_k_with_m(k, m) / (self._C ** (k + self._w - 1))
                        for m in range(1, self._C + 1)])

    # def length_at_least_m(self, k, n):
    #     return self.count_length_at_least_k_with_m(k, n) / (m ** (k + self._w - 1))


class LeftBoundedSegmentLengthDistribution:
    """Compute the probabilities \overline{p}_{k}^l:
    What is the probability for a segment starting at or before X_1, and ending at X_k?
    """

    def __init__(self, unbounded):
        """Initialize with an UnboundedSegmentLengthDistribution object.
        """
        self._p_unbounded = unbounded

    def length_at_least(self, k):
        """Compute the probability \overline{p}_{k}^l that a left-bounded
        segment has a length of at least k.

        Used to compute the left and right bounded segment length distribution.
        """
        p_k = self._p_unbounded.length_at_least(k)
        p_k_plus_one = self._p_unbounded.length_at_least(k + 1)
        return p_k - p_k_plus_one


class LeftRightSegmentLengthDistribution:
    """Compute the probabilities p_{k}^{lr}:
    What is the probability for a segment starting at X_1, and ending at X_k?
    """
    def __init__(self, left_bounded):
        self._p_left_bounded = left_bounded

    def length_exactly(self, k):
        p_k = self._p_left_bounded.length_at_least(k)
        p_k_plus_one = self._p_left_bounded.length_at_least(k + 1)
        return p_k - p_k_plus_one


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compute probabilities for the number of adjacent winnowing window positions with a common minimum.\n\nValues psi_{C, k})_{i=k}^{K} are written to stdout, the error value psi_{C, >K} is written to stderr."
    )
    parser.add_argument('--universe-size', '-C', dest='C', type=int)
    parser.add_argument('--window-size', '-w', dest='w', type=int)
    parser.add_argument('--max-sequence_length', '-K', dest='max_k', type=int)
    args = parser.parse_args()
    return (args.C, args.w, args.max_k)


if __name__ == '__main__':
    C, w, max_k = parse_args()

    # The goal is to compute the sequence $\Psi_{w, K}^{C} = (psi_{C, k})_{i=k}^{K} \cup psi_{C, >K}$
    # which can be computed from the probability p_{k}^{lr} (cf. Equation 6.11) normalized by the
    # probability to observe left bounded segments at any specific position.
    P = UnboundedSegmentLengthDistribution(C, w)
    Pl = LeftBoundedSegmentLengthDistribution(P)
    Plr = LeftRightSegmentLengthDistribution(Pl)

    # Probability to observe any left-bounded segment.
    # This is used to 'pin down' a starting point for
    # segments without picking a specific position for X_1.
    pl_1 = Pl.length_at_least(1)

    # Compute probability psi_{C, k} to observe segments of length k
    for k in range(1, max_k + 1):
        p = P.length_at_least(k)
        pl = Pl.length_at_least(k)
        plr = Plr.length_exactly(k)
        
        # Probability that any segment with length exactly k is observed (i.e. a
        # segment that is both left and right bounded), for the probability
        # that a left bounded segment exists at any specific point.
        psi_C_k = plr / pl_1
        print(psi_C_k)
    # Compute the error term psi_{C, >K} for all remaining segments that are cut off by K (max_k)
    psi_C_greater_K = Pl.length_at_least(max_k + 1) / pl_1
    print(psi_C_greater_K, file=sys.stderr)
