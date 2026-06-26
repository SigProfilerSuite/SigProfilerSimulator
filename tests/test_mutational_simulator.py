"""
Unit tests for pure functions in SigProfilerSimulator.
No genome reference files required.
Run with: pytest tests/test_mutational_simulator.py -v
"""

import sys
import os
import pytest
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from SigProfilerSimulator.mutational_simulator import (
    revcompl,
    revbias,
    intersections,
    update_chromosome,
    random_base,
    random_base_sub,
    noise,
)
from SigProfilerSimulator.SigProfilerSimulator import context_identifier


# ─── revcompl ─────────────────────────────────────────────────────────────────

class TestRevCompl:
    def test_single_A(self):
        assert revcompl("A") == "T"

    def test_single_T(self):
        assert revcompl("T") == "A"

    def test_single_C(self):
        assert revcompl("C") == "G"

    def test_single_G(self):
        assert revcompl("G") == "C"

    def test_single_N(self):
        assert revcompl("N") == "N"

    def test_palindrome(self):
        assert revcompl("ACGT") == "ACGT"

    def test_triplet(self):
        assert revcompl("ACG") == "CGT"

    def test_empty(self):
        assert revcompl("") == ""


# ─── revbias ──────────────────────────────────────────────────────────────────

class TestRevBias:
    def test_U_maps_to_T(self):
        assert revbias("U") == "T"

    def test_T_maps_to_U(self):
        assert revbias("T") == "U"

    def test_1_maps_to_2(self):
        assert revbias("1") == "2"

    def test_2_maps_to_1(self):
        assert revbias("2") == "1"

    def test_symmetric_B(self):
        assert revbias("B") == "B"

    def test_symmetric_N(self):
        assert revbias("N") == "N"

    def test_symmetric_0(self):
        assert revbias("0") == "0"

    def test_symmetric_3(self):
        assert revbias("3") == "3"

    def test_symmetric_Q(self):
        assert revbias("Q") == "Q"

    def test_reverse_order_UT(self):
        # 'U'->'T', 'T'->'U', reversed → 'UT'
        assert revbias("UT") == "UT"

    def test_reverse_order_12(self):
        # '1'->'2', '2'->'1', reversed → '12'
        assert revbias("12") == "12"


# ─── intersections ────────────────────────────────────────────────────────────

class TestIntersections:
    def test_both_empty(self):
        assert intersections([], []) == []

    def test_a_empty(self):
        assert intersections([], [[1, 5]]) == []

    def test_b_empty(self):
        assert intersections([[1, 5]], []) == []

    def test_no_overlap(self):
        assert intersections([[1, 5]], [[7, 10]]) == []

    def test_partial_overlap(self):
        assert intersections([[1, 10]], [[5, 15]]) == [[5, 10]]

    def test_one_inside_other(self):
        assert intersections([[1, 20]], [[5, 10]]) == [[5, 10]]

    def test_exact_match(self):
        assert intersections([[3, 7]], [[3, 7]]) == [[3, 7]]

    def test_adjacent_merged(self):
        # [1,5]∩[3,7]=[3,5] and [5,10]∩[3,7]=[5,7] → adjacent, merged to [3,7]
        assert intersections([[1, 5], [5, 10]], [[3, 7]]) == [[3, 7]]

    def test_multiple_non_adjacent(self):
        assert intersections([[1, 5], [8, 12]], [[3, 10]]) == [[3, 5], [8, 10]]


# ─── update_chromosome ────────────────────────────────────────────────────────

class TestUpdateChromosome:
    def test_snp_modifies_correct_position(self):
        chrom = [1, 2, 3, 4, 5]
        result = update_chromosome(chrom, 2, "9", "SNP")
        assert result[2] == 9

    def test_snp_leaves_other_positions_unchanged(self):
        chrom = [1, 2, 3, 4, 5]
        update_chromosome(chrom, 2, "9", "SNP")
        assert chrom[0] == 1 and chrom[1] == 2 and chrom[3] == 4 and chrom[4] == 5

    def test_del_writes_at_location(self):
        chrom = [1, 2, 3, 4, 5]
        result = update_chromosome(chrom, 1, "89", "Del")
        assert result[1] == 8
        assert result[2] == 9

    def test_ins_writes_at_location_plus_one(self):
        chrom = [1, 2, 3, 4, 5]
        result = update_chromosome(chrom, 1, "78", "Ins")
        assert result[2] == 7
        assert result[3] == 8

    def test_dinuc_writes_values_directly(self):
        chrom = [10, 20, 30, 40, 50]
        result = update_chromosome(chrom, 1, [99, 88], "DINUC")
        assert result[1] == 99
        assert result[2] == 88

    def test_returns_same_object(self):
        chrom = [1, 2, 3]
        result = update_chromosome(chrom, 0, "5", "SNP")
        assert result is chrom


# ─── random_base ──────────────────────────────────────────────────────────────

class TestRandomBase:
    VALID = set("ATCG")

    def test_full_range_always_valid(self):
        for _ in range(200):
            assert random_base(0, 3) in self.VALID

    def test_restricted_range(self):
        for _ in range(100):
            assert random_base(0, 1) in {"A", "T"}

    def test_single_value(self):
        for _ in range(20):
            assert random_base(2, 2) == "C"


# ─── random_base_sub ──────────────────────────────────────────────────────────

class TestRandomBaseSub:
    def test_returns_from_given_set(self):
        for _ in range(100):
            assert random_base_sub("ACG") in {"A", "C", "G"}

    def test_single_base(self):
        for _ in range(20):
            assert random_base_sub("T") == "T"

    def test_covers_all_chars_in_set(self):
        results = {random_base_sub("ACGT") for _ in range(500)}
        assert results == {"A", "C", "G", "T"}


# ─── noise ────────────────────────────────────────────────────────────────────

class TestNoise:
    def test_no_noise_returns_none(self):
        result = noise({"C>A": 100, "C>G": 50})
        assert result is None

    def test_poisson_returns_same_dict(self):
        s = {"C>A": 100}
        result = noise(s, noisePoisson=True)
        assert result is s

    def test_poisson_values_nonnegative(self):
        np.random.seed(0)
        s = {"C>A": 100, "C>G": 50, "C>T": 200}
        result = noise(s, noisePoisson=True)
        for val in result.values():
            assert val >= 0

    def test_poisson_mean_close_to_lambda(self):
        np.random.seed(42)
        draws = []
        for _ in range(1000):
            s = {"C>A": 100}
            noise(s, noisePoisson=True)
            draws.append(s["C>A"])
        assert abs(np.mean(draws) - 100) < 5

    def test_uniform_returns_same_dict(self):
        s = {"C>A": 1000}
        result = noise(s, noiseUniform=20)
        assert result is s

    def test_uniform_within_expected_bounds(self):
        np.random.seed(42)
        s = {"C>A": 1000}
        noise(s, noiseUniform=20)
        # ±10% → [900, 1100]
        assert 850 <= s["C>A"] <= 1150

    def test_uniform_never_negative(self):
        np.random.seed(0)
        for _ in range(200):
            s = {"C>A": 1}
            noise(s, noiseUniform=500)
            assert s["C>A"] >= 0


# ─── context_identifier ───────────────────────────────────────────────────────

class TestContextIdentifier:
    def test_sbs6(self):
        ctx, nuc = context_identifier("C>A")
        assert ctx == "6"
        assert nuc == "C"

    def test_sbs96(self):
        ctx, nuc = context_identifier("A[C>A]T")
        assert ctx == "96"
        assert nuc == "ACT"

    def test_sbs1536(self):
        ctx, nuc = context_identifier("AA[C>A]TT")
        assert ctx == "1536"
        assert nuc == "AACTT"

    def test_sbs384(self):
        ctx, nuc = context_identifier("T:A[C>A]T")
        assert ctx == "384"
        assert nuc == "T:ACT"

    def test_sbs6144(self):
        ctx, nuc = context_identifier("AAA[C>A]TTT")
        assert ctx == "6144"
        assert nuc == "AAACTTT"

    def test_dbs78(self):
        ctx, nuc = context_identifier("AC>GT")
        assert ctx == "DBS78"
        assert nuc == "AC"

    def test_dbs186(self):
        ctx, nuc = context_identifier("T:AC>GT")
        assert ctx == "DBS186"
        assert nuc == "T:AC"

    def test_invalid_mutation_exits(self):
        with pytest.raises(SystemExit):
            context_identifier("INVALID")
