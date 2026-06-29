"""
Tests for functions in mutational_simulator.py that require file I/O.
Uses tmp_path (synthetic BED/mask files) and real MAF fixtures from tests/fixtures/.
Run with: pytest tests/test_mutational_simulator_files.py -v
"""

import os
import shutil
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from SigProfilerSimulator.mutational_simulator import (
    bed_ranges,
    combine_simulation_files,
    probability_mask,
)

FIXTURES_MAF = os.path.join(os.path.dirname(__file__), "fixtures", "maf_intermediates")

MAF_HEADER = "\t".join([
    "Hugo_symbol", "Entrez_gene_ID", "Center", "Genome", "Chrom",
    "Start_position", "End_position", "Strand", "Variant_Classification",
    "Variant_Type", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2",
    "dbSNP_RS", "dbSNP_Val_Status", "Tumor_Sample_Barcode", "matGenClass",
])


# ─── Fixtures ─────────────────────────────────────────────────────────────────

@pytest.fixture
def bed_file(tmp_path):
    content = "chrom\tstart\tend\n1\t100\t105\n1\t200\t203\n"
    p = tmp_path / "test.bed"
    p.write_text(content)
    return str(p)


@pytest.fixture
def mask_file(tmp_path):
    content = (
        "Chrom\tStart\tEnd\tProbability\n"
        "1\t0\t999\t0.6\n"
        "1\t1000\t1999\t0.4\n"
        "2\t0\t4999\t1.0\n"
    )
    p = tmp_path / "mask.tsv"
    p.write_text(content)
    return str(p)


@pytest.fixture
def maf_output_dir(tmp_path):
    """Copies real intermediate MAF fixtures into a fresh tmp_path."""
    for fname in ("1_1.maf", "1_2.maf", "1_6.maf"):
        shutil.copy(os.path.join(FIXTURES_MAF, fname), tmp_path / fname)
    return str(tmp_path) + "/"


# ─── bed_ranges ───────────────────────────────────────────────────────────────

class TestBedRanges:
    def test_returns_correct_positions(self, bed_file):
        result = bed_ranges("1", bed_file, cushion=0)
        expected = sorted(set(range(100, 105)) | set(range(200, 203)))
        assert sorted(result) == expected

    def test_cushion_extends_range(self, bed_file):
        result = bed_ranges("1", bed_file, cushion=2)
        expected = sorted(set(range(98, 107)) | set(range(198, 205)))
        assert sorted(result) == expected

    def test_absent_chrom_returns_empty(self, bed_file):
        result = bed_ranges("9", bed_file, cushion=0)
        assert result == []

    def test_chr_prefix_stripped(self, tmp_path):
        content = "chrom\tstart\tend\nchr1\t100\t105\n"
        p = tmp_path / "chr.bed"
        p.write_text(content)
        result = bed_ranges("1", str(p), cushion=0)
        assert sorted(result) == list(range(100, 105))

    def test_result_sorted_and_unique(self, bed_file):
        result = bed_ranges("1", bed_file, cushion=0)
        assert result == sorted(set(result))


# ─── probability_mask ─────────────────────────────────────────────────────────

class TestProbabilityMask:
    def test_absent_chrom_returns_none_pair(self, mask_file):
        ranges, probs = probability_mask("9", mask_file)
        assert ranges is None
        assert probs is None

    def test_returns_correct_ranges(self, mask_file):
        ranges, _ = probability_mask("1", mask_file)
        # [[Start, End-Start+1], ...]
        assert ranges == [[0, 1000], [1000, 1000]]

    def test_returns_correct_probs(self, mask_file):
        _, probs = probability_mask("1", mask_file)
        assert probs == [0.6, 0.4]

    def test_single_region_full_prob(self, mask_file):
        ranges, probs = probability_mask("2", mask_file)
        assert ranges == [[0, 5000]]
        assert probs == [1.0]

    def test_invalid_probs_raises(self, tmp_path):
        bad = (
            "Chrom\tStart\tEnd\tProbability\n"
            "1\t0\t999\t0.3\n"
            "1\t1000\t1999\t0.3\n"
        )
        p = tmp_path / "bad_mask.tsv"
        p.write_text(bad)
        with pytest.raises(ValueError):
            probability_mask("1", str(p))


# ─── combine_simulation_files ─────────────────────────────────────────────────

class TestCombineSimulationFiles:
    def test_creates_combined_maf(self, maf_output_dir):
        combine_simulation_files([1], maf_output_dir, ["1", "2", "6"])
        assert os.path.exists(maf_output_dir + "1.maf")

    def test_combined_maf_has_correct_header(self, maf_output_dir):
        combine_simulation_files([1], maf_output_dir, ["1", "2", "6"])
        with open(maf_output_dir + "1.maf") as f:
            first_line = f.readline().strip()
        assert first_line == MAF_HEADER

    def test_combined_maf_has_all_mutations(self, maf_output_dir):
        combine_simulation_files([1], maf_output_dir, ["1", "2", "6"])
        with open(maf_output_dir + "1.maf") as f:
            lines = f.readlines()
        # header + 0 (chrom1 empty) + 2 (chrom2) + 3 (chrom6)
        assert len(lines) == 6

    def test_intermediate_files_deleted(self, maf_output_dir):
        combine_simulation_files([1], maf_output_dir, ["1", "2", "6"])
        assert not os.path.exists(maf_output_dir + "1_1.maf")
        assert not os.path.exists(maf_output_dir + "1_2.maf")
        assert not os.path.exists(maf_output_dir + "1_6.maf")

    def test_missing_chrom_file_does_not_raise(self, maf_output_dir):
        combine_simulation_files([1], maf_output_dir, ["1", "2", "6", "9"])
        assert os.path.exists(maf_output_dir + "1.maf")

    def test_chrom2_mutations_in_output(self, maf_output_dir):
        combine_simulation_files([1], maf_output_dir, ["1", "2", "6"])
        with open(maf_output_dir + "1.maf") as f:
            content = f.read()
        assert "215057995" in content
        assert "101049053" in content

    def test_chrom6_mutations_in_output(self, maf_output_dir):
        combine_simulation_files([1], maf_output_dir, ["1", "2", "6"])
        with open(maf_output_dir + "1.maf") as f:
            content = f.read()
        assert "90639784" in content
        assert "14214055" in content
        assert "26072434" in content
