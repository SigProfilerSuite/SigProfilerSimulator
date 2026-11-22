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
    simulator,
)

# Mirrors the tsb_ref / tsb_ref_rev encoding built in SigProfilerSimulator.py:
# byte value -> [transcriptional bias, nucleotide]. simulator() takes these as
# plain arguments rather than importing them, so tests must supply the same mapping.
TSB_REF = {0: ['N', 'A'], 1: ['N', 'C'], 2: ['N', 'G'], 3: ['N', 'T'],
           4: ['T', 'A'], 5: ['T', 'C'], 6: ['T', 'G'], 7: ['T', 'T'],
           8: ['U', 'A'], 9: ['U', 'C'], 10: ['U', 'G'], 11: ['U', 'T'],
           12: ['B', 'A'], 13: ['B', 'C'], 14: ['B', 'G'], 15: ['B', 'T'],
           16: ['N', 'N'], 17: ['T', 'N'], 18: ['U', 'N'], 19: ['B', 'N']}
TSB_REF_REV = {'N': {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 16},
               'T': {'A': 4, 'C': 5, 'G': 6, 'T': 7, 'N': 17},
               'U': {'A': 8, 'C': 9, 'G': 10, 'T': 11, 'N': 18},
               'B': {'A': 12, 'C': 13, 'G': 14, 'T': 15, 'N': 19}}

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


# ─── simulator: indel coordinate writing (audit sigprofilersimulator-2026-05 #35) ─

@pytest.fixture
def synthetic_chrom_dir(tmp_path):
    """A single synthetic chromosome '1': alternating C/T bytes so that any
    position landed on has mismatching neighbors (repeat_count is always 0),
    letting the test stay deterministic without depending on fastrand's
    internal sequence."""
    d = tmp_path / "chroms"
    d.mkdir()
    pattern = bytes([1, 3] * 40)  # 'N'/'C', 'N'/'T' repeated -> "CTCTCT..."
    (d / "1.txt").write_bytes(pattern)
    return str(d) + "/"


def _run_indel_simulator(synthetic_chrom_dir, tmp_path, vcf, label):
    sample = "S1"
    chrom = "1"
    mutation_tracker = {"ID": {sample: {"1:Del:C:0": {chrom: 1}}}}
    output_path = str(tmp_path / ("out_" + label)) + "/"
    if vcf:
        os.makedirs(output_path + sample + "/")
    log_file = str(tmp_path / ("log_" + label + ".txt"))

    simulator(
        sample_names=[sample],
        mutation_tracker=mutation_tracker,
        chromosome_string_path=synthetic_chrom_dir,
        tsb_ref=TSB_REF,
        tsb_ref_rev=TSB_REF_REV,
        simulation_number=1,
        seed=7,
        cushion=0,
        output_path=output_path,
        updating=False,
        chromosomes=[chrom],
        project="f35demo",
        genome="GRCh37",
        bed=False,
        bed_file=None,
        contexts=["ID"],
        overlap=False,
        project_path=str(tmp_path) + "/",
        seqInfo=False,
        log_file=log_file,
        spacing=1,
        noisePoisson=False,
        noiseUniform=0,
        vcf=vcf,
        mask=None,
    )

    if vcf:
        out_file = output_path + sample + "/" + sample + "_1_" + chrom + ".vcf"
    else:
        out_file = output_path + "1_" + chrom + ".maf"
    with open(out_file) as f:
        rows = [l.rstrip("\n").split("\t") for l in f if l.strip()]
    assert len(rows) == 1, "expected exactly one simulated deletion, got %r" % (rows,)
    return rows[0]


class TestSimulatorIndelCoordinates:
    def test_maf_start_position_matches_vcf_position(self, synthetic_chrom_dir, tmp_path):
        """Finding #35: default (MAF) mode must anchor Reference_Allele at the
        same genomic position that vcf=True mode reports, not one base downstream."""
        vcf_row = _run_indel_simulator(synthetic_chrom_dir, tmp_path, vcf=True, label="vcf")
        maf_row = _run_indel_simulator(synthetic_chrom_dir, tmp_path, vcf=False, label="maf")

        vcf_pos = int(vcf_row[1])
        maf_start_position = int(maf_row[5])

        assert maf_start_position == vcf_pos
