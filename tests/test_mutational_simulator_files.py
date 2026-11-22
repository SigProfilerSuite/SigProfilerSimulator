"""
Tests for functions in mutational_simulator.py that require file I/O.
Uses tmp_path (synthetic BED/mask files) and real MAF fixtures from tests/fixtures/.
Run with: pytest tests/test_mutational_simulator_files.py -v
"""

import json
import os
import pickle
import shutil
import subprocess
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from SigProfilerSimulator.mutational_simulator import (
    bed_ranges,
    chrom_proportions,
    chrom_proportions_BED,
    combine_simulation_files,
    mut_tracker,
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
    """Standard headerless BED, as produced by real tools (audit finding #8)."""
    content = "1\t100\t105\n1\t200\t203\n"
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

    def test_cushion_excludes_edges(self, tmp_path):
        """Finding #9: cushion excludes edge bases, it does not expand the range."""
        p = tmp_path / "cushion.bed"
        p.write_text("1\t100\t110\n")
        result = bed_ranges("1", str(p), cushion=2)
        assert sorted(result) == list(range(102, 108))

    def test_cushion_can_empty_a_short_interval(self, tmp_path):
        p = tmp_path / "cushion_short.bed"
        p.write_text("1\t100\t103\n")
        result = bed_ranges("1", str(p), cushion=2)
        assert result == []

    def test_absent_chrom_returns_empty(self, bed_file):
        result = bed_ranges("9", bed_file, cushion=0)
        assert result == []

    def test_chr_prefix_stripped(self, tmp_path):
        p = tmp_path / "chr.bed"
        p.write_text("chr1\t100\t105\n")
        result = bed_ranges("1", str(p), cushion=0)
        assert sorted(result) == list(range(100, 105))

    def test_result_sorted_and_unique(self, bed_file):
        result = bed_ranges("1", bed_file, cushion=0)
        assert result == sorted(set(result))

    def test_headerless_bed_is_not_truncated(self, tmp_path):
        """Finding #8: a leading data row must not be mistaken for a header."""
        p = tmp_path / "headerless.bed"
        p.write_text("1\t100\t103\n")
        result = bed_ranges("1", str(p), cushion=0)
        assert sorted(result) == [100, 101, 102]

    def test_blank_and_pragma_lines_are_skipped(self, tmp_path):
        p = tmp_path / "pragma.bed"
        p.write_text("track name=demo\n\n#comment\n1\t100\t103\n")
        result = bed_ranges("1", str(p), cushion=0)
        assert sorted(result) == [100, 101, 102]

    def test_target_chromosome_not_skipped_while_scanning(self, tmp_path):
        """Finding #8: searching past several non-matching rows must not skip
        rows via a double-advanced file iterator."""
        p = tmp_path / "scan.bed"
        p.write_text("2\t0\t5\n3\t0\t5\n4\t0\t5\n5\t0\t5\n1\t100\t103\n")
        result = bed_ranges("1", str(p), cushion=0)
        assert sorted(result) == [100, 101, 102]

    def test_non_contiguous_chromosome_rows_all_collected(self, tmp_path):
        """Finding #8: rows for the target chromosome need not be grouped
        together; a later unrelated chromosome must not stop the scan."""
        p = tmp_path / "interleaved.bed"
        p.write_text("1\t100\t103\n2\t0\t5\n1\t200\t203\n")
        result = bed_ranges("1", str(p), cushion=0)
        expected = sorted(set(range(100, 103)) | set(range(200, 203)))
        assert sorted(result) == expected


# ─── chrom_proportions_BED ────────────────────────────────────────────────────

class TestChromProportionsBED:
    def test_proportions_match_eligible_base_counts(self, tmp_path):
        """Finding #6: proportions must reflect each chromosome's share of the
        total eligible bases (reproduces the audit's 200/50 example)."""
        bed = tmp_path / "two_chrom.bed"
        bed.write_text("1\t0\t200\n2\t0\t50\n")
        chrom_dir = tmp_path / "chroms"
        chrom_dir.mkdir()

        chrom_proportions_BED(str(bed), str(chrom_dir) + "/", "GRCh37", ["1", "2"])

        with open(str(chrom_dir) + "/BED_GRCh37_proportions.txt", "rb") as f:
            probs = pickle.load(f)
        assert probs == pytest.approx([0.8, 0.2])

    def test_proportions_sum_to_one(self, tmp_path):
        bed = tmp_path / "three_chrom.bed"
        bed.write_text("1\t0\t30\n2\t0\t20\n3\t0\t50\n")
        chrom_dir = tmp_path / "chroms"
        chrom_dir.mkdir()

        chrom_proportions_BED(str(bed), str(chrom_dir) + "/", "GRCh37", ["1", "2", "3"])

        with open(str(chrom_dir) + "/BED_GRCh37_proportions.txt", "rb") as f:
            probs = pickle.load(f)
        assert sum(probs) == pytest.approx(1.0)

    def test_chr_prefixed_names_do_not_raise(self, tmp_path):
        """Finding #7: chr-prefixed names must not cause a KeyError."""
        bed = tmp_path / "chr_prefixed.bed"
        bed.write_text("chr1\t0\t100\nchr2\t0\t100\n")
        chrom_dir = tmp_path / "chroms"
        chrom_dir.mkdir()

        chrom_proportions_BED(str(bed), str(chrom_dir) + "/", "GRCh37", ["1", "2"])

        with open(str(chrom_dir) + "/BED_GRCh37_proportions.txt", "rb") as f:
            probs = pickle.load(f)
        assert probs == pytest.approx([0.5, 0.5])

    def test_multiple_rows_per_chromosome_accumulate(self, tmp_path):
        bed = tmp_path / "multi_row.bed"
        bed.write_text("1\t0\t50\n1\t100\t150\n2\t0\t100\n")
        chrom_dir = tmp_path / "chroms"
        chrom_dir.mkdir()

        chrom_proportions_BED(str(bed), str(chrom_dir) + "/", "GRCh37", ["1", "2"])

        with open(str(chrom_dir) + "/BED_GRCh37_proportions.txt", "rb") as f:
            probs = pickle.load(f)
        assert probs == pytest.approx([0.5, 0.5])


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


# ─── sample ordering determinism (audit sigprofilersimulator-2026-05 #4) ──────

SAMPLE_ORDER = ["PD1001a", "PD1002a", "PD1003a", "PD1004a", "PD1005a"]
NUC = "A[C>A]A"


@pytest.fixture
def catalogue_file(tmp_path):
    p = tmp_path / "catalogue.txt"
    header = "MutationType\t" + "\t".join(SAMPLE_ORDER) + "\n"
    row = NUC + "\t" + "\t".join(["1"] * len(SAMPLE_ORDER)) + "\n"
    p.write_text(header + row)
    return str(p)


def _sample_names_under_hashseed(seed, func_name, catalogue_file):
    script = (
        "import json, sys; "
        "sys.path.insert(0, {repo!r}); "
        "from SigProfilerSimulator.mutational_simulator import {func} as fn; "
        "catalogue_files = {{'96': {path!r}}}; "
        "log_file = {log!r}; "
        "{call}"
        "print(json.dumps(sample_names))"
    ).format(
        repo=os.path.join(os.path.dirname(__file__), ".."),
        func=func_name,
        path=catalogue_file,
        log=catalogue_file + ".log",
        call={
            "mutation_preparation": "sample_names, samples = fn(catalogue_files, log_file); ",
            "mutation_preparation_chromosomes": "sample_names, samples, tracker = fn(catalogue_files, '', ['1'], 'proj', log_file); ",
            "mutation_preparation_region": "sample_names, samples, tracker = fn(catalogue_files, '', 'proj', log_file, '1'); ",
        }[func_name],
    )
    env = dict(os.environ)
    env["PYTHONHASHSEED"] = str(seed)
    out = subprocess.run(
        [sys.executable, "-c", script], env=env, capture_output=True, text=True, check=True,
    )
    return json.loads(out.stdout.strip().splitlines()[-1])


@pytest.mark.parametrize("func_name", [
    "mutation_preparation",
    "mutation_preparation_chromosomes",
    "mutation_preparation_region",
])
class TestSampleOrderDeterminism:
    def test_order_independent_of_hashseed(self, catalogue_file, func_name):
        """Finding #4: sample order must not depend on PYTHONHASHSEED."""
        results = [_sample_names_under_hashseed(seed, func_name, catalogue_file) for seed in (0, 1, 2, 3)]
        assert all(r == results[0] for r in results), results

    def test_order_matches_file(self, catalogue_file, func_name):
        result = _sample_names_under_hashseed(0, func_name, catalogue_file)
        assert result == SAMPLE_ORDER


# ─── spacing / overlap enforcement (audit sigprofilersimulator-2026-05 #12, #13) ──

def _write_chrom(tmp_path, pattern_bytes):
    d = tmp_path / "chroms2"
    if not d.exists():
        d.mkdir()
    (d / "1.txt").write_bytes(pattern_bytes)
    return str(d) + "/"


def _run_sim(chrom_dir, tmp_path, mutation_tracker, contexts, spacing, label):
    sample = "S1"
    chrom = "1"
    output_path = str(tmp_path / ("sim_out_" + label)) + "/"
    log_file = str(tmp_path / ("sim_log_" + label + ".txt"))

    simulator(
        sample_names=[sample],
        mutation_tracker=mutation_tracker,
        chromosome_string_path=chrom_dir,
        tsb_ref=TSB_REF,
        tsb_ref_rev=TSB_REF_REV,
        simulation_number=1,
        seed=7,
        cushion=0,
        output_path=output_path,
        updating=False,
        chromosomes=[chrom],
        project="spacingdemo",
        genome="GRCh37",
        bed=False,
        bed_file=None,
        contexts=contexts,
        overlap=False,
        project_path=str(tmp_path) + "/",
        seqInfo=False,
        log_file=log_file,
        spacing=spacing,
        noisePoisson=False,
        noiseUniform=0,
        vcf=False,
        mask=None,
    )

    out_file = output_path + "1_" + chrom + ".maf"
    rows = []
    if os.path.exists(out_file):
        with open(out_file) as f:
            rows = [l.rstrip("\n").split("\t") for l in f if l.strip()]
    return rows


class TestSpacingAndOverlapEnforcement:
    def test_min_spacing_is_symmetric(self, tmp_path):
        """Finding #12: no two placed SBS6 events may end up exactly `spacing`
        apart, regardless of which was placed first."""
        n_pairs = 50
        pattern = bytes([1, 1, 3] * n_pairs)  # "CCT" repeated: adjacent C,C pairs separated by T
        chrom_dir = _write_chrom(tmp_path, pattern)
        mutation_tracker = {"6": {"S1": {"C>A": {"1": n_pairs}}}}

        rows = _run_sim(chrom_dir, tmp_path, mutation_tracker, ["6"], spacing=1, label="spacing")
        positions = sorted(int(r[5]) for r in rows if r[9] == "SNP")

        assert len(positions) == n_pairs
        for a, b in zip(positions, positions[1:]):
            assert b - a > 1, "two SBS6 events placed only %d apart: %r" % (b - a, positions)

    def test_dbs_second_base_is_reserved(self, tmp_path):
        """Finding #13: a DBS's second base must be reserved so a later event
        (here an SBS6 at spacing=0) cannot land on it."""
        length = 60
        anchor = 30
        pattern = bytearray([3]) * length  # all 'T'
        pattern[anchor] = 0      # 'A' (ref base of the DBS)
        pattern[anchor + 1] = 1  # 'C' (2nd DBS base; the ONLY 'C' in the chromosome)
        chrom_dir = _write_chrom(tmp_path, bytes(pattern))

        mutation_tracker = {
            "DBS": {"S1": {"AC>GT": {"1": 1}}},
            "6": {"S1": {"C>A": {"1": 1}}},
        }

        rows = _run_sim(chrom_dir, tmp_path, mutation_tracker, ["DBS", "6"], spacing=0, label="dbs")
        dbs_rows = [r for r in rows if r[9] == "DBS"]
        snp_rows = [r for r in rows if r[9] == "SNP"]

        assert len(dbs_rows) == 1
        dbs_pos = int(dbs_rows[0][5])
        assert not snp_rows, (
            "an SBS6 event was placed on the DBS's reserved second base: "
            "DBS at %r, SBS6 rows %r" % (dbs_pos, snp_rows)
        )


# ─── ID chromosome allocation (audit sigprofilersimulator-2026-05 #10, #11) ───

def _write_fake_chrom(chrom_dir, name, length):
    (chrom_dir / (name + ".txt")).write_bytes(b"A" * length)


class TestIDChromosomeAllocation:
    def test_male_and_female_caches_do_not_collide(self, tmp_path):
        """Finding #10: a male run's cached proportions (with Y) and a female
        run's (without Y) must coexist, not overwrite each other, so a stale
        vector is never read misaligned against the wrong chromosome list."""
        chrom_dir = tmp_path / "chroms"
        chrom_dir.mkdir()
        chrom_dir_str = str(chrom_dir) + "/"
        _write_fake_chrom(chrom_dir, "X", 40)
        _write_fake_chrom(chrom_dir, "Y", 20)
        _write_fake_chrom(chrom_dir, "1", 60)

        def existing_proportion_files():
            return sorted(
                p.name for p in chrom_dir.iterdir() if p.name.endswith("_proportions.txt")
            )

        # A male run caches proportions for X, Y, 1.
        chrom_proportions(chrom_dir_str, "FAKE", ["X", "Y", "1"])
        male_files = existing_proportion_files()
        male_contents = {name: (chrom_dir / name).read_bytes() for name in male_files}

        # A later female run computes (and caches) its own, over X, 1 only.
        chrom_proportions(chrom_dir_str, "FAKE", ["X", "1"])

        # The male cache must still be exactly what it was: a separate file,
        # not overwritten by the female run.
        for name, content in male_contents.items():
            assert (chrom_dir / name).read_bytes() == content, (
                "male cache file %r changed after a female run -- the two "
                "chromosome sets collided on one shared filename" % name
            )
        assert existing_proportion_files() != male_files, (
            "no new cache file was created for the female chromosome set -- "
            "it reused (and could misalign with) the male one"
        )

        samples = {"ID": {"S1": {"1:Del:C:0": 10}}}
        log_file = str(tmp_path / "log.txt")
        tracker = mut_tracker(
            ["S1"], samples, "S1", {"ID": "unused"}, chrom_dir_str, "FAKE",
            ["X", "1"], False, log_file,
        )

        allocation = tracker["ID"]["S1"]["1:Del:C:0"]
        assert allocation == {"X": 4, "1": 6}, allocation

    def test_residual_allocation_is_opportunity_weighted(self, tmp_path):
        """Finding #11: the leftover-mutation residual must favor chromosomes
        by their genomic opportunity, not be uniform."""
        chrom_dir = tmp_path / "chroms2"
        chrom_dir.mkdir()
        chrom_dir_str = str(chrom_dir) + "/"
        _write_fake_chrom(chrom_dir, "A", 999)
        _write_fake_chrom(chrom_dir, "B", 1)
        chrom_proportions(chrom_dir_str, "FAKE2", ["A", "B"])

        samples = {"ID": {"S1": {"1:Del:C:0": 500}}}
        log_file = str(tmp_path / "log2.txt")

        b_got_residual = 0
        trials = 30
        for _ in range(trials):
            tracker = mut_tracker(
                ["S1"], samples, "S1", {"ID": "unused"}, chrom_dir_str, "FAKE2",
                ["A", "B"], False, log_file,
            )
            allocation = tracker["ID"]["S1"]["1:Del:C:0"]
            assert allocation["A"] + allocation["B"] == 500
            if allocation["B"] == 1:
                b_got_residual += 1

        # B has ~0.1% of genomic opportunity: under uniform allocation it would
        # win the single residual slot roughly half the time; under
        # opportunity-weighted allocation, essentially never.
        assert b_got_residual <= 2, (
            "B (0.1%% of opportunity) won the residual slot %d/%d times -- "
            "looks uniform, not opportunity-weighted" % (b_got_residual, trials)
        )
