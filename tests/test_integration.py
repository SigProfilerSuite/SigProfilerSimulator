"""
Integration tests for SigProfilerSimulator.
Require GRCh37 genome installed via SigProfilerMatrixGenerator.
Skipped automatically when the genome is not available.
"""

import io
import os
import shutil
import tempfile
from contextlib import redirect_stdout

import pytest

from SigProfilerMatrixGenerator.scripts import MutationMatrixGenerator as matRef
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
from SigProfilerSimulator.SigProfilerSimulator import SigProfilerSimulator, probability

PROJECT = "test_sps"
GENOME = "GRCh37"
CONTEXTS = ["96"]
SIMULATIONS = 1
N_INPUT_MUTATIONS = 12

TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
VCF_PATH = os.path.join(TESTS_DIR, "data", "GRCh37_test.vcf")
SEED_FILE = os.path.join(
    os.path.dirname(TESTS_DIR), "SigProfilerSimulator", "seeds.txt"
)
REPRODUCIBILITY_FIXTURE = os.path.join(
    TESTS_DIR, "fixtures", "reproducibility", "GRCh37_96_region1_seed_fixed.maf"
)
DBS_VCF_PATH = os.path.join(TESTS_DIR, "data", "GRCh37_dbs_test.vcf")
INDEL_VCF_PATH = os.path.join(TESTS_DIR, "data", "GRCh37_indel_test.vcf")
BED_PATH = os.path.join(TESTS_DIR, "data", "GRCh37_test.bed")

MAF_HEADER = "\t".join([
    "Hugo_symbol", "Entrez_gene_ID", "Center", "Genome", "Chrom",
    "Start_position", "End_position", "Strand", "Variant_Classification",
    "Variant_Type", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2",
    "dbSNP_RS", "dbSNP_Val_Status", "Tumor_Sample_Barcode", "matGenClass",
])


# ─── Genome availability ──────────────────────────────────────────────────────

def _grch37_available():
    try:
        path, _ = matRef.reference_paths("GRCh37")
        return os.path.exists(path) and len(os.listdir(path)) > 2
    except Exception:
        return False


grch37 = pytest.mark.skipif(
    not _grch37_available(),
    reason="GRCh37 genome not installed",
)


# ─── Fixture helpers ──────────────────────────────────────────────────────────

def _setup_project(vcf_source, project=PROJECT):
    """Create a temp project directory with input VCF; return project_path."""
    project_path = tempfile.mkdtemp() + "/"
    input_dir = project_path + "input/"
    os.makedirs(input_dir)
    shutil.copy(vcf_source, input_dir)
    return project_path


def _run_matgen_and_sps(project_path, vcf_basename, genome=GENOME,
                        contexts=CONTEXTS, simulations=SIMULATIONS, **sps_kwargs):
    """Run matGen then SPS; return the output path root for assertions."""
    matGen.SigProfilerMatrixGeneratorFunc(PROJECT, genome, project_path, plot=False)
    SigProfilerSimulator(PROJECT, project_path, genome, contexts,
                         simulations=simulations, **sps_kwargs)
    return project_path


# ─── Standard end-to-end simulation ──────────────────────────────────────────

@pytest.fixture(scope="module")
def simulation_output():
    """Runs a full simulation once and returns the output MAF path."""
    project_path = _setup_project(VCF_PATH)
    try:
        matGen.SigProfilerMatrixGeneratorFunc(PROJECT, GENOME, project_path, plot=False)
        SigProfilerSimulator(PROJECT, project_path, GENOME, CONTEXTS, simulations=SIMULATIONS)
        context_str = "_".join(CONTEXTS)
        maf = (
            f"{project_path}output/simulations/"
            f"{PROJECT}_simulations_{GENOME}_{context_str}/{SIMULATIONS}.maf"
        )
        yield maf
    finally:
        shutil.rmtree(project_path, ignore_errors=True)


@grch37
class TestEndToEndSimulation:
    def test_output_maf_created(self, simulation_output):
        assert os.path.exists(simulation_output)

    def test_output_maf_header(self, simulation_output):
        with open(simulation_output) as f:
            header = f.readline().strip()
        assert header == MAF_HEADER

    def test_output_maf_mutation_count(self, simulation_output):
        with open(simulation_output) as f:
            n = sum(1 for _ in f) - 1
        assert n == N_INPUT_MUTATIONS

    def test_output_maf_genome_field(self, simulation_output):
        with open(simulation_output) as f:
            f.readline()
            first_data = f.readline().strip().split("\t")
        assert first_data[3] == GENOME

    def test_output_maf_variant_type_snp(self, simulation_output):
        with open(simulation_output) as f:
            lines = f.readlines()[1:]
        variant_types = {line.strip().split("\t")[9] for line in lines if line.strip()}
        assert variant_types == {"SNP"}


# ─── probability() function ───────────────────────────────────────────────────

@grch37
class TestProbabilityFunction:
    """Lines 81-128: probability() function."""

    def test_no_chromosome_raises(self):
        with pytest.raises(ValueError):
            probability(genome=GENOME)

    def test_no_position_raises(self):
        with pytest.raises(ValueError):
            probability(genome=GENOME, chromosome="1")

    def test_no_mutation_raises(self):
        with pytest.raises(ValueError):
            probability(genome=GENOME, chromosome="1", position="15832678")

    def test_computes_and_prints_value(self):
        buf = io.StringIO()
        with redirect_stdout(buf):
            probability(
                chromosome="1",
                position="15832678",
                mutation="A[T>C]C",
                genome=GENOME,
            )
        assert "probabilt" in buf.getvalue().lower()

    def test_exome_mode(self):
        buf = io.StringIO()
        with redirect_stdout(buf):
            probability(
                chromosome="1",
                position="15832678",
                mutation="A[T>C]C",
                genome=GENOME,
                exome=True,
            )
        assert "probabilt" in buf.getvalue().lower()


# ─── Path variants ────────────────────────────────────────────────────────────
# Lines 141 (no trailing slash), 205/208 (log files overwritten), 377-378 (output exists)

@pytest.fixture(scope="module")
def sim_path_variants():
    """Run simulation twice in the same project (second run: log files exist + output exists)."""
    project_path = _setup_project(VCF_PATH)
    try:
        # First run — also exercises line 141 (no trailing slash on project_path)
        path_no_slash = project_path.rstrip("/")  # line 141 fires
        matGen.SigProfilerMatrixGeneratorFunc(PROJECT, GENOME, project_path, plot=False)
        SigProfilerSimulator(PROJECT, path_no_slash, GENOME, CONTEXTS, simulations=SIMULATIONS)

        # Second run in same dir — log files exist (205, 208), output path exists (377-378)
        SigProfilerSimulator(PROJECT, project_path, GENOME, CONTEXTS, simulations=SIMULATIONS)

        context_str = "_".join(CONTEXTS)
        maf = (
            f"{project_path}output/simulations/"
            f"{PROJECT}_simulations_{GENOME}_{context_str}/{SIMULATIONS}.maf"
        )
        yield maf
    finally:
        shutil.rmtree(project_path, ignore_errors=True)


@grch37
class TestSimulatorPathVariants:
    def test_second_run_produces_maf(self, sim_path_variants):
        assert os.path.exists(sim_path_variants)

    def test_second_run_has_correct_mutation_count(self, sim_path_variants):
        with open(sim_path_variants) as f:
            n = sum(1 for _ in f) - 1
        assert n == N_INPUT_MUTATIONS


# ─── region parameter ─────────────────────────────────────────────────────────
# Lines 194, 392-393

@pytest.fixture(scope="module")
def sim_region():
    project_path = _setup_project(VCF_PATH)
    try:
        matGen.SigProfilerMatrixGeneratorFunc(PROJECT, GENOME, project_path, plot=False)
        SigProfilerSimulator(PROJECT, project_path, GENOME, CONTEXTS,
                             simulations=SIMULATIONS, region="1")
        context_str = "_".join(CONTEXTS)
        maf = (
            f"{project_path}output/simulations/"
            f"{PROJECT}_simulations_{GENOME}_{context_str}/{SIMULATIONS}.maf"
        )
        yield maf
    finally:
        shutil.rmtree(project_path, ignore_errors=True)


@grch37
class TestSimulatorRegion:
    def test_output_maf_created(self, sim_region):
        assert os.path.exists(sim_region)

    def test_only_chr1_mutations(self, sim_region):
        with open(sim_region) as f:
            lines = f.readlines()[1:]
        chroms = {l.strip().split("\t")[4] for l in lines if l.strip()}
        assert chroms == {"1"}


# ─── vcf=True output ──────────────────────────────────────────────────────────
# Lines 401-405

@pytest.fixture(scope="module")
def sim_vcf_output():
    project_path = _setup_project(VCF_PATH)
    try:
        matGen.SigProfilerMatrixGeneratorFunc(PROJECT, GENOME, project_path, plot=False)
        SigProfilerSimulator(PROJECT, project_path, GENOME, CONTEXTS,
                             simulations=SIMULATIONS, vcf=True, region="1")
        yield project_path
    finally:
        shutil.rmtree(project_path, ignore_errors=True)


@grch37
class TestSimulatorVcfOutput:
    def test_sample_vcf_dirs_created(self, sim_vcf_output):
        context_str = "_".join(CONTEXTS)
        sim_dir = (
            f"{sim_vcf_output}output/simulations/"
            f"{PROJECT}_simulations_{GENOME}_{context_str}/"
        )
        assert os.path.isdir(sim_dir)


# ─── seqInfo=True output ──────────────────────────────────────────────────────
# Lines 454-464

@pytest.fixture(scope="module")
def sim_seqinfo():
    project_path = _setup_project(VCF_PATH)
    try:
        matGen.SigProfilerMatrixGeneratorFunc(PROJECT, GENOME, project_path, plot=False)
        SigProfilerSimulator(PROJECT, project_path, GENOME, CONTEXTS,
                             simulations=SIMULATIONS, seqInfo=True, region="1")
        yield project_path
    finally:
        shutil.rmtree(project_path, ignore_errors=True)


@grch37
class TestSimulatorSeqInfo:
    def test_seqout_dir_created(self, sim_seqinfo):
        seq_dir = sim_seqinfo + "output/vcf_files/simulations/"
        assert os.path.isdir(seq_dir)

    def test_context_subdir_created(self, sim_seqinfo):
        ctx_dir = sim_seqinfo + "output/vcf_files/simulations/96/"
        assert os.path.isdir(ctx_dir)


# ─── exome=True ───────────────────────────────────────────────────────────────
# Lines 271, 312, 333, 371, 450-451

@pytest.fixture(scope="module")
def sim_exome():
    project_path = _setup_project(VCF_PATH)
    try:
        matGen.SigProfilerMatrixGeneratorFunc(PROJECT, GENOME, project_path,
                                              plot=False, exome=True)
        SigProfilerSimulator(PROJECT, project_path, GENOME, CONTEXTS,
                             simulations=SIMULATIONS, exome=True)
        context_str = "_".join(CONTEXTS)
        maf = (
            f"{project_path}output/simulations/"
            f"{PROJECT}_simulations_{GENOME}_{context_str}_exome/{SIMULATIONS}.maf"
        )
        yield maf
    finally:
        shutil.rmtree(project_path, ignore_errors=True)


@grch37
class TestSimulatorExome:
    def test_exome_maf_created(self, sim_exome):
        assert os.path.exists(sim_exome)

    def test_exome_output_path_contains_exome(self, sim_exome):
        assert "_exome" in sim_exome


# ─── bed_file + region ────────────────────────────────────────────────────────
# Lines 149, 243-245, 274, 327-328, 369
# Note: region= is required to avoid deleting shared context-distribution files

@pytest.fixture(scope="module")
def sim_bed():
    project_path = _setup_project(VCF_PATH)
    try:
        matGen.SigProfilerMatrixGeneratorFunc(PROJECT, GENOME, project_path,
                                              plot=False, bed_file=BED_PATH)
        SigProfilerSimulator(PROJECT, project_path, GENOME, CONTEXTS,
                             simulations=SIMULATIONS, bed_file=BED_PATH, region="1")
        context_str = "_".join(CONTEXTS)
        maf = (
            f"{project_path}output/simulations/"
            f"{PROJECT}_simulations_{GENOME}_{context_str}_BED/{SIMULATIONS}.maf"
        )
        yield maf
    finally:
        shutil.rmtree(project_path, ignore_errors=True)


@grch37
class TestSimulatorBED:
    def test_bed_maf_created(self, sim_bed):
        assert os.path.exists(sim_bed)

    def test_bed_output_path_contains_BED(self, sim_bed):
        assert "_BED" in sim_bed


# ─── chrom_based=True ─────────────────────────────────────────────────────────
# Lines 286-298, 389-390

@pytest.fixture(scope="module")
def sim_chrom_based():
    project_path = _setup_project(VCF_PATH)
    try:
        matGen.SigProfilerMatrixGeneratorFunc(PROJECT, GENOME, project_path,
                                              plot=False, chrom_based=True)
        SigProfilerSimulator(PROJECT, project_path, GENOME, CONTEXTS,
                             simulations=SIMULATIONS, chrom_based=True, region="1")
        context_str = "_".join(CONTEXTS)
        maf = (
            f"{project_path}output/simulations/"
            f"{PROJECT}_simulations_{GENOME}_{context_str}/{SIMULATIONS}.maf"
        )
        yield maf
    finally:
        shutil.rmtree(project_path, ignore_errors=True)


@grch37
class TestSimulatorChromBased:
    def test_chrom_based_maf_created(self, sim_chrom_based):
        assert os.path.exists(sim_chrom_based)


# ─── SPS without prior matGen (auto-creates catalogue) ────────────────────────
# Lines 304-308

@pytest.fixture(scope="module")
def sim_no_prior_matgen():
    """Call SPS directly — it detects missing catalogue and runs matGen itself."""
    project_path = _setup_project(VCF_PATH)
    try:
        # Do NOT call matGen first: SPS will detect missing catalogue (lines 304-308)
        SigProfilerSimulator(PROJECT, project_path, GENOME, CONTEXTS,
                             simulations=SIMULATIONS, region="1")
        context_str = "_".join(CONTEXTS)
        maf = (
            f"{project_path}output/simulations/"
            f"{PROJECT}_simulations_{GENOME}_{context_str}/{SIMULATIONS}.maf"
        )
        yield maf
    finally:
        shutil.rmtree(project_path, ignore_errors=True)


@grch37
class TestSimulatorNoMatGen:
    def test_sps_creates_catalogue_and_runs(self, sim_no_prior_matgen):
        assert os.path.exists(sim_no_prior_matgen)


# ─── DBS context ──────────────────────────────────────────────────────────────
# Lines 252-257, 392-393

@pytest.fixture(scope="module")
def sim_dbs():
    project_path = _setup_project(DBS_VCF_PATH)
    try:
        matGen.SigProfilerMatrixGeneratorFunc(PROJECT, GENOME, project_path, plot=False)
        SigProfilerSimulator(PROJECT, project_path, GENOME, ["DBS"],
                             simulations=SIMULATIONS, region="1")
        maf = (
            f"{project_path}output/simulations/"
            f"{PROJECT}_simulations_{GENOME}_DBS/{SIMULATIONS}.maf"
        )
        yield maf
    finally:
        shutil.rmtree(project_path, ignore_errors=True)


@grch37
class TestSimulatorDBS:
    def test_dbs_maf_created(self, sim_dbs):
        assert os.path.exists(sim_dbs)

    def test_dbs_variant_type(self, sim_dbs):
        with open(sim_dbs) as f:
            lines = f.readlines()[1:]
        if lines:
            types = {l.strip().split("\t")[9] for l in lines if l.strip()}
            assert types <= {"SNP", "DBS", "DNP"}


# ─── INDEL context ────────────────────────────────────────────────────────────
# Lines 259-264, 392-393

@pytest.fixture(scope="module")
def sim_indel():
    project_path = _setup_project(INDEL_VCF_PATH)
    try:
        matGen.SigProfilerMatrixGeneratorFunc(PROJECT, GENOME, project_path, plot=False)
        SigProfilerSimulator(PROJECT, project_path, GENOME, ["ID"],
                             simulations=SIMULATIONS, region="1")
        maf = (
            f"{project_path}output/simulations/"
            f"{PROJECT}_simulations_{GENOME}_ID/{SIMULATIONS}.maf"
        )
        yield maf
    finally:
        shutil.rmtree(project_path, ignore_errors=True)


@grch37
class TestSimulatorINDEL:
    def test_indel_maf_created(self, sim_indel):
        assert os.path.exists(sim_indel)


# ─── Multiple simulations (simulations=2) ─────────────────────────────────────
# Covers iter_bin reset (line 394) when simulations > max_seed with region="1"

@pytest.fixture(scope="module")
def sim_multi():
    project_path = _setup_project(VCF_PATH)
    try:
        matGen.SigProfilerMatrixGeneratorFunc(PROJECT, GENOME, project_path, plot=False)
        SigProfilerSimulator(PROJECT, project_path, GENOME, CONTEXTS,
                             simulations=2, region="1", seed_file=SEED_FILE)
        context_str = "_".join(CONTEXTS)
        sim_dir = (
            f"{project_path}output/simulations/"
            f"{PROJECT}_simulations_{GENOME}_{context_str}/"
        )
        yield sim_dir
    finally:
        shutil.rmtree(project_path, ignore_errors=True)


@grch37
class TestSimulatorMultipleSimulations:
    def test_first_maf_created(self, sim_multi):
        assert os.path.exists(sim_multi + "1.maf")

    def test_second_maf_created(self, sim_multi):
        assert os.path.exists(sim_multi + "2.maf")

    def test_each_maf_has_correct_mutation_count(self, sim_multi):
        for fname in ("1.maf", "2.maf"):
            with open(sim_multi + fname) as f:
                n = sum(1 for _ in f) - 1
            assert n == N_INPUT_MUTATIONS


# ─── Male gender ──────────────────────────────────────────────────────────────
# gender='male' keeps Y chromosome; with many chromosomes may cover chrom_bin
# reset (line 386) depending on available CPUs

@pytest.fixture(scope="module")
def sim_male():
    project_path = _setup_project(VCF_PATH)
    try:
        matGen.SigProfilerMatrixGeneratorFunc(PROJECT, GENOME, project_path, plot=False)
        SigProfilerSimulator(PROJECT, project_path, GENOME, CONTEXTS,
                             simulations=SIMULATIONS, region="1", gender="male",
                             seed_file=SEED_FILE)
        context_str = "_".join(CONTEXTS)
        maf = (
            f"{project_path}output/simulations/"
            f"{PROJECT}_simulations_{GENOME}_{context_str}/{SIMULATIONS}.maf"
        )
        yield maf
    finally:
        shutil.rmtree(project_path, ignore_errors=True)


@grch37
class TestSimulatorMaleGender:
    def test_male_maf_created(self, sim_male):
        assert os.path.exists(sim_male)

    def test_male_maf_has_correct_mutation_count(self, sim_male):
        with open(sim_male) as f:
            n = sum(1 for _ in f) - 1
        assert n == N_INPUT_MUTATIONS


# ─── Noise: Poisson and Uniform ───────────────────────────────────────────────

@pytest.fixture(scope="module")
def sim_noise_poisson():
    project_path = _setup_project(VCF_PATH)
    try:
        matGen.SigProfilerMatrixGeneratorFunc(PROJECT, GENOME, project_path, plot=False)
        SigProfilerSimulator(PROJECT, project_path, GENOME, CONTEXTS,
                             simulations=SIMULATIONS, region="1",
                             noisePoisson=True, seed_file=SEED_FILE)
        context_str = "_".join(CONTEXTS)
        maf = (
            f"{project_path}output/simulations/"
            f"{PROJECT}_simulations_{GENOME}_{context_str}/{SIMULATIONS}.maf"
        )
        yield maf
    finally:
        shutil.rmtree(project_path, ignore_errors=True)


@pytest.fixture(scope="module")
def sim_noise_uniform():
    project_path = _setup_project(VCF_PATH)
    try:
        matGen.SigProfilerMatrixGeneratorFunc(PROJECT, GENOME, project_path, plot=False)
        SigProfilerSimulator(PROJECT, project_path, GENOME, CONTEXTS,
                             simulations=SIMULATIONS, region="1",
                             noiseUniform=20, seed_file=SEED_FILE)
        context_str = "_".join(CONTEXTS)
        maf = (
            f"{project_path}output/simulations/"
            f"{PROJECT}_simulations_{GENOME}_{context_str}/{SIMULATIONS}.maf"
        )
        yield maf
    finally:
        shutil.rmtree(project_path, ignore_errors=True)


@grch37
class TestSimulatorNoise:
    def test_poisson_noise_maf_created(self, sim_noise_poisson):
        assert os.path.exists(sim_noise_poisson)

    def test_poisson_noise_maf_has_mutations(self, sim_noise_poisson):
        with open(sim_noise_poisson) as f:
            n = sum(1 for _ in f) - 1
        assert n > 0

    def test_uniform_noise_maf_created(self, sim_noise_uniform):
        assert os.path.exists(sim_noise_uniform)

    def test_uniform_noise_maf_has_mutations(self, sim_noise_uniform):
        with open(sim_noise_uniform) as f:
            n = sum(1 for _ in f) - 1
        assert n > 0


# ─── seqInfo second run (lines 430-432) ───────────────────────────────────────
# Second run with seqInfo=True in same dir: context subdir already exists,
# triggers rmtree + makedirs branch

@pytest.fixture(scope="module")
def sim_seqinfo_twice():
    project_path = _setup_project(VCF_PATH)
    try:
        matGen.SigProfilerMatrixGeneratorFunc(PROJECT, GENOME, project_path, plot=False)
        SigProfilerSimulator(PROJECT, project_path, GENOME, CONTEXTS,
                             simulations=SIMULATIONS, seqInfo=True, region="1",
                             seed_file=SEED_FILE)
        SigProfilerSimulator(PROJECT, project_path, GENOME, CONTEXTS,
                             simulations=SIMULATIONS, seqInfo=True, region="1",
                             seed_file=SEED_FILE)
        yield project_path
    finally:
        shutil.rmtree(project_path, ignore_errors=True)


@grch37
class TestSimulatorSeqInfoSecondRun:
    def test_seqinfo_second_run_context_dir_exists(self, sim_seqinfo_twice):
        ctx_dir = sim_seqinfo_twice + "output/vcf_files/simulations/96/"
        assert os.path.isdir(ctx_dir)


# ─── Non-standard DBS context (line 224) ──────────────────────────────────────
# context "DBS186" → file_name = '.' + context (not '.DBS78')

@pytest.fixture(scope="module")
def sim_dbs186():
    project_path = _setup_project(DBS_VCF_PATH)
    try:
        matGen.SigProfilerMatrixGeneratorFunc(PROJECT, GENOME, project_path, plot=False)
        SigProfilerSimulator(PROJECT, project_path, GENOME, ["DBS186"],
                             simulations=SIMULATIONS, region="1", seed_file=SEED_FILE)
        maf = (
            f"{project_path}output/simulations/"
            f"{PROJECT}_simulations_{GENOME}_DBS186/{SIMULATIONS}.maf"
        )
        yield maf
    finally:
        shutil.rmtree(project_path, ignore_errors=True)


@grch37
class TestSimulatorDBS186:
    def test_dbs186_maf_created(self, sim_dbs186):
        assert os.path.exists(sim_dbs186)


# ─── Non-standard ID context (line 231) ───────────────────────────────────────
# context "ID415" → file_name = '.' + context (not '.ID83')

@pytest.fixture(scope="module")
def sim_id415():
    project_path = _setup_project(INDEL_VCF_PATH)
    try:
        matGen.SigProfilerMatrixGeneratorFunc(PROJECT, GENOME, project_path, plot=False)
        SigProfilerSimulator(PROJECT, project_path, GENOME, ["ID415"],
                             simulations=SIMULATIONS, region="1", seed_file=SEED_FILE)
        maf = (
            f"{project_path}output/simulations/"
            f"{PROJECT}_simulations_{GENOME}_ID415/{SIMULATIONS}.maf"
        )
        yield maf
    finally:
        shutil.rmtree(project_path, ignore_errors=True)


@grch37
class TestSimulatorID415:
    def test_id415_maf_created(self, sim_id415):
        assert os.path.exists(sim_id415)


# ─── Reproducibility ──────────────────────────────────────────────────────────

@pytest.fixture(scope="module")
def sim_reproducibility():
    """Run two simulations with the same seed_file; return both MAF paths."""
    results = []
    for _ in range(2):
        project_path = _setup_project(VCF_PATH)
        try:
            matGen.SigProfilerMatrixGeneratorFunc(PROJECT, GENOME, project_path, plot=False)
            SigProfilerSimulator(PROJECT, project_path, GENOME, CONTEXTS,
                                 simulations=SIMULATIONS, region="1", seed_file=SEED_FILE)
            maf = (
                f"{project_path}output/simulations/"
                f"{PROJECT}_simulations_{GENOME}_96/{SIMULATIONS}.maf"
            )
            results.append((project_path, maf))
        except Exception:
            shutil.rmtree(project_path, ignore_errors=True)
            raise
    yield results[0][1], results[1][1]
    for project_path, _ in results:
        shutil.rmtree(project_path, ignore_errors=True)


@grch37
class TestReproducibility:
    def test_two_runs_same_seed_identical(self, sim_reproducibility):
        maf1, maf2 = sim_reproducibility
        with open(maf1) as f1, open(maf2) as f2:
            assert f1.read() == f2.read()

    def test_output_matches_stored_fixture(self, sim_reproducibility):
        maf1, _ = sim_reproducibility
        with open(maf1) as f, open(REPRODUCIBILITY_FIXTURE) as ref:
            assert f.read() == ref.read()
